#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  Policy Control Barrier Function (PCBF) — Circle Obstacle Avoidance
//
//  Based on:
//    Knoedler et al., "Safety on the Fly: Constructing Robust Safety Filters
//    via Policy Control Barrier Functions at Runtime", IEEE RA-L 2025.
//
//  HOW THIS DIFFERS FROM THE STANDARD CBF (cbf_circle_obstacles_qp_lim.c):
//
//  Standard CBF:
//    Uses h_i(x,y) directly as the barrier value and ∇h_i as the gradient.
//    This is instantaneous — it only looks at the current state.
//
//  PCBF (this file):
//    Instead of h_i(x,y), it computes V_T^{h,π} — the MAXIMUM value of h_i
//    along a finite-horizon rollout of a design policy π. This looks ahead in
//    time, so the filter activates earlier and handles high relative-degree
//    systems more gracefully.
//
//  THE MATH (from the paper, Section III.B):
//
//    Barrier function:
//      h_i(x,y) = (x-cx_i)^2 + (y-cy_i)^2 - r_i^2
//      Positive = safe (outside), negative = violation (inside).
//
//    Policy value function (finite horizon approximation, Eq. 10):
//      V_T^{h,π}(x0) = max_{0 <= k < H}  h_i(x_k)
//      where x_k is the state at step k when following policy π from x0.
//      This is the worst (maximum) constraint violation over the horizon.
//
//    CBF safety condition (Eq. 4b of paper):
//      ∇V^T · (f(x) + g(x)·u) + α·V ≤ 0
//      Rearranged into QP inequality:  A_i · u ≤ b_i
//
//    Gradient of V (finite-horizon, via chain rule through rollout):
//      ∇V = ∇h(x_k*) · ∂x_k* / ∂x0
//      where k* = argmax h_i(x_k) is the step where max occurs.
//      ∂x_k / ∂x0 is the state transition sensitivity (propagated forward).
//
//    QP solved each control step (Eq. CBF-QP of paper):
//      min   ½ (u - u₀)ᵀ W (u - u₀)        [stay close to nominal]
//      s.t.  A · u ≤ b                       [PCBF safety, one row per obstacle]
//            u_min ≤ u ≤ u_max              [input limits]
//
// =============================================================================


// ---- Obstacle definition ----------------------------------------------------
typedef struct {
    double cx;     // center x (world frame)
    double cy;     // center y (world frame)
    double r;      // safety radius
    double alpha;  // CBF gain α (larger = more aggressive near boundary)
    int    active; // 1 = on, 0 = off
} CircleCBF;

// ---- Edit obstacles here ----------------------------------------------------
#define N_OBS 4

const double alpha = 5.0;
const double buffer = 0.1;

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, alpha, 1 },
    { -4.0, -1.0, 0.35, alpha, 1 },
    { -2.0, -1.0, 0.35, alpha, 1 },
    { -4.0,  1.0, 0.35, alpha, 1 },
};
// -----------------------------------------------------------------------------


// ---- Rollout parameters -----------------------------------------------------
#define PCBF_H      50      // number of rollout steps (horizon length)
#define PCBF_DT     0.1     // rollout timestep (seconds)

// Design policy selector — change this single value to switch policies:
//   POLICY_ZERO_VELOCITY : π = [0, 0, 0]  (robot brakes to a stop)
//   POLICY_NOMINAL       : π = current nominal command held constant
#define POLICY_ZERO_VELOCITY  0
#define POLICY_NOMINAL        1
static int pcbf_policy = POLICY_ZERO_VELOCITY;   // <-- change here to switch
// -----------------------------------------------------------------------------


// ---- Input limits (keep in sync with control_params.h) ---------------------
static double u_max[3] = {  1.0,  0.3,  1.0 };
static double u_min[3] = { -1.0, -0.3, -1.0 };

// ---- Slack variable penalty -------------------------------------------------
// The CBF constraint is made SOFT via a slack variable δ >= 0:
//   A·u <= b + δ     (CBF, softened)
//   u_min <= u <= u_max  (input limits, still hard)
// The QP cost gains a term: p_slack * δ²
// Large p_slack = CBF is almost hard (penalize violations heavily).
// Small p_slack = more permissive but always feasible.
// This prevents primal infeasibility when input limits conflict with CBF.
static double p_slack = 1e4;
// -----------------------------------------------------------------------------


// ---- OSQP persistent state --------------------------------------------------
static OSQPSolver    *pcbf_solver   = NULL;
static OSQPSettings   pcbf_settings;

// Variables: u = [vx, vy, wz, δ]  (4 variables, δ = slack)
// Constraints:
//   Rows 0..nc-1   : A_cbf·[vx,vy,wz] - δ <= b_cbf   (soft CBF)
//   Rows nc..nc+2  : identity on [vx,vy,wz]           (input limits)
//   Row  nc+3      : δ >= 0  i.e.  -δ <= 0            (slack non-negative)
#define N_VARS  4               // [vx, vy, wz, delta]
#define N_ROWS  (N_OBS + 3 + 1) // CBF rows + box rows + slack >= 0 row

static OSQPFloat pcbf_P_x[N_VARS];          // diagonal of P (4 entries)
static OSQPInt   pcbf_P_i[N_VARS];
static OSQPInt   pcbf_P_p[N_VARS+1];
static OSQPCscMatrix pcbf_P_mat;

static OSQPFloat pcbf_A_x[N_ROWS * N_VARS];
static OSQPInt   pcbf_A_i[N_ROWS * N_VARS];
static OSQPInt   pcbf_A_p[N_VARS+1];
static OSQPCscMatrix pcbf_A_mat;

static OSQPFloat pcbf_q[N_VARS];
static OSQPFloat pcbf_l[N_ROWS];
static OSQPFloat pcbf_u[N_ROWS];
// -----------------------------------------------------------------------------


// ---- Unicycle model rollout -------------------------------------------------
// The high-level controller uses a unicycle (kinematic) model:
//   ẋ = cos(θ)·vx - sin(θ)·vy
//   ẏ = sin(θ)·vx + cos(θ)·vy
//   θ̇ = wz
//
// We roll this forward H steps using Euler integration at timestep PCBF_DT.
// State: [x, y, theta].
// Input: u = [vx, vy, wz] — held constant across the horizon (design policy).
//
// Sensitivity ∂x_k/∂x0 is a 3x3 matrix propagated via linearized dynamics.
// This lets us compute the gradient of V w.r.t. x0 via the chain rule.
// -----------------------------------------------------------------------------


// =============================================================================
//  pcbf_rollout_and_value()
//
//  Rolls out the design policy from (x0, y0, theta0) for H steps.
//  For obstacle i, computes:
//    V_i   = max_{k=0..H-1} h_i(x_k, y_k)          (Eq. 10 of paper)
//    dV_dx = gradient of V_i w.r.t. [x0, y0, theta0] (for CBF-QP row)
//
//  Outputs:
//    V_val[i]   : scalar barrier value for obstacle i
//    dV_dx0[i]  : 3-vector gradient [∂V/∂x, ∂V/∂y, ∂V/∂θ]
// =============================================================================
static void pcbf_rollout_and_value(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS],
    double dV_dx0[N_OBS][3])
{
    // Policy input held constant over horizon
    double ux, uy, uw;
    if (pcbf_policy == POLICY_ZERO_VELOCITY) {
        ux = 0.0;  uy = 0.0;  uw = 0.0;
    } else {
        ux = vx_nom;  uy = vy_nom;  uw = wz_nom;
    }

    // Initialize obstacle tracking — also evaluate h at initial state (k=0)
    double h_max[N_OBS];
    int    k_max[N_OBS];
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { h_max[i] = -1e9; k_max[i] = 0; continue; }
        double dx0 = x0 - obstacles[i].cx;
        double dy0 = y0 - obstacles[i].cy;
        h_max[i] = obstacles[i].r * obstacles[i].r - dx0*dx0 - dy0*dy0;
        k_max[i] = 0;
    }

    // State trajectory storage (need positions to compute gradient later)
    double xs[PCBF_H+1], ys[PCBF_H+1], ths[PCBF_H+1];
    xs[0]  = x0;
    ys[0]  = y0;
    ths[0] = theta0;

    // -------------------------------------------------------------------------
    // Forward rollout: integrate unicycle dynamics H steps
    // -------------------------------------------------------------------------
    for (int k = 0; k < PCBF_H; k++)
    {
        double c = cos(ths[k]);
        double s = sin(ths[k]);

        // Unicycle Euler step
        xs[k+1]  = xs[k]  + PCBF_DT * (c*ux - s*uy);
        ys[k+1]  = ys[k]  + PCBF_DT * (s*ux + c*uy);
        ths[k+1] = ths[k] + PCBF_DT * uw;

        // h_i using PAPER convention (Section III.A, Eq. 5):
        //   A = {x | h(x) > 0}  →  h > 0 means INSIDE obstacle (unsafe)
        //                           h < 0 means OUTSIDE obstacle (safe)
        //   h_i = r² - dx² - dy²
        //
        // This is the NEGATION of the standard CBF used in
        // cbf_circle_obstacles_qp_lim.c where h = dx²+dy²-r² (h>0 outside).
        // The sign flip propagates to V and ∇V, but the QP form is consistent.
        //
        // With this convention:
        //   V ≤ 0  → safe (max h along rollout is non-positive = stayed outside)
        //   V > 0  → unsafe (worst-case trajectory enters obstacle)
        //
        // Far from obstacle: h = r²-(large) << 0 → V << 0 → b = -α·V >> 0
        //   → constraint A·u ≤ large_positive → easy to satisfy → no correction
        // Near boundary:     h → 0 → V → 0 → b → 0 → constraint tightens
        // Inside obstacle:   h > 0 → V > 0 → b < 0 → correction forced
        for (int i = 0; i < N_OBS; i++) {
            if (!obstacles[i].active) continue;
            double dx = xs[k+1] - obstacles[i].cx;
            double dy = ys[k+1] - obstacles[i].cy;
            double h  = obstacles[i].r * obstacles[i].r - dx*dx - dy*dy;
            if (h > h_max[i]) {
                h_max[i] = h;
                k_max[i] = k + 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Gradient computation via sensitivity propagation
    //
    // We need: ∇_x0 V = ∇h(x_{k*}) · ∂x_{k*}/∂x0
    //
    // The sensitivity matrix S_k = ∂x_k/∂x0 evolves as:
    //   S_{k+1} = A_k · S_k
    // where A_k = ∂f/∂x (Jacobian of unicycle dynamics at step k):
    //
    //   A_k = I + dt * [  0   0   (-s*ux - c*uy) ]
    //                  [  0   0   ( c*ux - s*uy) ]
    //                  [  0   0        0          ]
    //
    // ∇h_i at state (x_k*, y_k*):
    //   ∇h = [2*(x_k* - cx_i),  2*(y_k* - cy_i),  0]
    //
    // Then: dV/dx0 = ∇h · S_{k*}
    // -------------------------------------------------------------------------

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) {
            V_val[i] = -1e9;  // safely inactive: V << 0, no constraint generated
            dV_dx0[i][0] = 0; dV_dx0[i][1] = 0; dV_dx0[i][2] = 0;
            continue;
        }

        V_val[i] = h_max[i];

        // Propagate sensitivity S from step 0 to k*
        // S starts as identity (∂x0/∂x0 = I)
        double S[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

        for (int k = 0; k < k_max[i]; k++)
        {
            double c = cos(ths[k]);
            double s = sin(ths[k]);

            // dtheta terms: partial derivatives of ẋ, ẏ w.r.t. theta
            double df_dtheta_x = -s*ux - c*uy;
            double df_dtheta_y =  c*ux - s*uy;

            // A_k = I + dt * Jacobian
            // Only the third column (theta sensitivity) is non-trivial
            double A[3][3] = {
                {1, 0, PCBF_DT * df_dtheta_x},
                {0, 1, PCBF_DT * df_dtheta_y},
                {0, 0, 1}
            };

            // S_new = A * S
            double S_new[3][3] = {{0}};
            for (int r = 0; r < 3; r++)
                for (int c2 = 0; c2 < 3; c2++)
                    for (int m = 0; m < 3; m++)
                        S_new[r][c2] += A[r][m] * S[m][c2];

            memcpy(S, S_new, sizeof(S));
        }

        // ∇h at the max step using paper convention h = r² - dx² - dy²:
        //   ∂h/∂x = -2*(x - cx) = -2*dx
        //   ∂h/∂y = -2*(y - cy) = -2*dy
        //   ∂h/∂θ = 0
        double dx_star = xs[k_max[i]] - obstacles[i].cx;
        double dy_star = ys[k_max[i]] - obstacles[i].cy;
        double grad_h[3] = {-2.0*dx_star, -2.0*dy_star, 0.0};

        // dV/dx0 = grad_h · S
        for (int j = 0; j < 3; j++) {
            dV_dx0[i][j] = 0.0;
            for (int m = 0; m < 3; m++)
                dV_dx0[i][j] += grad_h[m] * S[m][j];
        }
    }
}


// ---- OSQP setup (called once on first use) ----------------------------------
static int pcbf_init_solver(double wx, double wy, double ww, int nc)
{
    int nr = nc + 3 + 1;  // CBF rows + box rows + slack >= 0 row

    // ---- P matrix: 4x4 diagonal [wx, wy, ww, p_slack] ----------------------
    // Stored as upper-triangular CSC (diagonal only for diagonal matrix)
    pcbf_P_x[0] = (OSQPFloat)wx;
    pcbf_P_x[1] = (OSQPFloat)wy;
    pcbf_P_x[2] = (OSQPFloat)ww;
    pcbf_P_x[3] = (OSQPFloat)p_slack;
    for (int j = 0; j < N_VARS; j++) { pcbf_P_i[j] = j; pcbf_P_p[j] = j; }
    pcbf_P_p[N_VARS] = N_VARS;

    pcbf_P_mat.m = N_VARS; pcbf_P_mat.n = N_VARS;
    pcbf_P_mat.nz = -1;    pcbf_P_mat.nzmax = N_VARS;
    pcbf_P_mat.x = pcbf_P_x; pcbf_P_mat.i = pcbf_P_i; pcbf_P_mat.p = pcbf_P_p;

    // ---- A matrix: nr x 4, column-major CSC ---------------------------------
    //
    // Column layout (one column per variable):
    //
    //  col 0 (vx):  nc CBF rows (A_cbf[row][0]),  1 identity row (box vx)
    //  col 1 (vy):  nc CBF rows (A_cbf[row][1]),  1 identity row (box vy)
    //  col 2 (wz):  nc CBF rows (A_cbf[row][2]),  1 identity row (box wz)
    //  col 3 (δ):   nc CBF rows (-1 each),         1 row for δ >= 0 (-1)
    //
    // Row layout:
    //   rows 0..nc-1   : CBF constraints (A_cbf·u - δ <= b_cbf)
    //   rows nc..nc+2  : box constraints on vx, vy, wz (identity)
    //   row  nc+3      : δ >= 0  stored as  -δ <= 0

    int nnz = 0;

    // col 0: vx
    pcbf_A_p[0] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;  // CBF row (updated each step)
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 0; nnz++;   // box: vx

    // col 1: vy
    pcbf_A_p[1] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 1; nnz++;   // box: vy

    // col 2: wz
    pcbf_A_p[2] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 2; nnz++;   // box: wz

    // col 3: δ (slack)
    pcbf_A_p[3] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = -1.0f; pcbf_A_i[nnz] = row; nnz++; // -δ in CBF rows
    }
    pcbf_A_x[nnz] = -1.0f; pcbf_A_i[nnz] = nc + 3; nnz++;  // -δ <= 0 (δ >= 0)

    pcbf_A_p[4] = nnz;

    pcbf_A_mat.m = nr; pcbf_A_mat.n = N_VARS;
    pcbf_A_mat.nz = -1; pcbf_A_mat.nzmax = nnz;
    pcbf_A_mat.x = pcbf_A_x; pcbf_A_mat.i = pcbf_A_i; pcbf_A_mat.p = pcbf_A_p;

    // ---- Initial bounds -----------------------------------------------------
    for (int j = 0; j < N_VARS; j++) pcbf_q[j] = 0.0f;

    // CBF rows: l = -inf, u = b_cbf (updated each step)
    for (int i = 0; i < nc; i++) {
        pcbf_l[i] = (OSQPFloat)(-OSQP_INFTY);
        pcbf_u[i] = 0.0f;
    }
    // Box rows: u_min <= u_j <= u_max (fixed)
    for (int j = 0; j < 3; j++) {
        pcbf_l[nc + j] = (OSQPFloat)u_min[j];
        pcbf_u[nc + j] = (OSQPFloat)u_max[j];
    }
    // Slack row: -inf <= -δ <= 0  (i.e. δ >= 0, δ uncapped above)
    pcbf_l[nc + 3] = (OSQPFloat)(-OSQP_INFTY);
    pcbf_u[nc + 3] = 0.0f;

    osqp_set_default_settings(&pcbf_settings);
    pcbf_settings.verbose       = 0;
    pcbf_settings.warm_starting = 1;
    pcbf_settings.eps_abs       = 1e-3;
    pcbf_settings.eps_rel       = 1e-3;
    pcbf_settings.max_iter      = 1000;

    OSQPInt flag = osqp_setup(&pcbf_solver, &pcbf_P_mat, pcbf_q,
                               &pcbf_A_mat, pcbf_l, pcbf_u, nr, N_VARS, &pcbf_settings);
    if (flag != 0) {
        fprintf(stderr, "[PCBF] OSQP setup failed (code %d)\n", (int)flag);
        return 0;
    }
    return 1;
}


// =============================================================================
//  cbf_circle_obstacles_filter()   — drop-in replacement for the standard CBF
//
//  Same function signature as cbf_circle_obstacles_qp_lim.c so you only need
//  to change the #include line in my_controller.c.
//
//  Inputs:
//    x, y, theta         — robot world position and yaw
//    vx_cmd, vy_cmd, wz_cmd — nominal commands (body frame), modified in-place
//    wx, wy, ww          — QP cost weights
//
//  Returns 1 if commands were modified, 0 otherwise.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    // =========================================================================
    // STEP 1: Roll out the design policy and compute V_T^{h,π} for each
    //         obstacle (Eq. 10 of paper), plus its gradient w.r.t. [x, y, θ].
    //
    //  V_val[i]   = max_{k=0..H-1} h_i(x_k)       (scalar, the barrier value)
    //  dV_dx0[i]  = ∇_{x0} V_val[i]               (3-vector, for QP row)
    // =========================================================================

    double V_val[N_OBS];
    double dV_dx0[N_OBS][3];   // [∂V/∂x,  ∂V/∂y,  ∂V/∂θ]

    pcbf_rollout_and_value(x, y, theta,
                           *vx_cmd, *vy_cmd, *wz_cmd,
                           V_val, dV_dx0);

    // =========================================================================
    // STEP 2: Build QP constraint rows from the PCBF values and gradients.
    //
    //  Standard CBF row (instantaneous):
    //    A_i = -∇h_i · g(x)           b_i = α · h_i + ∇h_i · f(x)
    //
    //  PCBF row (finite-horizon, from Eq. 4b + chain rule):
    //    The gradient of V already encodes how the worst-case h changes
    //    w.r.t. the initial state. The CBF-QP row becomes:
    //
    //    dV/dx · ẋ + α · V ≤ 0
    //
    //  For the unicycle:  ẋ = [cos θ · vx - sin θ · vy,
    //                          sin θ · vx + cos θ · vy,
    //                          wz]
    //
    //  Expanding and isolating the control input u = [vx, vy, wz]:
    //    (∂V/∂x)(c·vx - s·vy) + (∂V/∂y)(s·vx + c·vy) + (∂V/∂θ)·wz ≤ -α·V
    //
    //  So:
    //    A_i[0] = (∂V/∂x)·c + (∂V/∂y)·s       coefficient of vx
    //    A_i[1] = -(∂V/∂x)·s + (∂V/∂y)·c      coefficient of vy
    //    A_i[2] = ∂V/∂θ                         coefficient of wz
    //    b_i    = -α · V                         (note negative sign)
    // =========================================================================

    const double c = cos(theta);
    const double s = sin(theta);

    double A_cbf[N_OBS][3];
    double b_cbf[N_OBS];
    int nc = 0;

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) continue;

        double dVdx = dV_dx0[i][0];
        double dVdy = dV_dx0[i][1];
        double dVdt = dV_dx0[i][2];

        A_cbf[nc][0] =  dVdx*c + dVdy*s;   // coefficient of vx
        A_cbf[nc][1] = -dVdx*s + dVdy*c;   // coefficient of vy
        A_cbf[nc][2] =  dVdt;               // coefficient of wz
        b_cbf[nc]    = -obstacles[i].alpha * V_val[i] - buffer;

        nc++;
    }

    if (nc == 0) return 0;

    // =========================================================================
    // STEP 3: Initialize OSQP on first call
    // =========================================================================

    if (pcbf_solver == NULL) {
        if (!pcbf_init_solver(wx, wy, ww, nc)) return 0;
    }

    // =========================================================================
    // STEP 4: Update OSQP data for this timestep
    //
    //  q = [-wx·vx0, -wy·vy0, -ww·wz0, 0]   (δ has no linear cost term)
    //  A CBF rows: columns 0-2 get A_cbf values, column 3 stays -1 (fixed)
    //  Upper bounds: b_cbf for CBF rows, u_max for box rows (box unchanged)
    // =========================================================================

    pcbf_q[0] = (OSQPFloat)(-wx * (*vx_cmd));
    pcbf_q[1] = (OSQPFloat)(-wy * (*vy_cmd));
    pcbf_q[2] = (OSQPFloat)(-ww * (*wz_cmd));
    pcbf_q[3] = 0.0f;   // no linear cost on slack δ

    // Update A matrix — only the CBF rows in columns 0,1,2 change each step.
    // Column 3 (slack entries: -1 for CBF rows, -1 for slack>=0 row) never change.
    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < nc; row++)
            pcbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        nnz++;  // skip the identity entry for box row (stays 1.0, fixed)
    }
    // column 3 (slack) entries never change — skip entirely in update

    for (int i = 0; i < nc; i++)
        pcbf_u[i] = (OSQPFloat)b_cbf[i];

    osqp_update_data_vec(pcbf_solver, pcbf_q, pcbf_l, pcbf_u);
    osqp_update_data_mat(pcbf_solver, NULL, NULL, 0, pcbf_A_x, NULL, pcbf_A_mat.nzmax);

    // =========================================================================
    // STEP 5: Solve and extract result
    //
    //  Solution is [vx*, vy*, wz*, δ*].
    //  δ* > 0 means the CBF was softened — the constraint was infeasible
    //  with hard limits and a small safety violation was allowed.
    //  Print a warning in that case so you know it happened.
    // =========================================================================

    osqp_solve(pcbf_solver);

    int modified = 0;

    if (pcbf_solver->info->status_val == OSQP_SOLVED ||
        pcbf_solver->info->status_val == OSQP_SOLVED_INACCURATE)
    {
        *vx_cmd = (double)pcbf_solver->solution->x[0];
        *vy_cmd = (double)pcbf_solver->solution->x[1];
        *wz_cmd = (double)pcbf_solver->solution->x[2];
        double delta = (double)pcbf_solver->solution->x[3];
        if (delta > 1e-3)
            fprintf(stderr, "[PCBF] slack active: delta=%.4f (CBF softened)\n", delta);
        modified = 1;
    }
    else
    {
        fprintf(stderr, "[PCBF] OSQP did not converge (%s), clamping to limits\n",
                pcbf_solver->info->status);
        *vx_cmd = fmax(u_min[0], fmin(u_max[0], *vx_cmd));
        *vy_cmd = fmax(u_min[1], fmin(u_max[1], *vy_cmd));
        *wz_cmd = fmax(u_min[2], fmin(u_max[2], *wz_cmd));
    }

    return modified;
}
