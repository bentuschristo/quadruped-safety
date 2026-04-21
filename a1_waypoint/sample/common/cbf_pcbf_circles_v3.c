#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  Policy Control Barrier Function (PCBF) — Circle Obstacle Avoidance  (v3)
//
//  Changes from v2:
//    - Added POLICY_PUSH_VX: open-loop constant forward velocity design policy.
//      Rolls out the robot moving at full vx in its current heading direction,
//      regardless of the nominal command. This predicts danger from obstacles
//      directly ahead and produces smoother corrections than POLICY_NOMINAL
//      because the rollout is independent of the time-varying nominal command.
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
//      h_i(x,y) = r_i^2 - (x-cx_i)^2 - (y-cy_i)^2  (positive inside obstacle)
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

// Logger
static FILE* pcbf_log_file = NULL;

static void pcbf_open_log(void)
{
    if (pcbf_log_file != NULL) return;  // already open
    pcbf_log_file = fopen("pcbf_log.csv", "w");
    if (pcbf_log_file == NULL) {
        fprintf(stderr, "[PCBF] WARNING: could not open pcbf_log.csv\n");
        return;
    }
    // CSV header
    fprintf(pcbf_log_file,
        "step,"
        "x,y,theta,"
        "h0,h1,h2,h3,"
        "V0,V1,V2,V3,"
        "dVdx0,dVdy0,dVdt0,"
        "dVdx1,dVdy1,dVdt1,"
        "dVdx2,dVdy2,dVdt2,"
        "dVdx3,dVdy3,dVdt3\n");
}






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

const double buffer = 0.00;
const double alpha = 1.0;

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, alpha, 1 },
    { -4.0, -1.0, 0.35, alpha, 1 },
    { -2.0, -1.0, 0.35, alpha, 1 },
    { -4.0,  1.0, 0.35, alpha, 1 },
};
// -----------------------------------------------------------------------------


// ---- Rollout parameters -----------------------------------------------------
#define PCBF_H      20      // number of rollout steps (horizon length)
#define PCBF_DT     0.1     // rollout timestep (seconds)

// Design policy selector — change this single value to switch policies:
//   POLICY_ZERO_VELOCITY : π = [0,    0,  0]   robot stays stationary
//   POLICY_NOMINAL       : π = [vx_nom, vy_nom, wz_nom]  held constant
//   POLICY_PUSH_VX       : π = [vx_push, 0, 0]  constant forward push
//                          Open-loop: always moves at vx_push in current
//                          heading, ignoring the nominal command entirely.
//                          Predicts danger ahead without the discontinuities
//                          caused by switching nominal commands.
#define POLICY_ZERO_VELOCITY  0
#define POLICY_NOMINAL        1
#define POLICY_PUSH_VX        2

static int    pcbf_policy = POLICY_PUSH_VX;   // <-- change here to switch
static double vx_push     = 0.5;              // forward push speed (m/s)
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
//  CUBIC SPLINE HELPERS  (Section III.C of paper)
//
//  Instead of taking max over discrete h values (naive, causes gradient errors),
//  we fit a natural cubic spline through {h_k} and find the continuous maximum
//  analytically. This resolves the piecewise-constant gradient discontinuity
//  shown in Fig. 2 of the paper.
//
//  Natural cubic spline on n+1 uniform points (spacing dt):
//    Second derivatives M_k satisfy a tridiagonal system (Thomas algorithm).
//    On each segment [k, k+1] with τ = (t - t_k)/dt ∈ [0,1]:
//      s_k(τ) = a_k + b_k·τ + c_k·τ² + d_k·τ³
//    where:
//      a_k = y_k
//      c_k = M_k · dt²/2
//      d_k = (M_{k+1} - M_k) · dt²/6
//      b_k = (y_{k+1} - y_k) - (2·M_k + M_{k+1}) · dt²/6
//
//  Maximum found per segment by solving s'_k(τ) = 0 (quadratic in τ).
//  Error bounds (Eq. 17 of paper): O(dt³) for gradient vs O(1) for naive.
// =============================================================================


// Solve tridiagonal system for natural cubic spline second derivatives M[].
// Input:  y[0..n],  n+1 data points, uniform spacing dt
// Output: M[0..n],  second derivatives (M[0]=M[n]=0 for natural spline)
static void spline_fit_natural(const double *y, int n, double dt, double *M)
{
    M[0] = 0.0;
    M[n] = 0.0;

    if (n < 2) return;

    // Fixed-size arrays — sized for PCBF_H which is the max n we ever pass
    double c_arr[PCBF_H];
    double d_arr[PCBF_H];
    double w_arr[PCBF_H];

    int m = n - 1;
    double inv_dt2 = 1.0 / (dt * dt);
    double scale   = 6.0 * inv_dt2;

    // Forward sweep
    c_arr[0] = 1.0 / 4.0;
    d_arr[0] = scale * (y[0] - 2.0*y[1] + y[2]) / 4.0;

    for (int i = 1; i < m; i++) {
        double denom = 4.0 - c_arr[i-1];
        c_arr[i] = 1.0 / denom;
        d_arr[i] = (scale * (y[i] - 2.0*y[i+1] + y[i+2]) - d_arr[i-1]) / denom;
    }

    // Back substitution
    w_arr[m-1] = d_arr[m-1];
    for (int i = m-2; i >= 0; i--)
        w_arr[i] = d_arr[i] - c_arr[i] * w_arr[i+1];

    for (int i = 0; i < m; i++)
        M[i+1] = w_arr[i];
}


// Find the maximum of the cubic spline over the entire domain [0, n*dt].
// Returns the maximum value and sets tau_star_out = continuous fractional
// position within the winning segment (0 <= tau_star_out <= 1),
// and k_star_out = winning segment index.
static double spline_find_max(
    const double *y, const double *M, int n, double dt,
    int *k_star_out, double *tau_star_out)
{
    double val_max = y[0];
    int    k_star  = 0;
    double tau_star = 0.0;

    double dt2 = dt * dt;

    for (int k = 0; k < n; k++)
    {
        // Spline coefficients on segment k
        double a = y[k];
        double c = M[k]   * dt2 / 2.0;
        double d = (M[k+1] - M[k]) * dt2 / 6.0;
        double b = (y[k+1] - y[k]) - c - d;

        // s'(τ) = b + 2c·τ + 3d·τ² = 0
        // Coefficients of quadratic: A·τ² + B·τ + C = 0
        double A = 3.0 * d;
        double B = 2.0 * c;
        double C = b;

        // Check endpoints of segment
        double v0 = a;
        double v1 = a + b + c + d;
        if (v0 > val_max) { val_max = v0; k_star = k; tau_star = 0.0; }
        if (v1 > val_max) { val_max = v1; k_star = k; tau_star = 1.0; }

        // Check interior critical points (roots of quadratic)
        if (fabs(A) > 1e-12)
        {
            double disc = B*B - 4.0*A*C;
            if (disc >= 0.0) {
                double sq = sqrt(disc);
                double t1 = (-B + sq) / (2.0*A);
                double t2 = (-B - sq) / (2.0*A);
                for (int r = 0; r < 2; r++) {
                    double tau = (r == 0) ? t1 : t2;
                    if (tau > 0.0 && tau < 1.0) {
                        double v = a + tau*(b + tau*(c + tau*d));
                        if (v > val_max) {
                            val_max = v; k_star = k; tau_star = tau;
                        }
                    }
                }
            }
        }
        else if (fabs(B) > 1e-12)
        {
            // Linear derivative: one root
            double tau = -C / B;
            if (tau > 0.0 && tau < 1.0) {
                double v = a + tau*(b + tau*(c + tau*d));
                if (v > val_max) {
                    val_max = v; k_star = k; tau_star = tau;
                }
            }
        }
    }

    *k_star_out  = k_star;
    *tau_star_out = tau_star;
    return val_max;
}


// =============================================================================
//  pcbf_rollout_and_value()  —  V2: uses cubic splines (Section III.C)
//
//  Rolls out the design policy for H steps.
//  For each obstacle i:
//    1. Collects discrete h values h_k and sensitivity matrices S_k
//    2. Fits a natural cubic spline to {h_k}
//    3. Finds continuous maximum analytically (no gradient discontinuity)
//    4. Interpolates S at the continuous max time for accurate gradient
//
//  Outputs:
//    V_val[i]   : barrier value (spline maximum of h_i along rollout)
//    dV_dx0[i]  : gradient ∇V w.r.t. [x0, y0, θ0]
// =============================================================================
static void pcbf_rollout_and_value(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS],
    double dV_dx0[N_OBS][3])
{
    // Design policy input
    double ux, uy, uw;
    if (pcbf_policy == POLICY_ZERO_VELOCITY) {
        // π = zero: robot stays put, rollout is stationary
        ux = 0.0;  uy = 0.0;  uw = 0.0;
    } else if (pcbf_policy == POLICY_NOMINAL) {
        // π = nominal: hold current nominal command constant over horizon
        ux = vx_nom;  uy = vy_nom;  uw = wz_nom;
    } else {
        // π = push_vx: open-loop constant forward push in current heading
        // No lateral or angular component — purely probes what is ahead.
        // This is smooth and time-invariant, avoiding the discontinuities
        // that POLICY_NOMINAL suffers when the nominal command switches.
        ux = vx_push;  uy = 0.0;  uw = 0.0;
    }

    // -------------------------------------------------------------------------
    // Storage for full trajectory and sensitivity matrices
    //
    // xs, ys, ths : state at each step 0..H
    // h_traj[i][k]: h value for obstacle i at step k  (for spline fitting)
    // S_all[k]    : sensitivity matrix ∂x_k/∂x0       (for gradient)
    // -------------------------------------------------------------------------
    double xs[PCBF_H+1], ys[PCBF_H+1], ths[PCBF_H+1];
    double h_traj[N_OBS][PCBF_H+1];
    double S_all[PCBF_H+1][3][3];

    xs[0] = x0;  ys[0] = y0;  ths[0] = theta0;

    // S_0 = I  (∂x0/∂x0 = identity)
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            S_all[0][r][c] = (r == c) ? 1.0 : 0.0;

    // Evaluate h at step 0
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { h_traj[i][0] = -1e9; continue; }
        double dx = x0 - obstacles[i].cx;
        double dy = y0 - obstacles[i].cy;
        h_traj[i][0] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
    }

    // -------------------------------------------------------------------------
    // Forward rollout: state + sensitivity + h values
    // -------------------------------------------------------------------------
    for (int k = 0; k < PCBF_H; k++)
    {
        double c = cos(ths[k]);
        double s = sin(ths[k]);

        // Unicycle Euler step
        xs[k+1]  = xs[k]  + PCBF_DT * (c*ux - s*uy);
        ys[k+1]  = ys[k]  + PCBF_DT * (s*ux + c*uy);
        ths[k+1] = ths[k] + PCBF_DT * uw;

        // Jacobian A_k = I + dt * ∂f/∂x
        double df_dth_x = -s*ux - c*uy;
        double df_dth_y =  c*ux - s*uy;

        double A[3][3] = {
            {1, 0, PCBF_DT * df_dth_x},
            {0, 1, PCBF_DT * df_dth_y},
            {0, 0, 1}
        };

        // S_{k+1} = A_k · S_k
        for (int r = 0; r < 3; r++)
            for (int cc = 0; cc < 3; cc++) {
                S_all[k+1][r][cc] = 0.0;
                for (int m = 0; m < 3; m++)
                    S_all[k+1][r][cc] += A[r][m] * S_all[k][m][cc];
            }

        // h_i at this step (paper convention: h > 0 inside obstacle)
        for (int i = 0; i < N_OBS; i++) {
            if (!obstacles[i].active) { h_traj[i][k+1] = -1e9; continue; }
            double dx = xs[k+1] - obstacles[i].cx;
            double dy = ys[k+1] - obstacles[i].cy;
            h_traj[i][k+1] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
        }
    }

    // -------------------------------------------------------------------------
    // For each obstacle: fit cubic spline, find continuous max, compute gradient
    // -------------------------------------------------------------------------
    double M_spline[PCBF_H+1];   // second derivatives of spline

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) {
            V_val[i] = -1e9;
            dV_dx0[i][0] = 0; dV_dx0[i][1] = 0; dV_dx0[i][2] = 0;
            continue;
        }

        // Step 1: Fit natural cubic spline to h_traj[i][0..H]
        spline_fit_natural(h_traj[i], PCBF_H, PCBF_DT, M_spline);

        // Step 2: Find continuous maximum analytically
        //         Returns V = max of spline, k* = segment, τ* ∈ [0,1]
        int    k_star;
        double tau_star;
        V_val[i] = spline_find_max(h_traj[i], M_spline, PCBF_H, PCBF_DT,
                                   &k_star, &tau_star);

        // Step 2b: Enforce both B(x) and h(x) — take worst of current state
        //          and future rollout (paper Theorem 1 requires both).
        //
        //   V_T(x0) >= h(x0) in theory, but numerical errors can break this.
        //   Explicitly taking fmax ensures:
        //     - If robot is currently inside obstacle (h(x0) > 0), QP activates
        //     - If rollout sees future danger (V_T > h(x0)), QP activates
        //   When h(x0) > V_T, k* resets to 0 so gradV = gradh(x0) = standard CBF.
        {
            double h_current = h_traj[i][0];   // h at x0, already computed
            if (h_current > V_val[i]) {
                V_val[i] = h_current;
                k_star   = 0;
                tau_star = 0.0;
            }
        }

        // Step 3: Interpolate sensitivity S at continuous time t* = (k* + τ*)*dt
        //         S(τ*) ≈ (1-τ*)·S_{k*} + τ*·S_{k*+1}
        //         This gives accurate gradient at the true continuous maximum.
        double S_star[3][3];
        for (int r = 0; r < 3; r++)
            for (int cc = 0; cc < 3; cc++)
                S_star[r][cc] = (1.0-tau_star)*S_all[k_star][r][cc]
                                + tau_star*S_all[k_star+1][r][cc];

        // Step 4: Interpolate state at t* for accurate ∇h evaluation
        double x_star = (1.0-tau_star)*xs[k_star] + tau_star*xs[k_star+1];
        double y_star = (1.0-tau_star)*ys[k_star] + tau_star*ys[k_star+1];

        // ∇h at continuous t* (paper convention: h = r² - dx² - dy²)
        //   ∂h/∂x = -2*(x - cx)
        //   ∂h/∂y = -2*(y - cy)
        double dx_star = x_star - obstacles[i].cx;
        double dy_star = y_star - obstacles[i].cy;
        double grad_h[3];
        grad_h[0] = -2.0*dx_star;
        grad_h[1] = -2.0*dy_star;
        grad_h[2] =  0.0;

        // Step 5: dV/dx0 = ∇h(x*) · S(τ*)
        for (int j = 0; j < 3; j++) {
            dV_dx0[i][j] = 0.0;
            for (int m = 0; m < 3; m++)
                dV_dx0[i][j] += grad_h[m] * S_star[m][j];
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
    pcbf_settings.max_iter      = 2000;

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
    
    // ---- LOGGER — place here, all variables in scope ----
    {
        static int pcbf_log_step = 0;
        pcbf_open_log();
        pcbf_log_step++;

        if (pcbf_log_file != NULL) {
            double h_cur[N_OBS];
            for (int i = 0; i < N_OBS; i++) {
                double dx0 = x - obstacles[i].cx;
                double dy0 = y - obstacles[i].cy;
                h_cur[i] = obstacles[i].active
                    ? obstacles[i].r*obstacles[i].r - dx0*dx0 - dy0*dy0
                    : 0.0;
            }
            fprintf(pcbf_log_file,
                "%d,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                pcbf_log_step, x, y, theta,
                h_cur[0], h_cur[1], h_cur[2], h_cur[3],
                V_val[0],  V_val[1],  V_val[2],  V_val[3],
                dV_dx0[0][0], dV_dx0[0][1], dV_dx0[0][2],
                dV_dx0[1][0], dV_dx0[1][1], dV_dx0[1][2],
                dV_dx0[2][0], dV_dx0[2][1], dV_dx0[2][2],
                dV_dx0[3][0], dV_dx0[3][1], dV_dx0[3][2]);
            fflush(pcbf_log_file);
        }
    }
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
