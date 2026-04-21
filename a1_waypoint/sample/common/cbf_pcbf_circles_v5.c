#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  Policy Control Barrier Function (PCBF) — Circle Obstacle Avoidance  (v5)
//
//  Changes from v4:
//    - Clarified that finite-difference gradient IS the correct equivalent of
//      the paper's autodiff: differentiates the full spline-max value functional
//      V w.r.t. x0, including the spline fitting and max-finding steps.
//    - Fixed buffer to static const to prevent multiple-include linker issues.
//    - Removed stale "sensitivity matrix" comment from unicycle section.
//    - Added explicit comment explaining why FD > sensitivity matrix here.
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
//      h_i(x,y) = r_i^2 - (x-cx_i)^2 - (y-cy_i)^2
//      Positive = violation (inside), negative = safe (outside).
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
//    Gradient of V (v4 implementation):
//      ∇V is estimated numerically by central finite differences of the full
//      spline-max value functional with respect to [x0, y0, θ0].
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
static const double buffer = 0.15;

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, 0.8, 1 },
    { -4.0, -1.0, 0.35, 0.8, 1 },
    { -2.0, -1.0, 0.35, 0.8, 1 },
    { -4.0,  1.0, 0.35, 0.8, 1 },
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
// Gradient of V is computed by central finite differences of the full
// spline-max value functional — this differentiates the exact quantity used
// in the QP, equivalent to automatic differentiation as in the paper.
// -----------------------------------------------------------------------------


// =============================================================================
//  CUBIC SPLINE HELPERS  (Section III.C of paper)
//
//  Instead of taking max over discrete h values (naive, causes gradient errors),
//  we fit a natural cubic spline through {h_k} and find the continuous maximum
//  analytically. This resolves the piecewise-constant gradient discontinuity
//  shown in Fig. 2 of the paper for the VALUE computation. In v4, the gradient
//  is then estimated by finite differences of that spline-max value.
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
//  pcbf_rollout_value_only()
//
//  Rolls out the design policy for H steps and computes the spline-smoothed
//  finite-horizon PCBF value V for each obstacle. This function computes only
//  the value, not its gradient.
// =============================================================================
static void pcbf_rollout_value_only(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS])
{
    // Design policy input
    double ux, uy, uw;
    if (pcbf_policy == POLICY_ZERO_VELOCITY) {
        ux = 0.0;  uy = 0.0;  uw = 0.0;
    } else if (pcbf_policy == POLICY_NOMINAL) {
        ux = vx_nom;  uy = vy_nom;  uw = wz_nom;
    } else {
        ux = vx_push;  uy = 0.0;  uw = 0.0;
    }

    double xs[PCBF_H+1], ys[PCBF_H+1], ths[PCBF_H+1];
    double h_traj[N_OBS][PCBF_H+1];
    double M_spline[PCBF_H+1];

    xs[0] = x0;  ys[0] = y0;  ths[0] = theta0;

    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { h_traj[i][0] = -1e9; continue; }
        double dx = x0 - obstacles[i].cx;
        double dy = y0 - obstacles[i].cy;
        h_traj[i][0] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
    }

    for (int k = 0; k < PCBF_H; k++)
    {
        double c = cos(ths[k]);
        double s = sin(ths[k]);

        xs[k+1]  = xs[k]  + PCBF_DT * (c*ux - s*uy);
        ys[k+1]  = ys[k]  + PCBF_DT * (s*ux + c*uy);
        ths[k+1] = ths[k] + PCBF_DT * uw;

        for (int i = 0; i < N_OBS; i++) {
            if (!obstacles[i].active) { h_traj[i][k+1] = -1e9; continue; }
            double dx = xs[k+1] - obstacles[i].cx;
            double dy = ys[k+1] - obstacles[i].cy;
            h_traj[i][k+1] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
        }
    }

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) {
            V_val[i] = -1e9;
            continue;
        }

        spline_fit_natural(h_traj[i], PCBF_H, PCBF_DT, M_spline);

        int k_star;
        double tau_star;
        V_val[i] = spline_find_max(h_traj[i], M_spline, PCBF_H, PCBF_DT,
                                   &k_star, &tau_star);

        // Preserve the current-state barrier value so V >= h(x0) numerically.
        {
            double h_current = h_traj[i][0];
            if (h_current > V_val[i]) {
                V_val[i] = h_current;
            }
        }
    }
}


// =============================================================================
//  pcbf_rollout_and_value()  — numerical autodiff via central finite differences
//
//  The paper (Algorithm 1, line 8) computes ∇V using automatic differentiation
//  through the rollout + spline pipeline. In C without an autodiff library,
//  central finite differences of the full value functional is the exact
//  equivalent:
//
//    ∂V/∂x0 ≈ [V(x0+ε,y0,θ0) - V(x0-ε,y0,θ0)] / 2ε
//    ∂V/∂y0 ≈ [V(x0,y0+ε,θ0) - V(x0,y0-ε,θ0)] / 2ε
//    ∂V/∂θ0 ≈ [V(x0,y0,θ0+ε) - V(x0,y0,θ0-ε)] / 2ε
//
//  This differentiates the SAME quantity used in the QP (spline-smoothed max),
//  not an approximation of it. Cost: 6 extra rollouts per control step.
//
//  Why this is better than the sensitivity matrix (v2/v3):
//    The sensitivity matrix approximates ∇V by ∇h(x*) · S(τ*), where S(τ*)
//    is linearly interpolated. This has O(dt) error from the interpolation.
//    FD differentiates V directly — the spline max itself — so it correctly
//    handles cases where the argmax t* moves as x0 changes.
// =============================================================================
static void pcbf_rollout_and_value(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS],
    double dV_dx0[N_OBS][3])
{
    // Base value
    pcbf_rollout_value_only(x0, y0, theta0,
                            vx_nom, vy_nom, wz_nom,
                            V_val);

    // Finite-difference step sizes.
    // Position eps should be small compared to workspace scale but not so small
    // that floating-point cancellation dominates. Theta eps is in radians.
    const double eps_pos = 1e-4;
    const double eps_th  = 1e-5;

    double V_plus[N_OBS], V_minus[N_OBS];

    // d/dx0
    pcbf_rollout_value_only(x0 + eps_pos, y0, theta0,
                            vx_nom, vy_nom, wz_nom, V_plus);
    pcbf_rollout_value_only(x0 - eps_pos, y0, theta0,
                            vx_nom, vy_nom, wz_nom, V_minus);
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { dV_dx0[i][0] = 0.0; continue; }
        dV_dx0[i][0] = (V_plus[i] - V_minus[i]) / (2.0 * eps_pos);
    }

    // d/dy0
    pcbf_rollout_value_only(x0, y0 + eps_pos, theta0,
                            vx_nom, vy_nom, wz_nom, V_plus);
    pcbf_rollout_value_only(x0, y0 - eps_pos, theta0,
                            vx_nom, vy_nom, wz_nom, V_minus);
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { dV_dx0[i][1] = 0.0; continue; }
        dV_dx0[i][1] = (V_plus[i] - V_minus[i]) / (2.0 * eps_pos);
    }

    // d/dtheta0
    pcbf_rollout_value_only(x0, y0, theta0 + eps_th,
                            vx_nom, vy_nom, wz_nom, V_plus);
    pcbf_rollout_value_only(x0, y0, theta0 - eps_th,
                            vx_nom, vy_nom, wz_nom, V_minus);
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { dV_dx0[i][2] = 0.0; continue; }
        dV_dx0[i][2] = (V_plus[i] - V_minus[i]) / (2.0 * eps_th);
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
