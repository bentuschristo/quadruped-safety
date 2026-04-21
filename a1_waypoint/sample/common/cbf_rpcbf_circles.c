#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  Robust Policy Control Barrier Function (RPCBF) — Circle Obstacle Avoidance
//
//  Based on:
//    Knoedler et al., "Safety on the Fly: Constructing Robust Safety Filters
//    via Policy Control Barrier Functions at Runtime", IEEE RA-L 2025.
//
//  HOW THIS EXTENDS PCBF (cbf_pcbf_circles_v2.c):
//
//  PCBF uses a single rollout under the design policy to compute V.
//  RPCBF runs N rollouts, each with a different disturbance trajectory sampled
//  from the bounded disturbance set D = [d_min, d_max]. The barrier value is
//  the WORST CASE across all N samples (Eq. 20 of paper):
//
//    V_T,N(x0) = max_{i=1..N}  sup_{0≤t<T}  h(x_t^i)
//
//  where x_t^i is the state under disturbance sample i.
//  The gradient ∇V_T,N comes from the worst-case sample (the one achieving max).
//
//  THE DISTURBED UNICYCLE MODEL:
//    ẋ = cos(θ)·vx - sin(θ)·vy + d_x
//    ẏ = sin(θ)·vx + cos(θ)·vy + d_y
//    θ̇ = wz + d_θ
//
//  where d = [d_x, d_y, d_θ] is bounded: d_min ≤ d ≤ d_max.
//  Disturbance models unmodelled dynamics, terrain effects, gait imperfections.
//
//  SAMPLING STRATEGY (paper Section III.D):
//  Mix of:
//    - Vertices of the disturbance box (worst-case corners) — always included
//    - Uniform random samples — fill remaining budget
//  For a 3D disturbance box there are 2³=8 vertices. If N_SAMPLES >= 8,
//  all vertices are included; remainder are random.
//
//  QP (same structure as PCBF, solved by OSQP):
//    min   ½ (u - u₀)ᵀ W (u - u₀)  +  p_slack · δ²
//    s.t.  A · u  ≤  b + δ          (RPCBF safety, soft via slack δ)
//          u_min ≤ u ≤ u_max        (input limits, hard)
//          δ ≥ 0
//
// =============================================================================


// ---- Obstacle definition ----------------------------------------------------
typedef struct {
    double cx;     // center x (world frame)
    double cy;     // center y (world frame)
    double r;      // safety radius
    double alpha;  // CBF gain α
    int    active; // 1 = on, 0 = off
} CircleCBF;

// ---- Edit obstacles here ----------------------------------------------------
#define N_OBS 4

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, 1.0, 1 },
    { -4.0, -1.0, 0.35, 1.0, 1 },
    { -2.0, -1.0, 0.35, 1.0, 1 },
    { -4.0,  1.0, 0.35, 1.0, 1 },
};
// -----------------------------------------------------------------------------


// ---- Rollout parameters -----------------------------------------------------
#define PCBF_H      20      // horizon steps
#define PCBF_DT     0.1     // rollout timestep (s)   →  T = H × dt = 2.0 s

#define POLICY_ZERO_VELOCITY  0
#define POLICY_NOMINAL        1
static int pcbf_policy = POLICY_NOMINAL;  // <-- change here to switch
// -----------------------------------------------------------------------------


// ---- Disturbance parameters (Section III.D of paper) -----------------------
//
//  d = [d_x, d_y, d_θ] enters the unicycle as additive world-frame noise:
//    ẋ += d_x,   ẏ += d_y,   θ̇ += d_θ
//
//  Set to zero to recover the non-robust PCBF behavior.
//  Tune to match your quadruped's observed position/heading drift.
//
//  N_SAMPLES = number of disturbance trajectories per control step.
//  The paper used N=64 for hardware. Start with 8 (the 8 box vertices).
//  More samples = more robust but more compute.
//
#define N_SAMPLES   8       // number of disturbance rollouts

static double d_min[3] = { -0.02, -0.02, -0.01 };  // [x, y, theta] lower bound
static double d_max[3] = {  0.02,  0.02,  0.01 };  // [x, y, theta] upper bound
// -----------------------------------------------------------------------------


// ---- Input limits -----------------------------------------------------------
static double u_max[3] = {  1.0,  0.3,  1.0 };
static double u_min[3] = { -1.0, -0.3, -1.0 };

// ---- Slack penalty ----------------------------------------------------------
static double p_slack = 1e4;
// -----------------------------------------------------------------------------


// ---- Simple LCG random number generator (no stdlib dependency) --------------
// Generates uniform float in [0, 1]. Seeded once at first call.
static unsigned int rpcbf_rng_state = 12345;
static double rpcbf_rand01(void)
{
    rpcbf_rng_state = rpcbf_rng_state * 1664525u + 1013904223u;
    return (double)(rpcbf_rng_state >> 1) / (double)(1u << 31);
}
// -----------------------------------------------------------------------------


// ---- OSQP persistent state --------------------------------------------------
static OSQPSolver    *pcbf_solver   = NULL;
static OSQPSettings   pcbf_settings;

#define N_VARS  4               // [vx, vy, wz, delta]
#define N_ROWS  (N_OBS + 3 + 1) // CBF rows + box rows + slack >= 0

static OSQPFloat pcbf_P_x[N_VARS];
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


// =============================================================================
//  CUBIC SPLINE HELPERS  (Section III.C of paper)
//  — unchanged from cbf_pcbf_circles_v2.c —
// =============================================================================

static void spline_fit_natural(const double *y, int n, double dt, double *M)
{
    M[0] = 0.0;
    M[n] = 0.0;
    if (n < 2) return;

    double c_arr[PCBF_H];
    double d_arr[PCBF_H];
    double w_arr[PCBF_H];

    int m = n - 1;
    double inv_dt2 = 1.0 / (dt * dt);
    double scale   = 6.0 * inv_dt2;

    c_arr[0] = 1.0 / 4.0;
    d_arr[0] = scale * (y[0] - 2.0*y[1] + y[2]) / 4.0;

    for (int i = 1; i < m; i++) {
        double denom = 4.0 - c_arr[i-1];
        c_arr[i] = 1.0 / denom;
        d_arr[i] = (scale * (y[i] - 2.0*y[i+1] + y[i+2]) - d_arr[i-1]) / denom;
    }

    w_arr[m-1] = d_arr[m-1];
    for (int i = m-2; i >= 0; i--)
        w_arr[i] = d_arr[i] - c_arr[i] * w_arr[i+1];

    for (int i = 0; i < m; i++)
        M[i+1] = w_arr[i];
}


static double spline_find_max(
    const double *y, const double *M, int n, double dt,
    int *k_star_out, double *tau_star_out)
{
    double val_max  = y[0];
    int    k_star   = 0;
    double tau_star = 0.0;
    double dt2 = dt * dt;

    for (int k = 0; k < n; k++)
    {
        double a = y[k];
        double c = M[k]   * dt2 / 2.0;
        double d = (M[k+1] - M[k]) * dt2 / 6.0;
        double b = (y[k+1] - y[k]) - c - d;

        double A = 3.0 * d;
        double B = 2.0 * c;
        double C = b;

        double v0 = a;
        double v1 = a + b + c + d;
        if (v0 > val_max) { val_max = v0; k_star = k; tau_star = 0.0; }
        if (v1 > val_max) { val_max = v1; k_star = k; tau_star = 1.0; }

        if (fabs(A) > 1e-12) {
            double disc = B*B - 4.0*A*C;
            if (disc >= 0.0) {
                double sq = sqrt(disc);
                double t1 = (-B + sq) / (2.0*A);
                double t2 = (-B - sq) / (2.0*A);
                for (int r = 0; r < 2; r++) {
                    double tau = (r == 0) ? t1 : t2;
                    if (tau > 0.0 && tau < 1.0) {
                        double v = a + tau*(b + tau*(c + tau*d));
                        if (v > val_max) { val_max = v; k_star = k; tau_star = tau; }
                    }
                }
            }
        } else if (fabs(B) > 1e-12) {
            double tau = -C / B;
            if (tau > 0.0 && tau < 1.0) {
                double v = a + tau*(b + tau*(c + tau*d));
                if (v > val_max) { val_max = v; k_star = k; tau_star = tau; }
            }
        }
    }

    *k_star_out   = k_star;
    *tau_star_out = tau_star;
    return val_max;
}


// =============================================================================
//  single_rollout()
//
//  Performs ONE disturbed rollout for H steps and fills:
//    h_out[i][k]    : h_i value at step k for all obstacles
//    S_out[k][r][c] : sensitivity matrix ∂x_k/∂x0 at step k
//    xs_out, ys_out : world-frame positions at each step
//
//  The disturbance d[k] = [d_x, d_y, d_θ] is applied additively each step:
//    ẋ += d_x,  ẏ += d_y,  θ̇ += d_θ
//
//  Note: sensitivity S does NOT include disturbance terms (disturbances are
//  treated as external — the gradient is w.r.t. initial state x0, not d).
//  This is consistent with the paper's formulation.
// =============================================================================
static void single_rollout(
    double x0, double y0, double theta0,
    double ux, double uy, double uw,
    const double d[PCBF_H][3],
    double h_out[N_OBS][PCBF_H+1],
    double S_out[PCBF_H+1][3][3],
    double xs_out[PCBF_H+1],
    double ys_out[PCBF_H+1])
{
    double ths[PCBF_H+1];
    xs_out[0] = x0;  ys_out[0] = y0;  ths[0] = theta0;

    // S_0 = I
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            S_out[0][r][c] = (r == c) ? 1.0 : 0.0;

    // Evaluate h at initial state
    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { h_out[i][0] = -1e9; continue; }
        double dx = x0 - obstacles[i].cx;
        double dy = y0 - obstacles[i].cy;
        h_out[i][0] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
    }

    for (int k = 0; k < PCBF_H; k++)
    {
        double c = cos(ths[k]);
        double s = sin(ths[k]);

        // Disturbed Euler step
        // ẋ = c·vx - s·vy + d_x
        // ẏ = s·vx + c·vy + d_y
        // θ̇ = wz + d_θ
        xs_out[k+1] = xs_out[k] + PCBF_DT * (c*ux - s*uy + d[k][0]);
        ys_out[k+1] = ys_out[k] + PCBF_DT * (s*ux + c*uy + d[k][1]);
        ths[k+1]    = ths[k]    + PCBF_DT * (uw + d[k][2]);

        // Sensitivity S_{k+1} = A_k · S_k
        // A_k = I + dt * ∂f/∂x  (disturbance does not affect Jacobian
        //                         since d enters as additive constant)
        double df_dth_x = -s*ux - c*uy;
        double df_dth_y =  c*ux - s*uy;

        double A[3][3] = {
            {1, 0, PCBF_DT * df_dth_x},
            {0, 1, PCBF_DT * df_dth_y},
            {0, 0, 1}
        };

        for (int r = 0; r < 3; r++)
            for (int cc = 0; cc < 3; cc++) {
                S_out[k+1][r][cc] = 0.0;
                for (int m = 0; m < 3; m++)
                    S_out[k+1][r][cc] += A[r][m] * S_out[k][m][cc];
            }

        // h at this step
        for (int i = 0; i < N_OBS; i++) {
            if (!obstacles[i].active) { h_out[i][k+1] = -1e9; continue; }
            double dx = xs_out[k+1] - obstacles[i].cx;
            double dy = ys_out[k+1] - obstacles[i].cy;
            h_out[i][k+1] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
        }
    }
}


// =============================================================================
//  rpcbf_rollout_and_value()
//
//  Runs N_SAMPLES disturbed rollouts and computes the RPCBF value function:
//
//    V_T,N(x0) = max_{i=1..N}  spline_max(h along rollout i)   (Eq. 20)
//
//  The gradient ∇V_T,N comes from the worst-case sample (argmax over i).
//
//  Disturbance sampling (paper Section III.D):
//    - First 8 samples: all vertices of the 3D disturbance box {d_min,d_max}^3
//      (constant disturbance = worst-case corner applied every step)
//    - Remaining samples (if N_SAMPLES > 8): uniform random from [d_min, d_max]
//      (constant random disturbance per sample, different each control step)
//
//  Outputs:
//    V_val[i]    : worst-case barrier value for obstacle i across all samples
//    dV_dx0[i]   : gradient of V_val[i] w.r.t. [x0, y0, θ0]
// =============================================================================
static void rpcbf_rollout_and_value(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS],
    double dV_dx0[N_OBS][3])
{
    // Design policy input
    double ux, uy, uw;
    if (pcbf_policy == POLICY_ZERO_VELOCITY) {
        ux = 0.0;  uy = 0.0;  uw = 0.0;
    } else {
        ux = vx_nom;  uy = vy_nom;  uw = wz_nom;
    }

    // Initialize worst-case tracking across all samples
    double V_best[N_OBS];      // worst-case V per obstacle (max over samples)
    int    best_sample[N_OBS]; // which sample achieved the worst case
    for (int i = 0; i < N_OBS; i++) {
        V_best[i]      = -1e9;
        best_sample[i] = 0;
    }

    // Storage for all samples' trajectories — we keep the worst-case one
    // for gradient computation at the end.
    // To save memory we only store the worst-case sample's data per obstacle.
    double h_best[N_OBS][PCBF_H+1];
    double S_best[N_OBS][PCBF_H+1][3][3];
    double xs_best[N_OBS][PCBF_H+1];
    double ys_best[N_OBS][PCBF_H+1];

    // Temporary storage for current sample
    double h_tmp[N_OBS][PCBF_H+1];
    double S_tmp[PCBF_H+1][3][3];
    double xs_tmp[PCBF_H+1];
    double ys_tmp[PCBF_H+1];
    double M_spline[PCBF_H+1];

    // -------------------------------------------------------------------------
    // Loop over N_SAMPLES disturbance trajectories
    // -------------------------------------------------------------------------
    for (int samp = 0; samp < N_SAMPLES; samp++)
    {
        // ---------------------------------------------------------------------
        // Build disturbance trajectory for this sample
        //
        // Strategy: constant disturbance applied at every step.
        // First 8 samples = all 8 vertices of [d_min, d_max]^3 box.
        // Remaining samples = uniform random from [d_min, d_max].
        //
        // Using constant disturbance per sample (not time-varying) is simpler
        // and still captures the worst-case corners as the paper suggests.
        // ---------------------------------------------------------------------
        double d_k[3];  // disturbance for this sample (constant over horizon)

        if (samp < 8) {
            // Vertex sampling: bit j of samp selects d_min or d_max for dim j
            for (int j = 0; j < 3; j++)
                d_k[j] = (samp & (1 << j)) ? d_max[j] : d_min[j];
        } else {
            // Uniform random from [d_min, d_max]
            for (int j = 0; j < 3; j++)
                d_k[j] = d_min[j] + rpcbf_rand01() * (d_max[j] - d_min[j]);
        }

        // Build constant disturbance array for single_rollout
        double d_traj[PCBF_H][3];
        for (int k = 0; k < PCBF_H; k++) {
            d_traj[k][0] = d_k[0];
            d_traj[k][1] = d_k[1];
            d_traj[k][2] = d_k[2];
        }

        // ---------------------------------------------------------------------
        // Run the rollout for this disturbance sample
        // ---------------------------------------------------------------------
        single_rollout(x0, y0, theta0, ux, uy, uw,
                       d_traj, h_tmp, S_tmp, xs_tmp, ys_tmp);

        // ---------------------------------------------------------------------
        // For each obstacle: fit spline, find max V, update worst case
        // ---------------------------------------------------------------------
        for (int i = 0; i < N_OBS; i++)
        {
            if (!obstacles[i].active) continue;

            spline_fit_natural(h_tmp[i], PCBF_H, PCBF_DT, M_spline);

            int    k_star_tmp;
            double tau_star_tmp;
            double V_tmp = spline_find_max(h_tmp[i], M_spline, PCBF_H, PCBF_DT,
                                           &k_star_tmp, &tau_star_tmp);

            // Keep this sample if it's the worst case for obstacle i
            if (V_tmp > V_best[i]) {
                V_best[i]      = V_tmp;
                best_sample[i] = samp;

                // Save this sample's trajectory for gradient computation
                memcpy(h_best[i],  h_tmp[i],  sizeof(h_tmp[i]));
                memcpy(xs_best[i], xs_tmp,     sizeof(xs_tmp));
                memcpy(ys_best[i], ys_tmp,     sizeof(ys_tmp));
                for (int k = 0; k <= PCBF_H; k++)
                    memcpy(S_best[i][k], S_tmp[k], sizeof(S_tmp[k]));
            }
        }
    }

    // =========================================================================
    // Gradient computation from worst-case sample
    //
    // For each obstacle i, re-run spline on the worst-case trajectory to get
    // k* and τ*, then compute ∇V = ∇h(x*) · S(τ*) as in the PCBF.
    // =========================================================================
    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) {
            V_val[i] = -1e9;
            dV_dx0[i][0] = 0; dV_dx0[i][1] = 0; dV_dx0[i][2] = 0;
            continue;
        }

        V_val[i] = V_best[i];

        // Re-fit spline on worst-case trajectory to recover k*, τ*
        spline_fit_natural(h_best[i], PCBF_H, PCBF_DT, M_spline);

        int    k_star;
        double tau_star;
        spline_find_max(h_best[i], M_spline, PCBF_H, PCBF_DT, &k_star, &tau_star);

        // Interpolate sensitivity and state at continuous time t*
        double S_star[3][3];
        for (int r = 0; r < 3; r++)
            for (int cc = 0; cc < 3; cc++)
                S_star[r][cc] = (1.0-tau_star)*S_best[i][k_star][r][cc]
                                + tau_star*S_best[i][k_star+1][r][cc];

        double x_star = (1.0-tau_star)*xs_best[i][k_star]
                        + tau_star*xs_best[i][k_star+1];
        double y_star = (1.0-tau_star)*ys_best[i][k_star]
                        + tau_star*ys_best[i][k_star+1];

        // ∇h at x* (paper convention: h = r² - dx² - dy²)
        double dx_star = x_star - obstacles[i].cx;
        double dy_star = y_star - obstacles[i].cy;
        double grad_h[3];
        grad_h[0] = -2.0*dx_star;
        grad_h[1] = -2.0*dy_star;
        grad_h[2] =  0.0;

        // ∇V = ∇h · S*
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
    int nr = nc + 3 + 1;

    pcbf_P_x[0] = (OSQPFloat)wx;
    pcbf_P_x[1] = (OSQPFloat)wy;
    pcbf_P_x[2] = (OSQPFloat)ww;
    pcbf_P_x[3] = (OSQPFloat)p_slack;
    for (int j = 0; j < N_VARS; j++) { pcbf_P_i[j] = j; pcbf_P_p[j] = j; }
    pcbf_P_p[N_VARS] = N_VARS;

    pcbf_P_mat.m = N_VARS; pcbf_P_mat.n = N_VARS;
    pcbf_P_mat.nz = -1;    pcbf_P_mat.nzmax = N_VARS;
    pcbf_P_mat.x = pcbf_P_x; pcbf_P_mat.i = pcbf_P_i; pcbf_P_mat.p = pcbf_P_p;

    int nnz = 0;
    // col 0 (vx): nc CBF rows + 1 box row
    pcbf_A_p[0] = nnz;
    for (int row = 0; row < nc; row++) { pcbf_A_x[nnz]=0.0f; pcbf_A_i[nnz]=row; nnz++; }
    pcbf_A_x[nnz]=1.0f; pcbf_A_i[nnz]=nc+0; nnz++;

    // col 1 (vy)
    pcbf_A_p[1] = nnz;
    for (int row = 0; row < nc; row++) { pcbf_A_x[nnz]=0.0f; pcbf_A_i[nnz]=row; nnz++; }
    pcbf_A_x[nnz]=1.0f; pcbf_A_i[nnz]=nc+1; nnz++;

    // col 2 (wz)
    pcbf_A_p[2] = nnz;
    for (int row = 0; row < nc; row++) { pcbf_A_x[nnz]=0.0f; pcbf_A_i[nnz]=row; nnz++; }
    pcbf_A_x[nnz]=1.0f; pcbf_A_i[nnz]=nc+2; nnz++;

    // col 3 (delta): -1 for each CBF row, -1 for slack >= 0 row
    pcbf_A_p[3] = nnz;
    for (int row = 0; row < nc; row++) { pcbf_A_x[nnz]=-1.0f; pcbf_A_i[nnz]=row; nnz++; }
    pcbf_A_x[nnz]=-1.0f; pcbf_A_i[nnz]=nc+3; nnz++;

    pcbf_A_p[4] = nnz;

    pcbf_A_mat.m = nr; pcbf_A_mat.n = N_VARS;
    pcbf_A_mat.nz = -1; pcbf_A_mat.nzmax = nnz;
    pcbf_A_mat.x = pcbf_A_x; pcbf_A_mat.i = pcbf_A_i; pcbf_A_mat.p = pcbf_A_p;

    for (int j = 0; j < N_VARS; j++) pcbf_q[j] = 0.0f;
    for (int i = 0; i < nc; i++)  { pcbf_l[i]    = (OSQPFloat)(-OSQP_INFTY); pcbf_u[i]    = 0.0f; }
    for (int j = 0; j < 3;  j++)  { pcbf_l[nc+j] = (OSQPFloat)u_min[j];      pcbf_u[nc+j] = (OSQPFloat)u_max[j]; }
    pcbf_l[nc+3] = (OSQPFloat)(-OSQP_INFTY);
    pcbf_u[nc+3] = 0.0f;

    osqp_set_default_settings(&pcbf_settings);
    pcbf_settings.verbose       = 0;
    pcbf_settings.warm_starting = 1;
    pcbf_settings.eps_abs       = 1e-3;
    pcbf_settings.eps_rel       = 1e-3;
    pcbf_settings.max_iter      = 1000;

    OSQPInt flag = osqp_setup(&pcbf_solver, &pcbf_P_mat, pcbf_q,
                               &pcbf_A_mat, pcbf_l, pcbf_u, nr, N_VARS, &pcbf_settings);
    if (flag != 0) {
        fprintf(stderr, "[RPCBF] OSQP setup failed (code %d)\n", (int)flag);
        return 0;
    }
    return 1;
}


// =============================================================================
//  cbf_circle_obstacles_filter()  — drop-in replacement
//
//  Same function signature as all previous CBF files.
//  Change #include to "cbf_rpcbf_circles.c" in my_controller.c.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    const double c = cos(theta);
    const double s = sin(theta);

    // =========================================================================
    // STEP 1: Run N disturbed rollouts, compute worst-case V and ∇V
    //         (Eq. 20 of paper — RPCBF approximation)
    // =========================================================================
    double V_val[N_OBS];
    double dV_dx0[N_OBS][3];

    rpcbf_rollout_and_value(x, y, theta,
                            *vx_cmd, *vy_cmd, *wz_cmd,
                            V_val, dV_dx0);

    // =========================================================================
    // STEP 2: Build QP constraint rows from worst-case V and ∇V
    //
    //  Same as PCBF — the robust extension only changes how V and ∇V are
    //  computed (worst case over N samples), not the QP structure.
    //
    //  CBF-QP row for obstacle i:
    //    A_i · u ≤ b_i
    //  where:
    //    A_i[0] =  dVdx·c + dVdy·s    (coefficient of vx)
    //    A_i[1] = -dVdx·s + dVdy·c    (coefficient of vy)
    //    A_i[2] =  dVdt               (coefficient of wz)
    //    b_i    = -alpha · V           (RHS — tighter = more conservative)
    // =========================================================================
    double A_cbf[N_OBS][3];
    double b_cbf[N_OBS];
    int nc = 0;

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) continue;

        double dVdx = dV_dx0[i][0];
        double dVdy = dV_dx0[i][1];
        double dVdt = dV_dx0[i][2];

        A_cbf[nc][0] =  dVdx*c + dVdy*s;
        A_cbf[nc][1] = -dVdx*s + dVdy*c;
        A_cbf[nc][2] =  dVdt;
        b_cbf[nc]    = -obstacles[i].alpha * V_val[i];

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
    // STEP 4: Update OSQP data
    // =========================================================================
    pcbf_q[0] = (OSQPFloat)(-wx * (*vx_cmd));
    pcbf_q[1] = (OSQPFloat)(-wy * (*vy_cmd));
    pcbf_q[2] = (OSQPFloat)(-ww * (*wz_cmd));
    pcbf_q[3] = 0.0f;

    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < nc; row++)
            pcbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        nnz++;  // skip identity entry
    }

    for (int i = 0; i < nc; i++)
        pcbf_u[i] = (OSQPFloat)b_cbf[i];

    osqp_update_data_vec(pcbf_solver, pcbf_q, pcbf_l, pcbf_u);
    osqp_update_data_mat(pcbf_solver, NULL, NULL, 0, pcbf_A_x, NULL, pcbf_A_mat.nzmax);

    // =========================================================================
    // STEP 5: Solve and extract result
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
            fprintf(stderr, "[RPCBF] slack active: delta=%.4f\n", delta);
        modified = 1;
    }
    else
    {
        fprintf(stderr, "[RPCBF] OSQP did not converge (%s), clamping\n",
                pcbf_solver->info->status);
        *vx_cmd = fmax(u_min[0], fmin(u_max[0], *vx_cmd));
        *vy_cmd = fmax(u_min[1], fmin(u_max[1], *vy_cmd));
        *wz_cmd = fmax(u_min[2], fmin(u_max[2], *wz_cmd));
    }

    return modified;
}
