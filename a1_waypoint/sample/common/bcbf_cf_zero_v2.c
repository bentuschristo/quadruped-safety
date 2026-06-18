#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>
#include <time.h>

// =============================================================================
// Closed-form LSE Backup-CBF filter with CONTROL DYNAMICS and AUGMENTED BACKUP
// ROLLOUT for the high-level quadruped model.
//
// Original plant state:        x = [px, py, theta]
// Actual plant input/state:    u = [vx, vy, omega]
// Surrogate decision input:    uhat = [hat_vx, hat_vy, hat_omega]
// Augmented state:             xhat = [px, py, theta, vx, vy, omega]
//
// Plant:
//      px_dot    = cos(theta) vx - sin(theta) vy
//      py_dot    = sin(theta) vx + cos(theta) vy
//      theta_dot = omega
//
// Control dynamics:
//      u_dot = -A u + A uhat
//
// Zero backup under control dynamics:
//      uhat_B = 0
//      u_dot_B = -A u_B
//
// IMPORTANT:
//      This file fixes the earlier inconsistency where zero backup was treated
//      as instantaneous u=0. Here the backup rollout is performed in the
//      augmented state [x;y;theta;vx;vy;omega], so the robot still moves while
//      the command u decays toward zero.
//
// Variant controlled by BCBF_USE_AUGMENTED_TERMINAL_SET:
//      0: terminal backup set only checks original state clearance.
//      1: terminal backup set checks original clearance AND near-stopped input.
//
// Use only one of the two generated files in common/my_controller.c.
// =============================================================================
// Timing note:
//      cbf_log_qp_time_sec intentionally measures only the closed-form
//      KKT update time (lambda, uhat, mu), not the full filter time.
//      This matches solver-only timing used by QP/OI baselines.
//

#ifndef BCBF_USE_AUGMENTED_TERMINAL_SET
#define BCBF_USE_AUGMENTED_TERMINAL_SET 0
#endif

static int cf_success_count = 0;
static int cf_fail_count = 0;

// Logging compatible with previous files.
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};
static double cbf_log_uhat_d[3]   = {0.0, 0.0, 0.0};
static double cbf_log_uhat[3]     = {0.0, 0.0, 0.0};
static double cbf_log_mu          = 0.0;
static double cbf_log_lambda      = 0.0;
static int    cbf_log_status      = 0;
static double cbf_log_qp_time_sec = 0.0;      // closed-form KKT update only [s]
static double cbf_log_filter_time_sec = 0.0;  // whole safety-filter call [s]
static double cbf_log_total_time_sec  = 0.0;  // accumulated whole-filter time [s]
static double cbf_log_avg_time_sec    = 0.0;  // running average whole-filter time [s]
static double cbf_log_max_time_sec    = 0.0;  // maximum observed whole-filter time [s]
static int    cbf_log_call_count      = 0;    // number of calls to cbf_circle_obstacles_filter
static double cbf_log_h_aug       = 0.0;
static double cbf_log_omega       = 0.0;
static double cbf_log_h_backup_lse = 0.0;
static double cbf_log_h_input_lse  = 0.0;

// --------------------------- Scenario parameters -----------------------------
typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} BCBFObstacle;

// #define BCBF_N_OBS 4
// static BCBFObstacle bcbf_obs[BCBF_N_OBS] = {
//     { -2.0,  1.0, 0.45, 0.5, 1 },
//     { -4.0, -1.0, 0.45, 0.5, 1 },
//     { -2.0, -1.0, 0.45, 0.5, 1 },
//     { -4.0,  1.0, 0.45, 0.5, 1 }
// };

#define BCBF_N_OBS 5

static BCBFObstacle bcbf_obs[BCBF_N_OBS] = {
    {  0.5, 1.0, 0.40, 0.5, 1 },
    { -0.5, 2.0, 0.40, 0.5, 1 },
    {  0.0, 4.0, 0.40, 0.5, 1 },
    {  0.4, 6.0, 0.40, 0.5, 1 },
    { -1.0, 5.5, 0.40, 0.5, 1 }
};

// Actual input limits for u = [vx, vy, omega].
static double bcbf_u_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_min[3] = { -1.0, -0.3, -1.0 };

// Backup set: enlarged-obstacle clearance set.
static const double bcbf_backup_margin = 0.1;

// Terminal near-stop set used only when BCBF_USE_AUGMENTED_TERMINAL_SET = 1.
static const double bcbf_terminal_u_eps[3] = {0.08, 0.04, 0.08};

// LSE gains. Use hierarchical LSE to avoid one giant 86-term soft-min dominating.
static const double bcbf_kappa_backup_lse = 70.0; // Original = 70.0
static const double bcbf_kappa_input_lse  = 70.0; // Original = 70.0
static const double bcbf_kappa_outer_lse  = 50.0; // Original = 50.0

// Outer CBF gain.
// The backup barrier is formed by LSE over the raw backup rollout values H_i,
// then one CBF condition is applied to the final composite barrier.
// This is closer to the original bCBF-QP, where each raw rollout/terminal
// term H_i has its own derivative constraint.
static const double bcbf_alpha_outer  = 0.5; //original = 0.5

// Closed-form slack penalty. Larger gamma discourages mu.
static const double bcbf_gamma_mu = 1.0e4; //original = 1.0e4

// Control dynamics: u_dot = -A u + A uhat.
static double bcbf_cd_a[3] = { 8.0, 8.0, 8.0 }; //original = 0.8

// Input tracking decay for u -> u_d through uhat_d.
static double bcbf_sigma0[3] = { 8.0, 8.0, 8.0 };  //original = 0.4

// Control-loop time step used to integrate u_dot and finite-difference u_d_dot.
static double bcbf_cd_dt = 0.001;

// Directional derivative step used to compute L_fhat H_i.
static const double bcbf_dir_eps = 1.0e-5;

// Final numerical guard. Theory should enforce bounds through barriers.
#define BCBF_CD_HARD_CLIP_OUTPUT 1

#define BCBF_N_SAMPLES 20
static const double bcbf_T_backup = 4.0;

#define BCBF_N_HORIZON_OBS_TERMS (BCBF_N_OBS * BCBF_N_SAMPLES)
#define BCBF_N_TERMINAL_CLEARANCE_TERMS BCBF_N_OBS
#define BCBF_N_TERMINAL_STOP_TERMS 1

#if BCBF_USE_AUGMENTED_TERMINAL_SET
#define BCBF_N_BACKUP_TERMS (BCBF_N_HORIZON_OBS_TERMS + BCBF_N_TERMINAL_CLEARANCE_TERMS + BCBF_N_TERMINAL_STOP_TERMS)
#else
#define BCBF_N_BACKUP_TERMS (BCBF_N_HORIZON_OBS_TERMS + BCBF_N_TERMINAL_CLEARANCE_TERMS)
#endif

#define BCBF_N_INPUT_TERMS  6

// ---------------------------- Utility functions ------------------------------
static inline double cbf_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static inline double bcbf_clip(double v, double lo, double hi)
{
    return fmax(lo, fmin(hi, v));
}

static inline double bcbf_wrap_pi(double a)
{
    while (a >  M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static inline void bcbf_dynamics(double x, double y, double th,
                                 const double u[3], double f[3])
{
    double c = cos(th);
    double s = sin(th);
    f[0] = c * u[0] - s * u[1];
    f[1] = s * u[0] + c * u[1];
    f[2] = u[2];
}

static inline void bcbf_fhat_drift(double x, double y, double theta,
                                   const double u[3], double fhat[6])
{
    double f_old[3];
    bcbf_dynamics(x, y, theta, u, f_old);
    fhat[0] = f_old[0];
    fhat[1] = f_old[1];
    fhat[2] = f_old[2];
    fhat[3] = -bcbf_cd_a[0] * u[0];
    fhat[4] = -bcbf_cd_a[1] * u[1];
    fhat[5] = -bcbf_cd_a[2] * u[2];
}

static inline void bcbf_backup_controller(double x, double y, double th, double ub[3])
{
    // Zero-action backup in surrogate-input space: uhat_B = 0.
    // The actual command u decays by u_dot = -A u during the backup rollout.
    (void)x; (void)y; (void)th;
    ub[0] = 0.0;
    ub[1] = 0.0;
    ub[2] = 0.0;
}

static inline double bcbf_h_obs(int i, double x, double y)
{
    double dx = x - bcbf_obs[i].cx;
    double dy = y - bcbf_obs[i].cy;
    return dx * dx + dy * dy - bcbf_obs[i].r * bcbf_obs[i].r;
}

static inline double bcbf_h_backup_obs(int i, double x, double y)
{
    double dx = x - bcbf_obs[i].cx;
    double dy = y - bcbf_obs[i].cy;
    double rb = bcbf_obs[i].r + bcbf_backup_margin;
    return dx * dx + dy * dy - rb * rb;
}

static inline double bcbf_h_terminal_stop(const double uT[3])
{
    double sx = uT[0] / fmax(bcbf_terminal_u_eps[0], 1e-9);
    double sy = uT[1] / fmax(bcbf_terminal_u_eps[1], 1e-9);
    double sw = uT[2] / fmax(bcbf_terminal_u_eps[2], 1e-9);
    // safe iff normalized squared command <= 1
    return 1.0 - (sx*sx + sy*sy + sw*sw);
}

static inline double bcbf_softmin_from_terms(const double *terms, int n, double kappa)
{
    double zmax = -1.0e300;
    for (int i = 0; i < n; i++) {
        double zi = -kappa * terms[i];
        if (zi > zmax) zmax = zi;
    }

    double sumexp = 0.0;
    for (int i = 0; i < n; i++) {
        sumexp += exp((-kappa * terms[i]) - zmax);
    }

    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int i = 1; i < n; i++) if (terms[i] < terms[imin]) imin = i;
        return terms[imin];
    }

    return -(log(sumexp) + zmax) / kappa;
}

// Jacobian of the augmented zero-backup dynamics f_B(z) = fhat(z).
// State ordering: z = [x, y, theta, vx, vy, omega].
static inline void bcbf_backup_jacobian_aug(const double z[6], double F[6][6])
{
    memset(F, 0, 36 * sizeof(double));

    const double th = z[2];
    const double vx = z[3];
    const double vy = z[4];
    const double c = cos(th);
    const double s = sin(th);

    F[0][2] = -s * vx - c * vy;
    F[0][3] =  c;
    F[0][4] = -s;

    F[1][2] =  c * vx - s * vy;
    F[1][3] =  s;
    F[1][4] =  c;

    F[2][5] = 1.0;

    F[3][3] = -bcbf_cd_a[0];
    F[4][4] = -bcbf_cd_a[1];
    F[5][5] = -bcbf_cd_a[2];
}

static inline void bcbf_matmul6(const double A[6][6],
                                const double B[6][6],
                                double C[6][6])
{
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double s = 0.0;
            for (int k = 0; k < 6; k++) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
    }
}

// Evaluate all augmented backup-flow and terminal terms together with their
// gradients with respect to the CURRENT augmented state z0.
//
// For each scalar term s_l(z0) = h(phi_B(tau_l,z0)),
//     grad s_l(z0) = Phi_B(tau_l,z0)^T grad h(z_B(tau_l)).
//
// This is the same sensitivity-chain-rule construction used by the
// bCBF-QP-CD and bCBF-LSE-QP-CD baselines.
static int bcbf_aug_backup_terms_and_grads(
    double x, double y, double theta,
    const double u0[3],
    double H[BCBF_N_BACKUP_TERMS],
    double gradH[BCBF_N_BACKUP_TERMS][6])
{
    double z[6] = {x, y, theta, u0[0], u0[1], u0[2]};
    double Phi[6][6] = {{0}};
    for (int i = 0; i < 6; i++) Phi[i][i] = 1.0;

    const double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);
    int term = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        for (int i = 0; i < BCBF_N_OBS; i++) {
            const double dx = z[0] - bcbf_obs[i].cx;
            const double dy = z[1] - bcbf_obs[i].cy;
            H[term] = dx * dx + dy * dy - bcbf_obs[i].r * bcbf_obs[i].r;

            const double grad_future[6] = {2.0*dx, 2.0*dy, 0.0, 0.0, 0.0, 0.0};
            for (int j = 0; j < 6; j++) {
                double s = 0.0;
                for (int m = 0; m < 6; m++) s += grad_future[m] * Phi[m][j];
                gradH[term][j] = s;
            }
            term++;
        }

        if (k < BCBF_N_SAMPLES - 1) {
            double fz[6], F[6][6], Phidot[6][6];
            double uu[3] = {z[3], z[4], z[5]};
            bcbf_fhat_drift(z[0], z[1], z[2], uu, fz);
            bcbf_backup_jacobian_aug(z, F);
            bcbf_matmul6(F, Phi, Phidot);

            // Euler propagation, intentionally matching the current CD baselines.
            for (int i = 0; i < 6; i++) z[i] += dt * fz[i];
            z[2] = bcbf_wrap_pi(z[2]);
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    Phi[i][j] += dt * Phidot[i][j];
        }
    }

    // Terminal enlarged-obstacle backup-set terms.
    for (int i = 0; i < BCBF_N_OBS; i++) {
        const double dx = z[0] - bcbf_obs[i].cx;
        const double dy = z[1] - bcbf_obs[i].cy;
        const double rb = bcbf_obs[i].r + bcbf_backup_margin;
        H[term] = dx * dx + dy * dy - rb * rb;

        const double grad_future[6] = {2.0*dx, 2.0*dy, 0.0, 0.0, 0.0, 0.0};
        for (int j = 0; j < 6; j++) {
            double s = 0.0;
            for (int m = 0; m < 6; m++) s += grad_future[m] * Phi[m][j];
            gradH[term][j] = s;
        }
        term++;
    }

#if BCBF_USE_AUGMENTED_TERMINAL_SET
    // Terminal near-stop term and its future-state gradient.
    double uT[3] = {z[3], z[4], z[5]};
    H[term] = bcbf_h_terminal_stop(uT);
    double grad_future[6] = {
        0.0, 0.0, 0.0,
        -2.0*z[3]/(bcbf_terminal_u_eps[0]*bcbf_terminal_u_eps[0]),
        -2.0*z[4]/(bcbf_terminal_u_eps[1]*bcbf_terminal_u_eps[1]),
        -2.0*z[5]/(bcbf_terminal_u_eps[2]*bcbf_terminal_u_eps[2])
    };
    for (int j = 0; j < 6; j++) {
        double s = 0.0;
        for (int m = 0; m < 6; m++) s += grad_future[m] * Phi[m][j];
        gradH[term][j] = s;
    }
    term++;
#endif

    return term;
}

// Stable soft-min value and analytical gradient.
static void bcbf_softmin_value_and_grad(
    const double *terms,
    const double grads[][6],
    int n,
    double kappa,
    double *value,
    double grad[6])
{
    double zmax = -1.0e300;
    for (int i = 0; i < n; i++) {
        const double zi = -kappa * terms[i];
        if (zi > zmax) zmax = zi;
    }

    double sumexp = 0.0;
    for (int j = 0; j < 6; j++) grad[j] = 0.0;

    for (int i = 0; i < n; i++) {
        const double e = exp((-kappa * terms[i]) - zmax);
        sumexp += e;
        for (int j = 0; j < 6; j++) grad[j] += e * grads[i][j];
    }

    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int i = 1; i < n; i++) if (terms[i] < terms[imin]) imin = i;
        *value = terms[imin];
        for (int j = 0; j < 6; j++) grad[j] = grads[imin][j];
        return;
    }

    *value = -(log(sumexp) + zmax) / kappa;
    for (int j = 0; j < 6; j++) grad[j] /= sumexp;
}

static void bcbf_backup_lse_value_and_grad(
    double x, double y, double theta,
    const double u[3],
    double *h_backup,
    double grad_backup[6])
{
    double H[BCBF_N_BACKUP_TERMS];
    double gradH[BCBF_N_BACKUP_TERMS][6];
    const int n = bcbf_aug_backup_terms_and_grads(x, y, theta, u, H, gradH);
    bcbf_softmin_value_and_grad(H, gradH, n, bcbf_kappa_backup_lse,
                                h_backup, grad_backup);
}

static void bcbf_input_lse_value_and_grad(
    const double u[3],
    double *h_input,
    double grad_input[6])
{
    double terms[BCBF_N_INPUT_TERMS];
    double grads[BCBF_N_INPUT_TERMS][6] = {{0}};
    int n = 0;

    terms[n] = bcbf_u_max[0] - u[0]; grads[n][3] = -1.0; n++;
    terms[n] = u[0] - bcbf_u_min[0]; grads[n][3] =  1.0; n++;
    terms[n] = bcbf_u_max[1] - u[1]; grads[n][4] = -1.0; n++;
    terms[n] = u[1] - bcbf_u_min[1]; grads[n][4] =  1.0; n++;
    terms[n] = bcbf_u_max[2] - u[2]; grads[n][5] = -1.0; n++;
    terms[n] = u[2] - bcbf_u_min[2]; grads[n][5] =  1.0; n++;

    bcbf_softmin_value_and_grad(terms, grads, n, bcbf_kappa_input_lse,
                                h_input, grad_input);
}

// Final hierarchical composite barrier and analytical gradient:
//   h_aug = softmin(h_backup, h_input).
static void bcbf_augmented_barrier_value_and_grad(
    double x, double y, double theta,
    const double u[3],
    double *h_aug,
    double grad_aug[6])
{
    double h_backup, h_input;
    double grad_backup[6], grad_input[6];

    bcbf_backup_lse_value_and_grad(x, y, theta, u,
                                   &h_backup, grad_backup);
    bcbf_input_lse_value_and_grad(u, &h_input, grad_input);

    cbf_log_h_backup_lse = h_backup;
    cbf_log_h_input_lse  = h_input;

    double outer_terms[2] = {h_backup, h_input};
    double outer_grads[2][6];
    for (int j = 0; j < 6; j++) {
        outer_grads[0][j] = grad_backup[j];
        outer_grads[1][j] = grad_input[j];
    }

    bcbf_softmin_value_and_grad(outer_terms, outer_grads, 2,
                                bcbf_kappa_outer_lse,
                                h_aug, grad_aug);
}

// Optional reset if you restart an experiment or teleport the robot.
int bcbf_cd_initialized = 0;
static double bcbf_cd_u_state[3]  = {0.0, 0.0, 0.0};
static double bcbf_cd_ud_prev[3]  = {0.0, 0.0, 0.0};

static inline void cbf_circle_obstacles_filter_reset(void)
{
    cf_success_count = 0;
    cf_fail_count = 0;
    cbf_log_status = 0;
    cbf_log_mu = 0.0;
    cbf_log_lambda = 0.0;
    cbf_log_qp_time_sec = 0.0;
    cbf_log_filter_time_sec = 0.0;
    cbf_log_total_time_sec = 0.0;
    cbf_log_avg_time_sec = 0.0;
    cbf_log_max_time_sec = 0.0;
    cbf_log_call_count = 0;
    bcbf_cd_initialized = 0;
}

// =============================================================================
// Public filter function. Same signature as the original CBF filter.
// Inputs *vx_cmd,*vy_cmd,*wz_cmd are interpreted as nominal u_d.
// Outputs are the integrated actual command u, not uhat.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    const double filter_t0 = cbf_now_sec();
    double ud[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    double W[3]  = {wx, wy, ww};
    for (int j = 0; j < 3; j++) if (W[j] <= 1e-9) W[j] = 1.0;

    cbf_log_u_nom[0] = ud[0];
    cbf_log_u_nom[1] = ud[1];
    cbf_log_u_nom[2] = ud[2];

    double ub_log[3];
    bcbf_backup_controller(x, y, theta, ub_log);
    cbf_log_u_backup[0] = ub_log[0];
    cbf_log_u_backup[1] = ub_log[1];
    cbf_log_u_backup[2] = ub_log[2];

    if (!bcbf_cd_initialized) {
        // Initialize actual input state inside admissible input set.
        for (int j = 0; j < 3; j++) {
            bcbf_cd_u_state[j] = bcbf_clip(ud[j], bcbf_u_min[j], bcbf_u_max[j]);
            bcbf_cd_ud_prev[j] = ud[j];
        }
        bcbf_cd_initialized = 1;
    }

    // Desired surrogate input uhat_d that would make u track ud.
    double ud_dot[3], uhat_d[3];
    for (int j = 0; j < 3; j++) {
        ud_dot[j] = (ud[j] - bcbf_cd_ud_prev[j]) / fmax(bcbf_cd_dt, 1e-6);
        uhat_d[j] = bcbf_cd_u_state[j]
                    + (1.0 / bcbf_cd_a[j]) * ud_dot[j]
                    + (bcbf_sigma0[j] / bcbf_cd_a[j]) * (ud[j] - bcbf_cd_u_state[j]);
        cbf_log_uhat_d[j] = uhat_d[j];
    }

    // Compute composite augmented barrier and its analytical gradient in one
    // augmented backup rollout using Phi_B, matching the CD QP baselines.
    double h_aug;
    double grad[6];
    bcbf_augmented_barrier_value_and_grad(
        x, y, theta, bcbf_cd_u_state, &h_aug, grad);
    double h_backup_dbg = cbf_log_h_backup_lse;
    double h_input_dbg  = cbf_log_h_input_lse;

    // Augmented drift fhat = [f_old(x,u); -A u].
    double fhat[6];
    bcbf_fhat_drift(x, y, theta, bcbf_cd_u_state, fhat);

    double Lf = 0.0;
    for (int j = 0; j < 6; j++) Lf += grad[j] * fhat[j];

    // Lg = grad_u * A because ghat = [0; A].
    double Lg[3];
    for (int j = 0; j < 3; j++) Lg[j] = grad[3+j] * bcbf_cd_a[j];

    // Closed-form KKT solution for:
    // min 0.5(uhat-uhat_d)'W(uhat-uhat_d) + 0.5 gamma mu^2
    // s.t. Lf + Lg uhat + alpha h + mu h >= 0.
    double omega = Lf + bcbf_alpha_outer * h_aug;
    for (int j = 0; j < 3; j++) omega += Lg[j] * uhat_d[j];

    double denom = (h_aug * h_aug) / bcbf_gamma_mu;
    for (int j = 0; j < 3; j++) denom += (Lg[j] * Lg[j]) / W[j];

    // Timing starts here to match solver-only timing used by QP/OI baselines.
    // This excludes constraint/barrier construction, rollout, and gradient computation.
    double qp_t0 = cbf_now_sec();

    double lambda = 0.0;
    if (omega < 0.0 && denom > 1e-12 && isfinite(denom)) {
        lambda = -omega / denom;
    }

    double uhat[3];
    for (int j = 0; j < 3; j++) {
        uhat[j] = uhat_d[j] + (Lg[j] / W[j]) * lambda;
        cbf_log_uhat[j] = uhat[j];
    }

    double mu = (h_aug * lambda) / bcbf_gamma_mu;

    cbf_log_qp_time_sec = cbf_now_sec() - qp_t0;

    // Debug print: keep while validating, remove or throttle later.
//     printf("[CD-bCBF-AUG%s] x=%.3f y=%.3f th=%.3f | "
//            "h=%.5f hB=%.5f hU=%.5f omega=%.5f den=%.3e lam=%.5f | "
//            "ud=[%.3f %.3f %.3f] uhatd=[%.3f %.3f %.3f] "
//            "uhat=[%.3f %.3f %.3f] u=[%.3f %.3f %.3f]\n",
// #if BCBF_USE_AUGMENTED_TERMINAL_SET
//            "-XHAT",
// #else
//            "-X",
// #endif
//            x, y, theta, h_aug, h_backup_dbg, h_input_dbg, omega, denom, lambda,
//            ud[0], ud[1], ud[2],
//            uhat_d[0], uhat_d[1], uhat_d[2],
//            uhat[0], uhat[1], uhat[2],
//            bcbf_cd_u_state[0], bcbf_cd_u_state[1], bcbf_cd_u_state[2]);

    // Integrate actual input state u_dot = -A u + A uhat.
    for (int j = 0; j < 3; j++) {
        bcbf_cd_u_state[j] += bcbf_cd_dt * (-bcbf_cd_a[j] * bcbf_cd_u_state[j]
                                            + bcbf_cd_a[j] * uhat[j]);
#if BCBF_CD_HARD_CLIP_OUTPUT
        bcbf_cd_u_state[j] = bcbf_clip(bcbf_cd_u_state[j], bcbf_u_min[j], bcbf_u_max[j]);
#endif
    }

    *vx_cmd = bcbf_cd_u_state[0];
    *vy_cmd = bcbf_cd_u_state[1];
    *wz_cmd = bcbf_cd_u_state[2];

    cbf_log_u_safe[0] = bcbf_cd_u_state[0];
    cbf_log_u_safe[1] = bcbf_cd_u_state[1];
    cbf_log_u_safe[2] = bcbf_cd_u_state[2];
    cbf_log_mu = mu;
    cbf_log_lambda = lambda;
    cbf_log_h_aug = h_aug;
    cbf_log_omega = omega;
    cbf_log_h_backup_lse = h_backup_dbg;
    cbf_log_h_input_lse = h_input_dbg;
    cbf_log_status = 1;

    for (int j = 0; j < 3; j++) bcbf_cd_ud_prev[j] = ud[j];

    if (isfinite(h_aug) && isfinite(lambda)) cf_success_count++;
    else cf_fail_count++;

    cbf_log_filter_time_sec = cbf_now_sec() - filter_t0;
    cbf_log_call_count++;
    cbf_log_total_time_sec += cbf_log_filter_time_sec;
    cbf_log_avg_time_sec = cbf_log_total_time_sec / (double)cbf_log_call_count;
    if (cbf_log_filter_time_sec > cbf_log_max_time_sec) {
        cbf_log_max_time_sec = cbf_log_filter_time_sec;
    }

    return 1;
}
