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
static double cbf_log_qp_time_sec = 0.0;
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

#define BCBF_N_OBS 4
static BCBFObstacle bcbf_obs[BCBF_N_OBS] = {
    { -2.0,  1.0, 0.35, 0.5, 1 },
    { -4.0, -1.0, 0.35, 0.5, 1 },
    { -2.0, -1.0, 0.35, 0.5, 1 },
    { -4.0,  1.0, 0.35, 0.5, 1 }
};

// Actual input limits for u = [vx, vy, omega].
static double bcbf_u_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_min[3] = { -1.0, -0.3, -1.0 };

// Backup set: enlarged-obstacle clearance set.
static const double bcbf_backup_margin = 0.15;

// Terminal near-stop set used only when BCBF_USE_AUGMENTED_TERMINAL_SET = 1.
static const double bcbf_terminal_u_eps[3] = {0.08, 0.04, 0.08};

// LSE gains. Use hierarchical LSE to avoid one giant 86-term soft-min dominating.
static const double bcbf_kappa_backup_lse = 50.0;
static const double bcbf_kappa_input_lse  = 50.0;
static const double bcbf_kappa_outer_lse  = 30.0;

// HOCBF/ECBF and outer R-CBF gains.
// For each backup rollout term H_i(xhat), construct
//      b_i = L_fhat H_i + alpha0 H_i,
// then LSE combines the b_i terms.
static const double bcbf_alpha0_hocbf = 0.5;
static const double bcbf_alpha_outer  = 0.5;

// Closed-form slack penalty. Larger gamma discourages mu.
static const double bcbf_gamma_mu = 1.0e4;

// Control dynamics: u_dot = -A u + A uhat.
static double bcbf_cd_a[3] = { 8.0, 8.0, 8.0 };

// Input tracking decay for u -> u_d through uhat_d.
static double bcbf_sigma0[3] = { 4.0, 4.0, 4.0 };

// Control-loop time step used to integrate u_dot and finite-difference u_d_dot.
static double bcbf_cd_dt = 0.02;

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

// Roll out the AUGMENTED zero-backup dynamics and return scalar backup terms:
//   - obstacle safety along the horizon,
//   - terminal enlarged-obstacle backup set,
//   - optionally terminal near-stop condition.
static int bcbf_aug_backup_terms(double x, double y, double theta,
                                 const double u0[3],
                                 double H[BCBF_N_BACKUP_TERMS])
{
    double z[6] = {x, y, theta, u0[0], u0[1], u0[2]};
    double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);
    int term = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        // Obstacle margins at augmented backup-flow sample.
        for (int i = 0; i < BCBF_N_OBS; i++) {
            H[term++] = bcbf_h_obs(i, z[0], z[1]);
        }

        if (k < BCBF_N_SAMPLES - 1) {
            double c = cos(z[2]);
            double s = sin(z[2]);

            double zdot[6];
            zdot[0] = c * z[3] - s * z[4];
            zdot[1] = s * z[3] + c * z[4];
            zdot[2] = z[5];

            // Zero surrogate backup: uhat_B = 0, so u_dot = -A u.
            zdot[3] = -bcbf_cd_a[0] * z[3];
            zdot[4] = -bcbf_cd_a[1] * z[4];
            zdot[5] = -bcbf_cd_a[2] * z[5];

            for (int j = 0; j < 6; j++) z[j] += dt * zdot[j];
            z[2] = bcbf_wrap_pi(z[2]);
        }
    }

    // Terminal backup set based on original state clearance.
    for (int i = 0; i < BCBF_N_OBS; i++) {
        H[term++] = bcbf_h_backup_obs(i, z[0], z[1]);
    }

#if BCBF_USE_AUGMENTED_TERMINAL_SET
    // Extra augmented terminal condition: command nearly stopped.
    double uT[3] = {z[3], z[4], z[5]};
    H[term++] = bcbf_h_terminal_stop(uT);
#endif

    return term;
}

// Composite backup barrier over HOCBF-prepared augmented backup rollout terms.
// For every H_i(xhat), compute b_i = L_fhat H_i + alpha0 H_i via directional
// finite difference along the current drift fhat. Then LSE over b_i.
static double bcbf_backup_hocbf_lse_value(double x, double y, double theta,
                                          const double u[3])
{
    double H0[BCBF_N_BACKUP_TERMS];
    double H1[BCBF_N_BACKUP_TERMS];

    int n0 = bcbf_aug_backup_terms(x, y, theta, u, H0);

    double fhat[6];
    bcbf_fhat_drift(x, y, theta, u, fhat);

    double xp = x     + bcbf_dir_eps * fhat[0];
    double yp = y     + bcbf_dir_eps * fhat[1];
    double tp = theta + bcbf_dir_eps * fhat[2];
    double up[3] = {
        u[0] + bcbf_dir_eps * fhat[3],
        u[1] + bcbf_dir_eps * fhat[4],
        u[2] + bcbf_dir_eps * fhat[5]
    };
    tp = bcbf_wrap_pi(tp);

    int n1 = bcbf_aug_backup_terms(xp, yp, tp, up, H1);
    int n = (n0 < n1) ? n0 : n1;

    double b_terms[BCBF_N_BACKUP_TERMS];
    for (int i = 0; i < n; i++) {
        double Hdot = (H1[i] - H0[i]) / bcbf_dir_eps;
        b_terms[i] = Hdot + bcbf_alpha0_hocbf * H0[i];
    }

    return bcbf_softmin_from_terms(b_terms, n, bcbf_kappa_backup_lse);
}

static double bcbf_input_lse_value(const double u[3])
{
    double terms[BCBF_N_INPUT_TERMS];
    int n = 0;
    terms[n++] = bcbf_u_max[0] - u[0];
    terms[n++] = u[0] - bcbf_u_min[0];
    terms[n++] = bcbf_u_max[1] - u[1];
    terms[n++] = u[1] - bcbf_u_min[1];
    terms[n++] = bcbf_u_max[2] - u[2];
    terms[n++] = u[2] - bcbf_u_min[2];
    return bcbf_softmin_from_terms(terms, n, bcbf_kappa_input_lse);
}

// Final hierarchical composite augmented barrier.
static double bcbf_augmented_barrier_value(double x, double y, double theta,
                                           const double u[3])
{
    double h_backup = bcbf_backup_hocbf_lse_value(x, y, theta, u);
    double h_input  = bcbf_input_lse_value(u);

    cbf_log_h_backup_lse = h_backup;
    cbf_log_h_input_lse  = h_input;

    double outer_terms[2] = {h_backup, h_input};
    return bcbf_softmin_from_terms(outer_terms, 2, bcbf_kappa_outer_lse);
}

// Finite-difference gradient of h_aug wrt xhat = [x,y,theta,vx,vy,w].
static void bcbf_augmented_barrier_grad_fd(double x, double y, double theta,
                                           const double u[3], double grad[6])
{
    const double eps[6] = {1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 1e-5};

    for (int j = 0; j < 6; j++) {
        double xp = x, yp = y, thp = theta, up[3] = {u[0], u[1], u[2]};
        double xm = x, ym = y, thm = theta, um[3] = {u[0], u[1], u[2]};

        if (j == 0) { xp += eps[j]; xm -= eps[j]; }
        else if (j == 1) { yp += eps[j]; ym -= eps[j]; }
        else if (j == 2) { thp = bcbf_wrap_pi(thp + eps[j]); thm = bcbf_wrap_pi(thm - eps[j]); }
        else { up[j-3] += eps[j]; um[j-3] -= eps[j]; }

        double hp = bcbf_augmented_barrier_value(xp, yp, thp, up);
        double hm = bcbf_augmented_barrier_value(xm, ym, thm, um);
        grad[j] = (hp - hm) / (2.0 * eps[j]);
    }
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
    double t0 = cbf_now_sec();

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

    // Compute composite augmented barrier and gradient.
    double h_aug = bcbf_augmented_barrier_value(x, y, theta, bcbf_cd_u_state);
    double h_backup_dbg = cbf_log_h_backup_lse;
    double h_input_dbg  = cbf_log_h_input_lse;

    double grad[6];
    bcbf_augmented_barrier_grad_fd(x, y, theta, bcbf_cd_u_state, grad);

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
    cbf_log_qp_time_sec = cbf_now_sec() - t0;

    for (int j = 0; j < 3; j++) bcbf_cd_ud_prev[j] = ud[j];

    if (isfinite(h_aug) && isfinite(lambda)) cf_success_count++;
    else cf_fail_count++;

    return 1;
}
