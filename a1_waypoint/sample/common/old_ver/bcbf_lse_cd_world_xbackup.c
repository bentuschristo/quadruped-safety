#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>
#include <time.h>

// =============================================================================
// Closed-form LSE Backup-CBF filter with WORLD-FRAME control dynamics.
//
// Interface remains the same as the old filter:
//   input/output commands are BODY-FRAME [vx, vy, omega].
//
// Internally, the control-dynamics state is WORLD-FRAME:
//   nu = [vX, vY, omega],
//   [vX; vY] = R(theta) [vx; vy].
//
// Control dynamics:
//   nu_dot = -A nu + A nuhat.
//
// Plant:
//   x_dot     = vX
//   y_dot     = vY
//   theta_dot = omega.
//
// Output back to the robot:
//   [vx; vy] = R(theta)^T [vX; vY].
//
// Variant controlled by BCBF_BACKUP_FLOW_MODE:
//   0: backup flow uses original state x only.
//      This is the simpler/older-style backup flow. Zero backup means the
//      backup position is fixed during the backup rollout.
//   1: backup flow uses augmented state xhat=[x,y,theta,vX,vY,omega].
//      Zero backup means nu decays by nu_dot=-A nu during the backup rollout.
//
// Variant controlled by BCBF_USE_AUGMENTED_TERMINAL_SET:
//   0: terminal backup set checks only original-state enlarged obstacle clearance.
//   1: terminal backup set also checks near-stopped terminal command.
// =============================================================================

#ifndef BCBF_BACKUP_FLOW_MODE
#define BCBF_BACKUP_FLOW_MODE 0
#endif

#ifndef BCBF_USE_AUGMENTED_TERMINAL_SET
#define BCBF_USE_AUGMENTED_TERMINAL_SET 0
#endif

static int cf_success_count = 0;
static int cf_fail_count = 0;

// Logging compatible with previous files.
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};  // body-frame nominal
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};  // body-frame backup log
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};  // body-frame output
static double cbf_log_uhat_d[3]   = {0.0, 0.0, 0.0};  // world-frame desired surrogate
static double cbf_log_uhat[3]     = {0.0, 0.0, 0.0};  // world-frame filtered surrogate
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
    { -2.0,  1.0, 0.45, 0.5, 1 },
    { -4.0, -1.0, 0.45, 0.5, 1 },
    { -2.0, -1.0, 0.45, 0.5, 1 },
    { -4.0,  1.0, 0.45, 0.5, 1 }
};

// Actual ROBOT BODY-FRAME command limits.
// Even though the internal control-dynamics state is world-frame nu, the command
// sent to the robot is body-frame, so these barriers are applied to R(theta)^T nu.
static double bcbf_u_body_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_body_min[3] = { -1.0, -0.3, -1.0 };

// Backup set: enlarged-obstacle clearance set.
static const double bcbf_backup_margin = 0.15;

// Terminal near-stop set used only when BCBF_USE_AUGMENTED_TERMINAL_SET = 1.
static const double bcbf_terminal_nu_eps[3] = {0.08, 0.04, 0.08};

// LSE gains. Hierarchical LSE keeps the input-barrier soft-min from swallowing
// the whole backup certificate.
static const double bcbf_kappa_backup_lse = 50.0;
static const double bcbf_kappa_input_lse  = 50.0;
static const double bcbf_kappa_outer_lse  = 30.0;

// HOCBF/ECBF and outer R-CBF gains.
static const double bcbf_alpha0_hocbf = 1.0;
static const double bcbf_alpha_outer  = 0.5;

// Closed-form slack penalty.
static const double bcbf_gamma_mu = 1.0e4;

// WORLD-FRAME control dynamics: nu_dot = -A nu + A nuhat.
static double bcbf_cd_a[3] = { 8.0, 8.0, 8.0 };

// Tracking decay for nu -> nu_d through nuhat_d.
static double bcbf_sigma0[3] = { 4.0, 4.0, 4.0 };

// Control-loop time step used to integrate nu_dot and finite-difference nu_d_dot.
static double bcbf_cd_dt = 0.02;

// Directional derivative step used to compute L_fhat H_i.
static const double bcbf_dir_eps = 1.0e-5;

// Final numerical guard. This clips the final BODY-FRAME output command.
#define BCBF_CD_HARD_CLIP_OUTPUT 1

#define BCBF_N_SAMPLES 20
static const double bcbf_T_backup = 4.0;

#define BCBF_N_HORIZON_OBS_TERMS (BCBF_N_OBS * BCBF_N_SAMPLES)
#define BCBF_N_TERMINAL_CLEARANCE_TERMS BCBF_N_OBS

#if BCBF_USE_AUGMENTED_TERMINAL_SET
#define BCBF_N_TERMINAL_STOP_TERMS 1
#else
#define BCBF_N_TERMINAL_STOP_TERMS 0
#endif

#define BCBF_N_BACKUP_TERMS (BCBF_N_HORIZON_OBS_TERMS + BCBF_N_TERMINAL_CLEARANCE_TERMS + BCBF_N_TERMINAL_STOP_TERMS)
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

static inline void bcbf_body_to_world(double theta, const double u_body[3], double nu_world[3])
{
    double c = cos(theta);
    double s = sin(theta);
    nu_world[0] = c * u_body[0] - s * u_body[1];
    nu_world[1] = s * u_body[0] + c * u_body[1];
    nu_world[2] = u_body[2];
}

static inline void bcbf_world_to_body(double theta, const double nu_world[3], double u_body[3])
{
    double c = cos(theta);
    double s = sin(theta);
    u_body[0] =  c * nu_world[0] + s * nu_world[1];
    u_body[1] = -s * nu_world[0] + c * nu_world[1];
    u_body[2] =  nu_world[2];
}

static inline void bcbf_fhat_drift_world(double x, double y, double theta,
                                         const double nu[3], double fhat[6])
{
    (void)x; (void)y;
    fhat[0] = nu[0];
    fhat[1] = nu[1];
    fhat[2] = nu[2];
    fhat[3] = -bcbf_cd_a[0] * nu[0];
    fhat[4] = -bcbf_cd_a[1] * nu[1];
    fhat[5] = -bcbf_cd_a[2] * nu[2];
}

static inline void bcbf_backup_controller_body(double x, double y, double th, double ub[3])
{
    // Zero backup for logging. In world-frame backup, zero body command is
    // equivalent to zero world velocity command.
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

static inline double bcbf_h_terminal_stop(const double nuT[3])
{
    double sx = nuT[0] / fmax(bcbf_terminal_nu_eps[0], 1e-9);
    double sy = nuT[1] / fmax(bcbf_terminal_nu_eps[1], 1e-9);
    double sw = nuT[2] / fmax(bcbf_terminal_nu_eps[2], 1e-9);
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

// Backup terms for WORLD-FRAME control dynamics.
// Mode 0: x-only backup flow: zero backup freezes x,y,theta.
// Mode 1: xhat backup flow: zero backup means nu decays while x,y,theta keep moving.
static int bcbf_backup_terms_world(double x, double y, double theta,
                                   const double nu0[3],
                                   double H[BCBF_N_BACKUP_TERMS])
{
    double z[6] = {x, y, theta, nu0[0], nu0[1], nu0[2]};
    double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);
    int term = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        for (int i = 0; i < BCBF_N_OBS; i++) {
            H[term++] = bcbf_h_obs(i, z[0], z[1]);
        }

        if (k < BCBF_N_SAMPLES - 1) {
#if BCBF_BACKUP_FLOW_MODE == 1
            // Augmented zero-backup rollout:
            // nuhat_B=0, so nu_dot=-A nu; x_dot=nu_x, y_dot=nu_y.
            double zdot[6];
            zdot[0] = z[3];
            zdot[1] = z[4];
            zdot[2] = z[5];
            zdot[3] = -bcbf_cd_a[0] * z[3];
            zdot[4] = -bcbf_cd_a[1] * z[4];
            zdot[5] = -bcbf_cd_a[2] * z[5];

            for (int j = 0; j < 6; j++) z[j] += dt * zdot[j];
            z[2] = bcbf_wrap_pi(z[2]);
#else
            // Original-state x-only zero-backup rollout:
            // zero actual command freezes the backup position.
            // Keep z unchanged.
            (void)dt;
#endif
        }
    }

    for (int i = 0; i < BCBF_N_OBS; i++) {
        H[term++] = bcbf_h_backup_obs(i, z[0], z[1]);
    }

#if BCBF_USE_AUGMENTED_TERMINAL_SET
    double nuT[3] = {z[3], z[4], z[5]};
    H[term++] = bcbf_h_terminal_stop(nuT);
#endif

    return term;
}

// HOCBF-prepared backup LSE: b_i = L_fhat H_i + alpha0 H_i.
// L_fhat H_i is computed by directional finite difference along the current
// world-frame augmented drift.
static double bcbf_backup_hocbf_lse_value_world(double x, double y, double theta,
                                                const double nu[3])
{
    double H0[BCBF_N_BACKUP_TERMS];
    double H1[BCBF_N_BACKUP_TERMS];

    int n0 = bcbf_backup_terms_world(x, y, theta, nu, H0);

    double fhat[6];
    bcbf_fhat_drift_world(x, y, theta, nu, fhat);

    double xp = x     + bcbf_dir_eps * fhat[0];
    double yp = y     + bcbf_dir_eps * fhat[1];
    double tp = theta + bcbf_dir_eps * fhat[2];
    double nup[3] = {
        nu[0] + bcbf_dir_eps * fhat[3],
        nu[1] + bcbf_dir_eps * fhat[4],
        nu[2] + bcbf_dir_eps * fhat[5]
    };
    tp = bcbf_wrap_pi(tp);

    int n1 = bcbf_backup_terms_world(xp, yp, tp, nup, H1);
    int n = (n0 < n1) ? n0 : n1;

    double b_terms[BCBF_N_BACKUP_TERMS];
    for (int i = 0; i < n; i++) {
        double Hdot = (H1[i] - H0[i]) / bcbf_dir_eps;
        b_terms[i] = Hdot + bcbf_alpha0_hocbf * H0[i];
    }

    return bcbf_softmin_from_terms(b_terms, n, bcbf_kappa_backup_lse);
}

// Input barriers are on the BODY-FRAME commands sent to the robot:
//   u_body = R(theta)^T nu_world.
static double bcbf_input_lse_value_world(double theta, const double nu[3])
{
    double ub[3];
    bcbf_world_to_body(theta, nu, ub);

    double terms[BCBF_N_INPUT_TERMS];
    int n = 0;
    terms[n++] = bcbf_u_body_max[0] - ub[0];
    terms[n++] = ub[0] - bcbf_u_body_min[0];
    terms[n++] = bcbf_u_body_max[1] - ub[1];
    terms[n++] = ub[1] - bcbf_u_body_min[1];
    terms[n++] = bcbf_u_body_max[2] - ub[2];
    terms[n++] = ub[2] - bcbf_u_body_min[2];
    return bcbf_softmin_from_terms(terms, n, bcbf_kappa_input_lse);
}

static double bcbf_augmented_barrier_value_world(double x, double y, double theta,
                                                 const double nu[3])
{
    double h_backup = bcbf_backup_hocbf_lse_value_world(x, y, theta, nu);
    double h_input  = bcbf_input_lse_value_world(theta, nu);

    cbf_log_h_backup_lse = h_backup;
    cbf_log_h_input_lse  = h_input;

    double outer_terms[2] = {h_backup, h_input};
    return bcbf_softmin_from_terms(outer_terms, 2, bcbf_kappa_outer_lse);
}

// Finite-difference gradient of h_aug wrt xhat = [x,y,theta,vX,vY,w].
static void bcbf_augmented_barrier_grad_fd_world(double x, double y, double theta,
                                                 const double nu[3], double grad[6])
{
    const double eps[6] = {1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 1e-5};

    for (int j = 0; j < 6; j++) {
        double xp = x, yp = y, thp = theta, nup[3] = {nu[0], nu[1], nu[2]};
        double xm = x, ym = y, thm = theta, num[3] = {nu[0], nu[1], nu[2]};

        if (j == 0) { xp += eps[j]; xm -= eps[j]; }
        else if (j == 1) { yp += eps[j]; ym -= eps[j]; }
        else if (j == 2) { thp = bcbf_wrap_pi(thp + eps[j]); thm = bcbf_wrap_pi(thm - eps[j]); }
        else { nup[j-3] += eps[j]; num[j-3] -= eps[j]; }

        double hp = bcbf_augmented_barrier_value_world(xp, yp, thp, nup);
        double hm = bcbf_augmented_barrier_value_world(xm, ym, thm, num);
        grad[j] = (hp - hm) / (2.0 * eps[j]);
    }
}

// Persistent world-frame velocity command state nu.
int bcbf_cd_initialized = 0;
static double bcbf_cd_nu_state[3]  = {0.0, 0.0, 0.0};
static double bcbf_cd_nu_d_prev[3] = {0.0, 0.0, 0.0};

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
// Inputs *vx_cmd,*vy_cmd,*wz_cmd are BODY-FRAME nominal commands.
// Internally converts them to WORLD-FRAME nu_d, filters nu, then converts back.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    double t0 = cbf_now_sec();

    double u_body_nom[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    double nu_d[3];
    bcbf_body_to_world(theta, u_body_nom, nu_d);

    // Keep user-provided weights, interpreted in world-frame surrogate coordinates.
    double W[3] = {wx, wy, ww};
    for (int j = 0; j < 3; j++) if (W[j] <= 1e-9) W[j] = 1.0;

    cbf_log_u_nom[0] = u_body_nom[0];
    cbf_log_u_nom[1] = u_body_nom[1];
    cbf_log_u_nom[2] = u_body_nom[2];

    double ub_log[3];
    bcbf_backup_controller_body(x, y, theta, ub_log);
    cbf_log_u_backup[0] = ub_log[0];
    cbf_log_u_backup[1] = ub_log[1];
    cbf_log_u_backup[2] = ub_log[2];

    if (!bcbf_cd_initialized) {
        // Initialize world-frame command state from nominal body command.
        for (int j = 0; j < 3; j++) {
            bcbf_cd_nu_state[j] = nu_d[j];
            bcbf_cd_nu_d_prev[j] = nu_d[j];
        }

        // Ensure initialized body output respects limits.
        double u_init_body[3];
        bcbf_world_to_body(theta, bcbf_cd_nu_state, u_init_body);
        u_init_body[0] = bcbf_clip(u_init_body[0], bcbf_u_body_min[0], bcbf_u_body_max[0]);
        u_init_body[1] = bcbf_clip(u_init_body[1], bcbf_u_body_min[1], bcbf_u_body_max[1]);
        u_init_body[2] = bcbf_clip(u_init_body[2], bcbf_u_body_min[2], bcbf_u_body_max[2]);
        bcbf_body_to_world(theta, u_init_body, bcbf_cd_nu_state);

        bcbf_cd_initialized = 1;
    }

    // Desired surrogate input nuhat_d that would make nu track nu_d.
    double nu_d_dot[3], nuhat_d[3];
    for (int j = 0; j < 3; j++) {
        nu_d_dot[j] = (nu_d[j] - bcbf_cd_nu_d_prev[j]) / fmax(bcbf_cd_dt, 1e-6);
        nuhat_d[j] = bcbf_cd_nu_state[j]
                     + (1.0 / bcbf_cd_a[j]) * nu_d_dot[j]
                     + (bcbf_sigma0[j] / bcbf_cd_a[j]) * (nu_d[j] - bcbf_cd_nu_state[j]);
        cbf_log_uhat_d[j] = nuhat_d[j];
    }

    double h_aug = bcbf_augmented_barrier_value_world(x, y, theta, bcbf_cd_nu_state);
    double h_backup_dbg = cbf_log_h_backup_lse;
    double h_input_dbg  = cbf_log_h_input_lse;

    double grad[6];
    bcbf_augmented_barrier_grad_fd_world(x, y, theta, bcbf_cd_nu_state, grad);

    double fhat[6];
    bcbf_fhat_drift_world(x, y, theta, bcbf_cd_nu_state, fhat);

    double Lf = 0.0;
    for (int j = 0; j < 6; j++) Lf += grad[j] * fhat[j];

    // ghat = [0; A], so Lg = grad_nu * A.
    double Lg[3];
    for (int j = 0; j < 3; j++) Lg[j] = grad[3+j] * bcbf_cd_a[j];

    double omega = Lf + bcbf_alpha_outer * h_aug;
    for (int j = 0; j < 3; j++) omega += Lg[j] * nuhat_d[j];

    double denom = (h_aug * h_aug) / bcbf_gamma_mu;
    for (int j = 0; j < 3; j++) denom += (Lg[j] * Lg[j]) / W[j];

    double lambda = 0.0;
    if (omega < 0.0 && denom > 1e-12 && isfinite(denom)) {
        lambda = -omega / denom;
    }

    double nuhat[3];
    for (int j = 0; j < 3; j++) {
        nuhat[j] = nuhat_d[j] + (Lg[j] / W[j]) * lambda;
        cbf_log_uhat[j] = nuhat[j];
    }

    double mu = (h_aug * lambda) / bcbf_gamma_mu;

    printf("[CD-bCBF-WORLD-%s%s] x=%.3f y=%.3f th=%.3f | "
           "h=%.5f hB=%.5f hU=%.5f omega=%.5f den=%.3e lam=%.5f | "
           "u_nom_body=[%.3f %.3f %.3f] nu_d=[%.3f %.3f %.3f] "
           "nuhat=[%.3f %.3f %.3f] nu=[%.3f %.3f %.3f]\n",
#if BCBF_BACKUP_FLOW_MODE == 1
           "XHAT",
#else
           "X",
#endif
#if BCBF_USE_AUGMENTED_TERMINAL_SET
           "-TERM",
#else
           "",
#endif
           x, y, theta, h_aug, h_backup_dbg, h_input_dbg, omega, denom, lambda,
           u_body_nom[0], u_body_nom[1], u_body_nom[2],
           nu_d[0], nu_d[1], nu_d[2],
           nuhat[0], nuhat[1], nuhat[2],
           bcbf_cd_nu_state[0], bcbf_cd_nu_state[1], bcbf_cd_nu_state[2]);

    // Integrate world-frame command state nu_dot = -A nu + A nuhat.
    for (int j = 0; j < 3; j++) {
        bcbf_cd_nu_state[j] += bcbf_cd_dt * (-bcbf_cd_a[j] * bcbf_cd_nu_state[j]
                                             + bcbf_cd_a[j] * nuhat[j]);
    }

    // Convert filtered world-frame command back to body frame.
    double u_body_safe[3];
    bcbf_world_to_body(theta, bcbf_cd_nu_state, u_body_safe);

#if BCBF_CD_HARD_CLIP_OUTPUT
    u_body_safe[0] = bcbf_clip(u_body_safe[0], bcbf_u_body_min[0], bcbf_u_body_max[0]);
    u_body_safe[1] = bcbf_clip(u_body_safe[1], bcbf_u_body_min[1], bcbf_u_body_max[1]);
    u_body_safe[2] = bcbf_clip(u_body_safe[2], bcbf_u_body_min[2], bcbf_u_body_max[2]);

    // After clipping body-frame output, reset internal world state to match
    // what was actually sent. This avoids internal/external command mismatch.
    bcbf_body_to_world(theta, u_body_safe, bcbf_cd_nu_state);
#endif

    *vx_cmd = u_body_safe[0];
    *vy_cmd = u_body_safe[1];
    *wz_cmd = u_body_safe[2];

    cbf_log_u_safe[0] = u_body_safe[0];
    cbf_log_u_safe[1] = u_body_safe[1];
    cbf_log_u_safe[2] = u_body_safe[2];

    cbf_log_mu = mu;
    cbf_log_lambda = lambda;
    cbf_log_h_aug = h_aug;
    cbf_log_omega = omega;
    cbf_log_h_backup_lse = h_backup_dbg;
    cbf_log_h_input_lse = h_input_dbg;
    cbf_log_status = 1;
    cbf_log_qp_time_sec = cbf_now_sec() - t0;

    for (int j = 0; j < 3; j++) bcbf_cd_nu_d_prev[j] = nu_d[j];

    if (isfinite(h_aug) && isfinite(lambda)) cf_success_count++;
    else cf_fail_count++;

    return 1;
}
