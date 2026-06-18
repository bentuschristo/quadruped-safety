#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>
#include <time.h>

// =============================================================================
// Closed-form LSE Backup-CBF filter with control dynamics for high-level quadruped
//
// Original plant state:        x = [px, py, theta]
// Actual plant input/state:    u = [vx, vy, omega]
// Surrogate decision input:    uhat = [hat_vx, hat_vy, hat_omega]
// Augmented state:             xhat = [px, py, theta, vx, vy, omega]
//
// Plant dynamics:
//      px_dot    = cos(theta) vx - sin(theta) vy
//      py_dot    = sin(theta) vx + cos(theta) vy
//      theta_dot = omega
//
// Control dynamics:
//      u_dot = -A u + A uhat,   A = diag(a_vx, a_vy, a_w)
//
// This file removes the OSQP call. The filter computes uhat from the closed-form
// one-constraint CBF solution, then integrates u_dot to produce the bounded
// command u sent to the robot.
//
// Based on the user's bcbf_lse_circle_obstacles_v2.c zero-action backup setup.
// =============================================================================

static int cf_success_count = 0;
static int cf_fail_count = 0;

// Logging compatible with previous files.
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};  // nominal u_d
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};  // backup command
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};  // actual integrated u
static double cbf_log_uhat_d[3]   = {0.0, 0.0, 0.0};  // desired surrogate input
static double cbf_log_uhat[3]     = {0.0, 0.0, 0.0};  // filtered surrogate input
static double cbf_log_mu          = 0.0;
static double cbf_log_lambda      = 0.0;
static int    cbf_log_status      = 0;
static double cbf_log_qp_time_sec = 0.0;              // closed-form compute time
static double cbf_log_h_aug       = 0.0;
static double cbf_log_omega       = 0.0;

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

// Actual input limits for u = [vx, vy, omega]. These are enforced as augmented
// state barriers, not as direct QP box constraints.
static double bcbf_u_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_min[3] = { -1.0, -0.3, -1.0 };

// Backup set, v2: enlarged-obstacle clearance set.
static const double bcbf_backup_margin = 0.15;
static const double bcbf_kappa_lse_B = 50.0;
static const double bcbf_kappa_lse   = 50.0;

// HOCBF/ECBF and outer R-CBF gains.
// For each backup rollout term H_i(x), construct b_i = Hdot_i + alpha0 H_i.
// Then LSE combines the b_i and the input-bound barriers.
static const double bcbf_alpha0_hocbf = 1.0;
static const double bcbf_alpha_outer  = 0.5;

// Closed-form slack penalty. Larger gamma discourages mu.
static const double bcbf_gamma_mu = 1.0e4;

// Control dynamics: u_dot = -A u + A uhat.
static double bcbf_cd_a[3] = { 8.0, 8.0, 8.0 };

// Input tracking decay for u -> u_d through uhat_d.
static double bcbf_sigma0[3] = { 4.0, 4.0, 4.0 };

// Control-loop time step used to integrate u_dot and finite-difference u_d_dot.
// Set this to your actual loop period, e.g., 0.02 for 50 Hz or 0.01 for 100 Hz.
static double bcbf_cd_dt = 0.02;

// Final numerical guard. Theory should enforce bounds through barriers; this clip
// is only a hardware safety guard against finite-difference/discretization error.
#define BCBF_CD_HARD_CLIP_OUTPUT 1

#define BCBF_N_SAMPLES 20
static const double bcbf_T_backup = 4.0;
#define BCBF_N_BACKUP_TERMS (BCBF_N_OBS * BCBF_N_SAMPLES + BCBF_N_OBS)
#define BCBF_N_INPUT_TERMS  6
#define BCBF_N_AUG_TERMS    (BCBF_N_BACKUP_TERMS + BCBF_N_INPUT_TERMS)

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

static inline void bcbf_backup_controller(double x, double y, double th, double ub[3])
{
    // v2 zero-action backup.
    (void)x; (void)y; (void)th;
    ub[0] = 0.0;
    ub[1] = 0.0;
    ub[2] = 0.0;
}

static inline void bcbf_closed_loop(double z[3], double zdot[3])
{
    double ub[3];
    bcbf_backup_controller(z[0], z[1], z[2], ub);
    bcbf_dynamics(z[0], z[1], z[2], ub, zdot);
}

static void bcbf_cl_jacobian_fd(const double z[3], double J[3][3])
{
    const double eps = 1e-5;
    double zp[3], zm[3], fp[3], fm[3];

    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) { zp[k] = z[k]; zm[k] = z[k]; }
        zp[j] += eps;
        zm[j] -= eps;
        bcbf_closed_loop(zp, fp);
        bcbf_closed_loop(zm, fm);
        for (int i = 0; i < 3; i++) {
            J[i][j] = (fp[i] - fm[i]) / (2.0 * eps);
        }
    }
}

static inline void bcbf_matmul3(double A[3][3], double B[3][3], double C[3][3])
{
    double T[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            T[i][j] = 0.0;
            for (int k = 0; k < 3; k++) T[i][j] += A[i][k] * B[k][j];
        }
    }
    memcpy(C, T, sizeof(T));
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

// Compute all backup rollout/terminal terms H_i(x) and their gradients wrt the
// current original state x0 = [x,y,theta]. This is the same sensitivity-based
// structure as the user's v2 LSE-bCBF code.
static int bcbf_backup_terms_and_grads(double x, double y, double theta,
                                       double H[BCBF_N_BACKUP_TERMS],
                                       double gradH[BCBF_N_BACKUP_TERMS][3])
{
    double z[3] = {x, y, theta};
    double Phi[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);
    int term = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        for (int i = 0; i < BCBF_N_OBS; i++) {
            double dx = z[0] - bcbf_obs[i].cx;
            double dy = z[1] - bcbf_obs[i].cy;
            double grad_z[3] = {2.0 * dx, 2.0 * dy, 0.0};

            H[term] = bcbf_h_obs(i, z[0], z[1]);
            gradH[term][0] = gradH[term][1] = gradH[term][2] = 0.0;
            for (int j = 0; j < 3; j++) {
                for (int m = 0; m < 3; m++) {
                    gradH[term][j] += grad_z[m] * Phi[m][j];
                }
            }
            term++;
        }

        if (k < BCBF_N_SAMPLES - 1) {
            double fz[3], J[3][3], Phidot[3][3];
            bcbf_closed_loop(z, fz);
            bcbf_cl_jacobian_fd(z, J);
            bcbf_matmul3(J, Phi, Phidot);

            for (int i = 0; i < 3; i++) z[i] += dt * fz[i];
            z[2] = bcbf_wrap_pi(z[2]);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) Phi[i][j] += dt * Phidot[i][j];
            }
        }
    }

    for (int i = 0; i < BCBF_N_OBS; i++) {
        double dx = z[0] - bcbf_obs[i].cx;
        double dy = z[1] - bcbf_obs[i].cy;
        double grad_z[3] = {2.0 * dx, 2.0 * dy, 0.0};

        H[term] = bcbf_h_backup_obs(i, z[0], z[1]);
        gradH[term][0] = gradH[term][1] = gradH[term][2] = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int m = 0; m < 3; m++) {
                gradH[term][j] += grad_z[m] * Phi[m][j];
            }
        }
        term++;
    }

    return term;
}

// Composite augmented barrier h_aug(x,u):
//   1) For each backup rollout term H_i(x), construct the relative-degree-2
//      HOCBF pre-control term b_i = Hdot_i + alpha0 H_i, where
//      Hdot_i = gradH_i(x) * f_old(x,u).
//   2) Add input-bound barriers phi(u).
//   3) LSE/soft-min all these relative-degree-1 terms.
static double bcbf_augmented_barrier_value(double x, double y, double theta,
                                           const double u[3])
{
    double H[BCBF_N_BACKUP_TERMS];
    double gradH[BCBF_N_BACKUP_TERMS][3];
    int nH = bcbf_backup_terms_and_grads(x, y, theta, H, gradH);

    double f_old[3];
    bcbf_dynamics(x, y, theta, u, f_old);

    double terms[BCBF_N_AUG_TERMS];
    int n = 0;

    for (int i = 0; i < nH; i++) {
        double Hdot = 0.0;
        for (int j = 0; j < 3; j++) Hdot += gradH[i][j] * f_old[j];
        terms[n++] = Hdot + bcbf_alpha0_hocbf * H[i];
    }

    // Input-bound barriers, relative degree 1 wrt uhat through u_dot.
    terms[n++] = bcbf_u_max[0] - u[0];
    terms[n++] = u[0] - bcbf_u_min[0];
    terms[n++] = bcbf_u_max[1] - u[1];
    terms[n++] = u[1] - bcbf_u_min[1];
    terms[n++] = bcbf_u_max[2] - u[2];
    terms[n++] = u[2] - bcbf_u_min[2];

    return bcbf_softmin_from_terms(terms, n, bcbf_kappa_lse);
}

// Finite-difference gradient of h_aug wrt xhat = [x,y,theta,vx,vy,w].
// This avoids hand-coding Hessian/sensitivity derivatives of the HOCBF/LSE term.
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
static inline void cbf_circle_obstacles_filter_reset(void)
{
    cf_success_count = 0;
    cf_fail_count = 0;
    cbf_log_status = 0;
    cbf_log_mu = 0.0;
    cbf_log_lambda = 0.0;
    // state reset happens lazily on next filter call
    extern int bcbf_cd_initialized;
    bcbf_cd_initialized = 0;
}

// Persistent actual input state u. This is the controller-state x_c for the
// simple choice x_c = u.
int bcbf_cd_initialized = 0;
static double bcbf_cd_u_state[3]  = {0.0, 0.0, 0.0};
static double bcbf_cd_ud_prev[3]  = {0.0, 0.0, 0.0};

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
    double grad[6];
    bcbf_augmented_barrier_grad_fd(x, y, theta, bcbf_cd_u_state, grad);

    // Augmented drift fhat = [f_old(x,u); -A u].
    double f_old[3];
    bcbf_dynamics(x, y, theta, bcbf_cd_u_state, f_old);

    double Lf = 0.0;
    for (int j = 0; j < 3; j++) Lf += grad[j] * f_old[j];
    for (int j = 0; j < 3; j++) Lf += grad[3+j] * (-bcbf_cd_a[j] * bcbf_cd_u_state[j]);

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

    printf("[CD-bCBF-Y] x=%.3f y=%.3f th=%.3f | "
       "h_aug=%.6f omega=%.6f denom=%.6e lambda=%.6f | "
       "ud=[%.3f %.3f %.3f] uhatd=[%.3f %.3f %.3f] "
       "uhat=[%.3f %.3f %.3f] u=[%.3f %.3f %.3f]\n",
       x, y, theta,
       h_aug, omega, denom, lambda,
       ud[0], ud[1], ud[2],
       uhat_d[0], uhat_d[1], uhat_d[2],
       uhat[0], uhat[1], uhat[2],
       bcbf_cd_u_state[0], bcbf_cd_u_state[1], bcbf_cd_u_state[2]);


    double mu = (h_aug * lambda) / bcbf_gamma_mu;

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
    cbf_log_status = 1;
    cbf_log_qp_time_sec = cbf_now_sec() - t0;

    for (int j = 0; j < 3; j++) bcbf_cd_ud_prev[j] = ud[j];

    if (isfinite(h_aug) && isfinite(lambda)) cf_success_count++;
    else cf_fail_count++;

    return 1;
}
