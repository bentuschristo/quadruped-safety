#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <string.h>
#include <time.h>

// =============================================================================
// Closed-form CBF-QP with bounded inputs using CONTROL DYNAMICS.
// Single-LSE version.
//
// This is intentionally simpler than the previous version:
//
//   1) ONE LSE combines all constraints into one composite barrier.
//   2) No clipping of nominal u_d before computing uhat_d.
//   3) No final output clipping by default.
//   4) Safety is enforced on the augmented state
//
//        xhat = [x, y, theta, vx, vy, omega]^T
//
//      where x_c = u = [vx, vy, omega]^T.
//
// Plant:
//      x_dot     = cos(theta) vx - sin(theta) vy
//      y_dot     = sin(theta) vx + cos(theta) vy
//      theta_dot = omega
//
// Control dynamics:
//      u_dot = -A u + A uhat
//
// Since obstacle h_i(x,y) has relative degree 2 wrt uhat, define
//
//      b_i(xhat) = hdot_i(xhat) + alpha0 h_i(x)
//
// and combine b_i with input-bound barriers using ONE soft-min/LSE:
//
//      hhat(xhat) = softmin( b_1,...,b_N,
//                            u_max-u, u-u_min )
//
// Then apply the single closed-form CBF condition:
//
//      L_fhat hhat + L_ghat hhat uhat + alpha hhat + mu hhat >= 0.
//
// =============================================================================

// ------------------------- Logging variables ---------------------------------
static int cf_success_count = 0;
static int cf_fail_count = 0;

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
static double cbf_log_min_obs_h   = 0.0;
static double cbf_log_min_term    = 0.0;

// --------------------------- Scenario parameters -----------------------------
typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} CircleCBFClosedForm;

#define CBF_N_OBS 4
#define CBF_N_INPUT_TERMS 6
#define CBF_N_TOTAL_TERMS (CBF_N_OBS + CBF_N_INPUT_TERMS)

// Same obstacle positions as cbf_circle_obstacles_qp_lim.c.
static CircleCBFClosedForm obstacles[CBF_N_OBS] = {
    { -2.0,  1.0, 0.35, 0.5, 1 },
    { -4.0, -1.0, 0.35, 0.5, 1 },
    { -2.0, -1.0, 0.35, 0.5, 1 },
    { -4.0,  1.0, 0.35, 0.5, 1 },
};

// Body-frame input limits.
static double u_max[3] = {  1.0,  0.3,  1.0 };
static double u_min[3] = { -1.0, -0.3, -1.0 };

// HOCBF gain for b_i = hdot_i + alpha0 h_i.
static const double cbf_alpha0_hocbf = 0.5;

// Outer CBF gain on composite hhat.
static const double cbf_alpha_outer = 0.5;

// No extra margin by default; match raw obstacle radii.
// Increase this if you want earlier intervention.
static const double cbf_safety_margin = 0.0;

// Single LSE gain over all terms.
// Larger = closer to true min, less conservative, sharper gradients.
static const double cbf_kappa_lse = 80.0;

// Slack penalty for relaxed CBF.
static const double cbf_gamma_mu = 1.0e5;

// Control dynamics u_dot = -A u + A uhat.
// Larger A makes actual command u track uhat faster.
static double cbf_cd_a[3] = { 20.0, 20.0, 20.0 };

// Tracking decay for u -> u_d through uhat_d.
static double cbf_sigma0[3] = { 8.0, 8.0, 8.0 };

// Loop period used for finite-difference ud_dot and integration.
static double cbf_cd_dt = 0.001;

// No output clipping by default. The point here is to test whether the xhat CBF
// enforces input bounds through input-barrier terms.
#define CBF_CD_HARD_CLIP_OUTPUT 0

// Debug print. Set to 0 when running long experiments.
#define CBF_CD_DEBUG_PRINT 1

// ---------------------------- Utility functions ------------------------------
static inline double cbf_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static inline double cbf_clip(double v, double lo, double hi)
{
    return fmax(lo, fmin(hi, v));
}

static inline double cbf_wrap_pi(double a)
{
    while (a >  M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static inline void cbf_body_kinematics(double x, double y, double theta,
                                       const double u[3], double xdot[3])
{
    (void)x; (void)y;
    const double c = cos(theta);
    const double s = sin(theta);
    xdot[0] = c * u[0] - s * u[1];
    xdot[1] = s * u[0] + c * u[1];
    xdot[2] = u[2];
}

static inline void cbf_aug_drift(double x, double y, double theta,
                                 const double u[3], double fhat[6])
{
    double f_old[3];
    cbf_body_kinematics(x, y, theta, u, f_old);

    fhat[0] = f_old[0];
    fhat[1] = f_old[1];
    fhat[2] = f_old[2];

    // Drift part of u_dot = -A u + A uhat.
    fhat[3] = -cbf_cd_a[0] * u[0];
    fhat[4] = -cbf_cd_a[1] * u[1];
    fhat[5] = -cbf_cd_a[2] * u[2];
}

static inline double cbf_h_obs(int i, double x, double y)
{
    double dx = x - obstacles[i].cx;
    double dy = y - obstacles[i].cy;
    return dx * dx + dy * dy - obstacles[i].r * obstacles[i].r - cbf_safety_margin;
}

static inline double cbf_softmin_from_terms(const double *terms, int n, double kappa)
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
        for (int i = 1; i < n; i++) {
            if (terms[i] < terms[imin]) imin = i;
        }
        return terms[imin];
    }

    return -(log(sumexp) + zmax) / kappa;
}

// Build the ONE-LSE terms:
//   obstacle HOCBF terms b_i = hdot_i + alpha0 h_i
//   input barriers.
static int cbf_all_terms(double x, double y, double theta,
                         const double u[3],
                         double terms[CBF_N_TOTAL_TERMS],
                         double *min_raw_obs_h,
                         double *min_term)
{
    double f_old[3];
    cbf_body_kinematics(x, y, theta, u, f_old);

    int n = 0;
    double min_h = 1.0e300;
    double mt = 1.0e300;

    for (int i = 0; i < CBF_N_OBS; i++) {
        if (!obstacles[i].active) continue;

        double dx = x - obstacles[i].cx;
        double dy = y - obstacles[i].cy;
        double h  = cbf_h_obs(i, x, y);
        if (h < min_h) min_h = h;

        // grad h = [2dx, 2dy, 0]
        double hdot = 2.0 * dx * f_old[0] + 2.0 * dy * f_old[1];

        // HOCBF pre-control term.
        terms[n] = hdot + cbf_alpha0_hocbf * h;
        if (terms[n] < mt) mt = terms[n];
        n++;
    }

    // Input-bound barriers on x_c = u.
    terms[n] = u_max[0] - u[0]; if (terms[n] < mt) mt = terms[n]; n++;
    terms[n] = u[0] - u_min[0]; if (terms[n] < mt) mt = terms[n]; n++;

    terms[n] = u_max[1] - u[1]; if (terms[n] < mt) mt = terms[n]; n++;
    terms[n] = u[1] - u_min[1]; if (terms[n] < mt) mt = terms[n]; n++;

    terms[n] = u_max[2] - u[2]; if (terms[n] < mt) mt = terms[n]; n++;
    terms[n] = u[2] - u_min[2]; if (terms[n] < mt) mt = terms[n]; n++;

    if (min_raw_obs_h) *min_raw_obs_h = min_h;
    if (min_term) *min_term = mt;

    return n;
}

static double cbf_augmented_barrier_value(double x, double y, double theta,
                                          const double u[3])
{
    double terms[CBF_N_TOTAL_TERMS];
    double min_h = 0.0;
    double min_t = 0.0;
    int n = cbf_all_terms(x, y, theta, u, terms, &min_h, &min_t);

    cbf_log_min_obs_h = min_h;
    cbf_log_min_term = min_t;

    return cbf_softmin_from_terms(terms, n, cbf_kappa_lse);
}

// Finite-difference gradient of hhat wrt xhat=[x,y,theta,vx,vy,omega].
static void cbf_augmented_barrier_grad_fd(double x, double y, double theta,
                                          const double u[3], double grad[6])
{
    const double eps[6] = {1e-4, 1e-4, 1e-5, 1e-5, 1e-5, 1e-5};

    for (int j = 0; j < 6; j++) {
        double xp = x, yp = y, thp = theta, up[3] = {u[0], u[1], u[2]};
        double xm = x, ym = y, thm = theta, um[3] = {u[0], u[1], u[2]};

        if (j == 0) {
            xp += eps[j]; xm -= eps[j];
        } else if (j == 1) {
            yp += eps[j]; ym -= eps[j];
        } else if (j == 2) {
            thp = cbf_wrap_pi(thp + eps[j]);
            thm = cbf_wrap_pi(thm - eps[j]);
        } else {
            up[j-3] += eps[j];
            um[j-3] -= eps[j];
        }

        double hp = cbf_augmented_barrier_value(xp, yp, thp, up);
        double hm = cbf_augmented_barrier_value(xm, ym, thm, um);
        grad[j] = (hp - hm) / (2.0 * eps[j]);
    }
}

// Persistent actual command state x_c = u.
int cbf_cd_initialized = 0;
static double cbf_cd_u_state[3] = {0.0, 0.0, 0.0};
static double cbf_cd_ud_prev[3] = {0.0, 0.0, 0.0};

static inline void cbf_circle_obstacles_filter_reset(void)
{
    cf_success_count = 0;
    cf_fail_count = 0;
    cbf_log_status = 0;
    cbf_log_mu = 0.0;
    cbf_log_lambda = 0.0;
    cbf_cd_initialized = 0;
}

// =============================================================================
// Public filter function. Same signature as cbf_circle_obstacles_qp_lim.c.
//
// Inputs:
//   x,y,theta: world pose and yaw
//   vx_cmd,vy_cmd,wz_cmd: BODY-FRAME nominal command, modified in-place
//   wx,wy,ww: cost weights
//
// Output:
//   integrated actual command u=[vx,vy,omega], modified in-place.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    double t0 = cbf_now_sec();

    // Nominal body-frame command. Do NOT clip this before uhat_d.
    double ud[3] = {*vx_cmd, *vy_cmd, *wz_cmd};

    double W[3] = {wx, wy, ww};
    for (int j = 0; j < 3; j++) {
        if (W[j] <= 1e-9) W[j] = 1.0;
    }

    cbf_log_u_nom[0] = ud[0];
    cbf_log_u_nom[1] = ud[1];
    cbf_log_u_nom[2] = ud[2];

    cbf_log_u_backup[0] = 0.0;
    cbf_log_u_backup[1] = 0.0;
    cbf_log_u_backup[2] = 0.0;

    if (!cbf_cd_initialized) {
        // Important assumption for forward invariance:
        // x_c(0)=u(0) should already be inside the input safe set.
        // We do not clip here because the user requested no clipping;
        // if the first nominal command is outside the bounds, the augmented
        // state starts unsafe and CBF invariance is not guaranteed.
        for (int j = 0; j < 3; j++) {
            cbf_cd_u_state[j] = ud[j];
            cbf_cd_ud_prev[j] = ud[j];
        }
        cbf_cd_initialized = 1;
    }

    // Desired surrogate input uhat_d that would make actual u track nominal ud.
    // No clipping of ud or uhat_d.
    double ud_dot[3], uhat_d[3];
    for (int j = 0; j < 3; j++) {
        ud_dot[j] = (ud[j] - cbf_cd_ud_prev[j]) / fmax(cbf_cd_dt, 1e-6);
        uhat_d[j] = cbf_cd_u_state[j]
                    + (1.0 / cbf_cd_a[j]) * ud_dot[j]
                    + (cbf_sigma0[j] / cbf_cd_a[j]) * (ud[j] - cbf_cd_u_state[j]);

        cbf_log_uhat_d[j] = uhat_d[j];
    }

    // Composite augmented barrier and gradient.
    double h_aug = cbf_augmented_barrier_value(x, y, theta, cbf_cd_u_state);
    double min_obs_dbg = cbf_log_min_obs_h;
    double min_term_dbg = cbf_log_min_term;

    double grad[6];
    cbf_augmented_barrier_grad_fd(x, y, theta, cbf_cd_u_state, grad);

    // Lf = grad h * fhat, where fhat=[kinematics; -Au].
    double fhat[6];
    cbf_aug_drift(x, y, theta, cbf_cd_u_state, fhat);

    double Lf = 0.0;
    for (int j = 0; j < 6; j++) {
        Lf += grad[j] * fhat[j];
    }

    // Lg = grad_u * A because ghat=[0; A].
    double Lg[3];
    for (int j = 0; j < 3; j++) {
        Lg[j] = grad[3+j] * cbf_cd_a[j];
    }

    // Closed-form KKT solution:
    // min 0.5(uhat-uhat_d)'W(uhat-uhat_d) + 0.5 gamma mu^2
    // s.t. Lf + Lg uhat + alpha h + mu h >= 0.
    double omega = Lf + cbf_alpha_outer * h_aug;
    for (int j = 0; j < 3; j++) {
        omega += Lg[j] * uhat_d[j];
    }

    double denom = (h_aug * h_aug) / cbf_gamma_mu;
    for (int j = 0; j < 3; j++) {
        denom += (Lg[j] * Lg[j]) / W[j];
    }

    double lambda = 0.0;
    if (omega < 0.0 && denom > 1e-12 && isfinite(denom)) {
        lambda = -omega / denom;
    }

    double uhat[3];
    for (int j = 0; j < 3; j++) {
        uhat[j] = uhat_d[j] + (Lg[j] / W[j]) * lambda;
        cbf_log_uhat[j] = uhat[j];
    }

    double mu = (h_aug * lambda) / cbf_gamma_mu;

#if CBF_CD_DEBUG_PRINT
    printf("[CF-CBF-CD-LSE1] x=%.3f y=%.3f th=%.3f | "
           "h=%.6f minTerm=%.6f minObs=%.6f "
           "omega=%.6f den=%.3e lam=%.6f | "
           "ud=[%.3f %.3f %.3f] uhatd=[%.3f %.3f %.3f] "
           "uhat=[%.3f %.3f %.3f] u=[%.3f %.3f %.3f]\n",
           x, y, theta,
           h_aug, min_term_dbg, min_obs_dbg,
           omega, denom, lambda,
           ud[0], ud[1], ud[2],
           uhat_d[0], uhat_d[1], uhat_d[2],
           uhat[0], uhat[1], uhat[2],
           cbf_cd_u_state[0], cbf_cd_u_state[1], cbf_cd_u_state[2]);
#endif

    // Integrate actual command state u_dot=-A u + A uhat.
    for (int j = 0; j < 3; j++) {
        cbf_cd_u_state[j] += cbf_cd_dt * (-cbf_cd_a[j] * cbf_cd_u_state[j]
                                          + cbf_cd_a[j] * uhat[j]);

#if CBF_CD_HARD_CLIP_OUTPUT
        cbf_cd_u_state[j] = cbf_clip(cbf_cd_u_state[j], u_min[j], u_max[j]);
#endif
    }

    *vx_cmd = cbf_cd_u_state[0];
    *vy_cmd = cbf_cd_u_state[1];
    *wz_cmd = cbf_cd_u_state[2];

    cbf_log_u_safe[0] = cbf_cd_u_state[0];
    cbf_log_u_safe[1] = cbf_cd_u_state[1];
    cbf_log_u_safe[2] = cbf_cd_u_state[2];

    cbf_log_mu = mu;
    cbf_log_lambda = lambda;
    cbf_log_h_aug = h_aug;
    cbf_log_omega = omega;
    cbf_log_status = 1;
    cbf_log_qp_time_sec = cbf_now_sec() - t0;

    for (int j = 0; j < 3; j++) {
        cbf_cd_ud_prev[j] = ud[j];
    }

    if (isfinite(h_aug) && isfinite(lambda)) cf_success_count++;
    else cf_fail_count++;

    return 1;
}
