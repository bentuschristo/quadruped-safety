#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <time.h>

// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure
static double cbf_log_qp_time_sec = 0.0;
static double cbf_log_filter_time_sec = 0.0;  // wall time for the whole safety filter call

static inline double cbf_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

// =============================================================================
// Hand-picked blending-function safety filter, v2 zero-action backup for the high-level quadruped model.
// Same public function signature as cbf_circle_obstacles_qp_lim.c.
//
// State: X = [x, y, theta]
// Input: U = [vx, vy, omega] body-frame command
// Dynamics:
//      xdot     = cos(theta) vx - sin(theta) vy
//      ydot     = sin(theta) vx + cos(theta) vy
//      thetadot = omega
//
// Safe-set convention:
//      C_S = { X : h_i(X) >= 0 for every obstacle i }
//      h_i(X) = ||p - p_obs,i||^2 - r_i^2
//
// Backup-set convention, v2:
//      k_b(X) = 0
//      C_B = { X : h_{b,i}(X) >= 0 for every obstacle i }
//      h_{b,i}(X) = ||p - p_obs,i||^2 - (r_i + d_B)^2
//
// Blending law:
//      u = (1 - mu) u_nom + mu u_b
// where mu is hand-picked from the implicit backup-safe margin h_I.
// =============================================================================

typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} BLObstacle;

// #define BL_N_OBS 4
// static BLObstacle bl_obs[BL_N_OBS] = {
//     { -2.0,  1.0, 0.45, 1.0, 1 },
//     { -4.0, -1.0, 0.45, 1.0, 1 },
//     { -2.0, -1.0, 0.45, 1.0, 1 },
//     { -4.0,  1.0, 0.45, 1.0, 1 }
// };

#define BL_N_OBS 5
static BLObstacle bl_obs[BL_N_OBS] = {
    {  0.5, 1.0, 0.40, 0.5, 1 },
    { -0.5, 2.0, 0.40, 0.5, 1 },
    {  0.0, 4.0, 0.40, 0.5, 1 },
    {  0.4, 6.0, 0.40, 0.5, 1 },
    { -1.0, 5.5, 0.40, 0.5, 1 }
};

static double bl_u_max[3] = {  1.0,  0.3,  1.0 };
static double bl_u_min[3] = { -1.0, -0.3, -1.0 };

// ----------------------- Control-dynamics extension --------------------------
// The blending law is applied to the surrogate command uhat. The applied command
// u is an internal state: u_dot = -A u + A uhat. Input limits enter the blending
// margin as current command-state safety terms instead of output saturation.
static double bl_cd_a[3]       = {8.0, 8.0, 8.0};
static double bl_sigma0[3]     = {8.0, 8.0, 8.0};
static double bl_cd_dt         = 0.001;
static int    bl_cd_init       = 0;
static double bl_cd_u[3]       = {0.0, 0.0, 0.0};
static double bl_cd_ud_prev[3] = {0.0, 0.0, 0.0};

static inline void bl_nominal_surrogate(const double ud[3], double uhat_d[3])
{
    for (int j = 0; j < 3; j++) {
        double ud_dot = (ud[j] - bl_cd_ud_prev[j]) / fmax(bl_cd_dt, 1.0e-9);
        uhat_d[j] = bl_cd_u[j]
                    + (ud_dot / bl_cd_a[j])
                    + (bl_sigma0[j] / bl_cd_a[j]) * (ud[j] - bl_cd_u[j]);
    }
}

static inline void bl_update_command_state(const double uhat[3])
{
    for (int j = 0; j < 3; j++) {
        double rho = exp(-bl_cd_a[j] * bl_cd_dt);
        bl_cd_u[j] = rho * bl_cd_u[j] + (1.0 - rho) * uhat[j];
    }
}

static const double bl_backup_margin = 0.15;
static const double bl_kappa_lse = 10.0;  // larger -> closer to hard min

#define BL_N_SAMPLES 20
static const double bl_T_backup = 4.0;
#define BL_MAX_HI_VALUES (BL_N_OBS * BL_N_SAMPLES + BL_N_OBS)


static inline double bl_clip(double v, double lo, double hi)
{
    return fmax(lo, fmin(hi, v));
}

static inline double bl_wrap_pi(double a)
{
    while (a >  M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static inline void bl_dynamics(double x, double y, double th,
                               const double u[3], double f[3])
{
    double c = cos(th);
    double s = sin(th);
    f[0] = c * u[0] - s * u[1];
    f[1] = s * u[0] + c * u[1];
    f[2] = u[2];
}

static inline void bl_backup_controller(double x, double y, double th, double ub[3])
{
    // Zero-action backup: commanded stop.
    (void)x; (void)y; (void)th;
    ub[0] = 0.0;
    ub[1] = 0.0;
    ub[2] = 0.0;
}

static inline double bl_h_obs(int i, double x, double y)
{
    double dx = x - bl_obs[i].cx;
    double dy = y - bl_obs[i].cy;
    return dx * dx + dy * dy - bl_obs[i].r * bl_obs[i].r; // safe iff >= 0
}

static inline double bl_h_backup_obs(int i, double x, double y)
{
    double dx = x - bl_obs[i].cx;
    double dy = y - bl_obs[i].cy;
    double rb = bl_obs[i].r + bl_backup_margin;
    return dx * dx + dy * dy - rb * rb;
}


static double bl_softmin_values(const double *vals, int n)
{
    const double kappa = bl_kappa_lse;
    if (n <= 0) return 1.0e100;

    double zmax = -1.0e300;
    for (int i = 0; i < n; i++) {
        double zi = -kappa * vals[i];
        if (zi > zmax) zmax = zi;
    }

    double sumexp = 0.0;
    for (int i = 0; i < n; i++) sumexp += exp((-kappa * vals[i]) - zmax);

    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        double hmin = vals[0];
        for (int i = 1; i < n; i++) if (vals[i] < hmin) hmin = vals[i];
        return hmin;
    }

    return -(log(sumexp) + zmax) / kappa;
}

static double bl_compute_hI(double x, double y, double theta, const double u_state[3])
{
    double vals[BL_MAX_HI_VALUES + 6];
    int nvals = 0;

    double z[6] = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double dt = bl_T_backup / (double)(BL_N_SAMPLES - 1);

    for (int k = 0; k < BL_N_SAMPLES; k++) {
        for (int i = 0; i < BL_N_OBS; i++) {
            if (!bl_obs[i].active) continue;
            if (nvals < BL_MAX_HI_VALUES + 6) vals[nvals++] = bl_h_obs(i, z[0], z[1]);
        }

        if (k < BL_N_SAMPLES - 1) {
            double c = cos(z[2]);
            double s = sin(z[2]);
            double zdot[6];
            zdot[0] = c * z[3] - s * z[4];
            zdot[1] = s * z[3] + c * z[4];
            zdot[2] = z[5];
            zdot[3] = -bl_cd_a[0] * z[3];
            zdot[4] = -bl_cd_a[1] * z[4];
            zdot[5] = -bl_cd_a[2] * z[5];
            for (int m = 0; m < 6; m++) z[m] += dt * zdot[m];
            z[2] = bl_wrap_pi(z[2]);
        }
    }

    for (int i = 0; i < BL_N_OBS; i++) {
        if (!bl_obs[i].active) continue;
        if (nvals < BL_MAX_HI_VALUES + 6) vals[nvals++] = bl_h_backup_obs(i, z[0], z[1]);
    }

    // Current command-state input-limit margins.
    for (int j = 0; j < 3; j++) {
        if (nvals < BL_MAX_HI_VALUES + 6) vals[nvals++] = u_state[j] - bl_u_min[j];
        if (nvals < BL_MAX_HI_VALUES + 6) vals[nvals++] = bl_u_max[j] - u_state[j];
    }

    return bl_softmin_values(vals, nvals);
}

static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    const double cbf_filter_t0 = cbf_now_sec();
    cbf_log_filter_time_sec = 0.0;
    (void)wx; (void)wy; (void)ww;

    double ud[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    double uhat_nom[3];
    double u_b[3];
    bl_backup_controller(x, y, theta, u_b);

    if (!bl_cd_init) {
        for (int j = 0; j < 3; j++) {
            bl_cd_u[j] = bl_clip(ud[j], bl_u_min[j], bl_u_max[j]);
            bl_cd_ud_prev[j] = ud[j];
        }
        bl_cd_init = 1;
    }

    bl_nominal_surrogate(ud, uhat_nom);
    double hI = bl_compute_hI(x, y, theta, bl_cd_u);

    const double eta = 2.0;
    double hpos = fmax(hI, 0.0);
    double mu = exp(-eta * hpos);
    mu = bl_clip(mu, 0.0, 1.0);

    double uhat[3];
    for (int j = 0; j < 3; j++) {
        uhat[j] = (1.0 - mu) * uhat_nom[j] + mu * u_b[j];
    }
    bl_update_command_state(uhat);

    *vx_cmd = bl_cd_u[0];
    *vy_cmd = bl_cd_u[1];
    *wz_cmd = bl_cd_u[2];

    // Logging
    for (int j = 0; j < 3; j++) {
        cbf_log_u_nom[j]    = ud[j];
        cbf_log_u_backup[j] = u_b[j];
        cbf_log_u_safe[j]   = bl_cd_u[j];
    }
    cbf_log_mu = mu;
    cbf_log_status = 1;

    for (int j = 0; j < 3; j++) bl_cd_ud_prev[j] = ud[j];
    cbf_log_qp_time_sec = 0.0;
    cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
    return (mu > 1e-6);
}
