#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>

// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure

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

#define BL_N_OBS 6
static BLObstacle bl_obs[BL_N_OBS] = {
    { -2.0,  1.0, 0.42, 1.0, 1 },
    { -4.0, -1.0, 0.42, 1.0, 1 },
    { -2.0, -1.0, 0.42, 1.0, 1 },
    { -4.0,  1.0, 0.42, 1.0, 1 },
    { -2.8,  0.0, 0.30, 1.0, 1 },
    { -3.3,  0.0, 0.30, 1.0, 1 }
};

static double bl_u_max[3] = {  1.0,  0.3,  1.0 };
static double bl_u_min[3] = { -1.0, -0.3, -1.0 };

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

static double bl_compute_hI(double x, double y, double th)
{
    double z[3] = {x, y, th};
    double dt = bl_T_backup / (double)(BL_N_SAMPLES - 1);
    double vals[BL_MAX_HI_VALUES];
    int nvals = 0;

    for (int k = 0; k < BL_N_SAMPLES; k++) {
        // Obstacle margins along the backup rollout.
        for (int i = 0; i < BL_N_OBS; i++) {
            if (!bl_obs[i].active) continue;
            if (nvals < BL_MAX_HI_VALUES) vals[nvals++] = bl_h_obs(i, z[0], z[1]);
        }

        if (k < BL_N_SAMPLES - 1) {
            double ub[3], fz[3];
            bl_backup_controller(z[0], z[1], z[2], ub);
            bl_dynamics(z[0], z[1], z[2], ub, fz);
            for (int m = 0; m < 3; m++) z[m] += dt * fz[m];
            z[2] = bl_wrap_pi(z[2]);
        }
    }

    // Terminal backup-set margins: enlarged-obstacle clearances.
    // The multiple backup-set conditions are combined using LSE/soft-min.
    for (int i = 0; i < BL_N_OBS; i++) {
        if (!bl_obs[i].active) continue;
        if (nvals < BL_MAX_HI_VALUES) vals[nvals++] = bl_h_backup_obs(i, z[0], z[1]);
    }

    return bl_softmin_values(vals, nvals);
}

static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    (void)wx; (void)wy; (void)ww;

    double u_nom[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    for (int j = 0; j < 3; j++) u_nom[j] = bl_clip(u_nom[j], bl_u_min[j], bl_u_max[j]);

    double u_b[3];
    bl_backup_controller(x, y, theta, u_b);

    double hI = bl_compute_hI(x, y, theta);

    // Paper-style exponential blending:
    // mu = exp(-eta * max(hI, 0))
    // If hI <= 0, mu = 1 and the backup controller is used.
    // If hI > 0, mu decays toward 0 as the system moves deeper inside C_I.
    const double eta = 2.0;   // tune this
    double hpos = fmax(hI, 0.0);
    double mu = exp(-eta * hpos);
    mu = bl_clip(mu, 0.0, 1.0);

    double u[3];
    for (int j = 0; j < 3; j++) {
        u[j] = (1.0 - mu) * u_nom[j] + mu * u_b[j];
        u[j] = bl_clip(u[j], bl_u_min[j], bl_u_max[j]);
    }

    *vx_cmd = u[0];
    *vy_cmd = u[1];
    *wz_cmd = u[2];

    // Logging
    for (int j = 0; j < 3; j++) {
        cbf_log_u_nom[j]    = u_nom[j];
        cbf_log_u_backup[j] = u_b[j];
        cbf_log_u_safe[j]   = u[j];
    }
    cbf_log_mu = mu;
    cbf_log_status = 1;

    return (mu > 1e-6);
}
