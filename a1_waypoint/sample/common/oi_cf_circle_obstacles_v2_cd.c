#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <time.h>
#include <string.h>


// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure
static double cbf_log_qp_time_sec = 0.0;
static double cbf_log_filter_time_sec = 0.0;  // wall time for the whole safety filter call


// =============================================================================
// Closed-form Optimal Interpolation safety filter, v2 zero-action backup for the high-level quadruped model.
// Same public function signature as cbf_circle_obstacles_qp_lim.c.
//
// It uses the same backup-horizon constraints as bCBF, but restricts the applied
// command to
//      u(mu) = u_nom + mu (u_b - u_nom),  mu in [0,1].
// Then each constraint becomes scalar-affine in mu:
//      a_i + b_i mu >= 0.
// The code implements the paper-style closed-form 1D OI solution by intersecting interval bounds.
// =============================================================================

typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} OIObstacle;

// #define OI_N_OBS 4
// static OIObstacle oi_obs[OI_N_OBS] = {
//     { -2.0,  1.1, 0.45, 1.0, 1 },
//     { -4.0, -1.1, 0.45, 1.0, 1 },
//     { -2.0, -1.1, 0.45, 1.0, 1 },
//     { -4.0,  1.1, 0.45, 1.0, 1 }
// };

#define OI_N_OBS 5
static OIObstacle oi_obs[OI_N_OBS] = {
    {  0.5, 1.0, 0.40, 0.5, 1 },
    { -0.5, 2.0, 0.40, 0.5, 1 },
    {  0.0, 4.0, 0.40, 0.5, 1 },
    {  0.4, 6.0, 0.40, 0.5, 1 },
    { -1.0, 5.5, 0.40, 0.5, 1 }
};

static double oi_u_max[3] = {  1.0,  0.3,  1.0 };
static double oi_u_min[3] = { -1.0, -0.3, -1.0 };

static const double oi_backup_margin = 0.1;
static const double oi_alpha_B  = 0.8;
static const double oi_kappa_lse_B = 10.0;  // larger -> closer to hard min over backup-set margins

#define OI_N_SAMPLES 20
static const double oi_T_backup = 4.0;
#define OI_N_CONS (OI_N_OBS * OI_N_SAMPLES + OI_N_OBS + 6)

static inline double cbf_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static inline double oi_clip(double v, double lo, double hi)
{
    return fmax(lo, fmin(hi, v));
}

static inline double oi_wrap_pi(double a)
{
    while (a >  M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static inline void oi_dynamics(double x, double y, double th,
                               const double u[3], double f[3])
{
    double c = cos(th);
    double s = sin(th);
    f[0] = c * u[0] - s * u[1];
    f[1] = s * u[0] + c * u[1];
    f[2] = u[2];
}

static inline void oi_input_matrix(double th, double G[3][3])
{
    double c = cos(th);
    double s = sin(th);
    G[0][0] =  c; G[0][1] = -s; G[0][2] = 0.0;
    G[1][0] =  s; G[1][1] =  c; G[1][2] = 0.0;
    G[2][0] = 0.0; G[2][1] = 0.0; G[2][2] = 1.0;
}

// ----------------------- Control-dynamics extension --------------------------
// OI now interpolates the surrogate input uhat, while the actually applied
// command u is an internal state with u_dot = -A u + A uhat.  Input limits are
// enforced as current command-state CBF constraints in the scalar mu problem.
static double oi_cd_a[3]       = {8.0, 8.0, 8.0};
static double oi_sigma0[3]     = {8.0, 8.0, 8.0};
static double oi_alpha_u       = 2.0;
static double oi_cd_dt         = 0.001;
static int    oi_cd_init       = 0;
static double oi_cd_u[3]       = {0.0, 0.0, 0.0};
static double oi_cd_ud_prev[3] = {0.0, 0.0, 0.0};

static inline void oi_fhat_drift_aug(const double z[6], double fhat[6])
{
    double c = cos(z[2]);
    double s = sin(z[2]);
    fhat[0] = c * z[3] - s * z[4];
    fhat[1] = s * z[3] + c * z[4];
    fhat[2] = z[5];
    fhat[3] = -oi_cd_a[0] * z[3];
    fhat[4] = -oi_cd_a[1] * z[4];
    fhat[5] = -oi_cd_a[2] * z[5];
}

static inline void oi_backup_jacobian_aug(const double z[6], double F[6][6])
{
    memset(F, 0, sizeof(double) * 36);
    double th = z[2], vx = z[3], vy = z[4];
    double c = cos(th), s = sin(th);
    F[0][2] = -s * vx - c * vy;
    F[0][3] =  c;
    F[0][4] = -s;
    F[1][2] =  c * vx - s * vy;
    F[1][3] =  s;
    F[1][4] =  c;
    F[2][5] =  1.0;
    F[3][3] = -oi_cd_a[0];
    F[4][4] = -oi_cd_a[1];
    F[5][5] = -oi_cd_a[2];
}

static inline void oi_matmul6(double A[6][6], double B[6][6], double C[6][6])
{
    double T[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            T[i][j] = 0.0;
            for (int k = 0; k < 6; k++) T[i][j] += A[i][k] * B[k][j];
        }
    }
    memcpy(C, T, sizeof(T));
}

static inline void oi_nominal_surrogate(const double ud[3], double uhat_d[3])
{
    for (int j = 0; j < 3; j++) {
        double ud_dot = (ud[j] - oi_cd_ud_prev[j]) / fmax(oi_cd_dt, 1.0e-9);
        uhat_d[j] = oi_cd_u[j]
                    + (ud_dot / oi_cd_a[j])
                    + (oi_sigma0[j] / oi_cd_a[j]) * (ud[j] - oi_cd_u[j]);
    }
}

static inline void oi_update_command_state(const double uhat[3])
{
    for (int j = 0; j < 3; j++) {
        double rho = exp(-oi_cd_a[j] * oi_cd_dt);
        oi_cd_u[j] = rho * oi_cd_u[j] + (1.0 - rho) * uhat[j];
    }
}

static inline void oi_aug_row_to_mu(const double grad0[6], const double z0[6],
                                    double hval, double alpha,
                                    const double uhat_nom[3], const double uhat_b[3],
                                    double *a_mu, double *b_mu)
{
    double fhat0[6];
    oi_fhat_drift_aug(z0, fhat0);
    double Lf = 0.0;
    for (int j = 0; j < 6; j++) Lf += grad0[j] * fhat0[j];
    double Lg[3] = {grad0[3] * oi_cd_a[0], grad0[4] * oi_cd_a[1], grad0[5] * oi_cd_a[2]};
    double duhat[3] = {uhat_b[0] - uhat_nom[0], uhat_b[1] - uhat_nom[1], uhat_b[2] - uhat_nom[2]};
    *a_mu = Lf + alpha * hval
          + Lg[0] * uhat_nom[0] + Lg[1] * uhat_nom[1] + Lg[2] * uhat_nom[2];
    *b_mu = Lg[0] * duhat[0] + Lg[1] * duhat[1] + Lg[2] * duhat[2];
}

static inline void oi_add_input_limit_mu_rows(const double z0[6],
                                              const double uhat_nom[3], const double uhat_b[3],
                                              double a_mu[6], double b_mu[6])
{
    double duhat[3] = {uhat_b[0] - uhat_nom[0], uhat_b[1] - uhat_nom[1], uhat_b[2] - uhat_nom[2]};
    for (int j = 0; j < 3; j++) {
        // lower: -a*u + a*uhat(mu) + alpha*(u-u_min) >= 0
        a_mu[2*j] = -oi_cd_a[j] * z0[3+j] + oi_cd_a[j] * uhat_nom[j]
                    + oi_alpha_u * (z0[3+j] - oi_u_min[j]);
        b_mu[2*j] = oi_cd_a[j] * duhat[j];
        // upper: a*u - a*uhat(mu) + alpha*(u_max-u) >= 0
        a_mu[2*j+1] = oi_cd_a[j] * z0[3+j] - oi_cd_a[j] * uhat_nom[j]
                      + oi_alpha_u * (oi_u_max[j] - z0[3+j]);
        b_mu[2*j+1] = -oi_cd_a[j] * duhat[j];
    }
}


static inline void oi_fhat_with_uhat_aug(const double z[6], const double uhat[3], double zdot[6])
{
    oi_fhat_drift_aug(z, zdot);
    zdot[3] += oi_cd_a[0] * uhat[0];
    zdot[4] += oi_cd_a[1] * uhat[1];
    zdot[5] += oi_cd_a[2] * uhat[2];
}

static inline void oi_matvec6(double A[6][6], const double x[6], double y[6])
{
    double tmp[6];
    for (int i = 0; i < 6; i++) {
        tmp[i] = 0.0;
        for (int j = 0; j < 6; j++) tmp[i] += A[i][j] * x[j];
    }
    for (int i = 0; i < 6; i++) y[i] = tmp[i];
}

static inline void oi_directional_rhs(const double z[6], const double q_p[6], const double q_b[6],
                                      double zdot[6], double qpdot[6], double qbdot[6])
{
    double F[6][6];
    oi_fhat_drift_aug(z, zdot);          // backup surrogate uhat_B = 0
    oi_backup_jacobian_aug(z, F);
    oi_matvec6(F, q_p, qpdot);
    oi_matvec6(F, q_b, qbdot);
}

static inline void oi_rk4_step_backup_state_and_directionals(double z[6], double q_p[6], double q_b[6], double dt)
{
    double k1z[6], k1p[6], k1b[6];
    double k2z[6], k2p[6], k2b[6];
    double k3z[6], k3p[6], k3b[6];
    double k4z[6], k4p[6], k4b[6];
    double ztmp[6], ptmp[6], btmp[6];

    oi_directional_rhs(z, q_p, q_b, k1z, k1p, k1b);

    for (int i = 0; i < 6; i++) {
        ztmp[i] = z[i] + 0.5 * dt * k1z[i];
        ptmp[i] = q_p[i] + 0.5 * dt * k1p[i];
        btmp[i] = q_b[i] + 0.5 * dt * k1b[i];
    }
    ztmp[2] = oi_wrap_pi(ztmp[2]);
    oi_directional_rhs(ztmp, ptmp, btmp, k2z, k2p, k2b);

    for (int i = 0; i < 6; i++) {
        ztmp[i] = z[i] + 0.5 * dt * k2z[i];
        ptmp[i] = q_p[i] + 0.5 * dt * k2p[i];
        btmp[i] = q_b[i] + 0.5 * dt * k2b[i];
    }
    ztmp[2] = oi_wrap_pi(ztmp[2]);
    oi_directional_rhs(ztmp, ptmp, btmp, k3z, k3p, k3b);

    for (int i = 0; i < 6; i++) {
        ztmp[i] = z[i] + dt * k3z[i];
        ptmp[i] = q_p[i] + dt * k3p[i];
        btmp[i] = q_b[i] + dt * k3b[i];
    }
    ztmp[2] = oi_wrap_pi(ztmp[2]);
    oi_directional_rhs(ztmp, ptmp, btmp, k4z, k4p, k4b);

    for (int i = 0; i < 6; i++) {
        z[i]   += (dt / 6.0) * (k1z[i] + 2.0 * k2z[i] + 2.0 * k3z[i] + k4z[i]);
        q_p[i] += (dt / 6.0) * (k1p[i] + 2.0 * k2p[i] + 2.0 * k3p[i] + k4p[i]);
        q_b[i] += (dt / 6.0) * (k1b[i] + 2.0 * k2b[i] + 2.0 * k3b[i] + k4b[i]);
    }
    z[2] = oi_wrap_pi(z[2]);
}

static inline void oi_backup_controller(double x, double y, double th, double ub[3])
{
    // Zero-action backup: commanded stop.
    (void)x; (void)y; (void)th;
    ub[0] = 0.0;
    ub[1] = 0.0;
    ub[2] = 0.0;
}

static inline void oi_closed_loop(double z[3], double zdot[3])
{
    double ub[3];
    oi_backup_controller(z[0], z[1], z[2], ub);
    oi_dynamics(z[0], z[1], z[2], ub, zdot);
}

static void oi_cl_jacobian_fd(const double z[3], double J[3][3])
{
    const double eps = 1e-5;
    double zp[3], zm[3], fp[3], fm[3];

    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) { zp[k] = z[k]; zm[k] = z[k]; }
        zp[j] += eps;
        zm[j] -= eps;
        oi_closed_loop(zp, fp);
        oi_closed_loop(zm, fm);
        for (int i = 0; i < 3; i++) {
            J[i][j] = (fp[i] - fm[i]) / (2.0 * eps);
        }
    }
}

static inline void oi_matmul3(double A[3][3], double B[3][3], double C[3][3])
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

static inline double oi_h_obs(int i, double x, double y)
{
    double dx = x - oi_obs[i].cx;
    double dy = y - oi_obs[i].cy;
    return dx * dx + dy * dy - oi_obs[i].r * oi_obs[i].r; // safe iff >= 0
}

static inline double oi_h_backup_obs(int i, double x, double y)
{
    double dx = x - oi_obs[i].cx;
    double dy = y - oi_obs[i].cy;
    double rb = oi_obs[i].r + oi_backup_margin;
    return dx * dx + dy * dy - rb * rb;
}


static inline void oi_h_backup_lse_and_gradx(double z[3], double Phi[3][3],
                                             double *h_lse, double grad_x[3])
{
    const double kappa = oi_kappa_lse_B;
    double hvals[OI_N_OBS];
    double zmax = -1.0e300;

    for (int i = 0; i < OI_N_OBS; i++) {
        hvals[i] = oi_h_backup_obs(i, z[0], z[1]);
        double zi = -kappa * hvals[i];
        if (zi > zmax) zmax = zi;
    }

    double sumexp = 0.0;
    double grad_z_lse[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < OI_N_OBS; i++) {
        double wi = exp((-kappa * hvals[i]) - zmax);
        sumexp += wi;
        double dx = z[0] - oi_obs[i].cx;
        double dy = z[1] - oi_obs[i].cy;
        grad_z_lse[0] += wi * 2.0 * dx;
        grad_z_lse[1] += wi * 2.0 * dy;
    }

    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int i = 1; i < OI_N_OBS; i++) if (hvals[i] < hvals[imin]) imin = i;
        *h_lse = hvals[imin];
        double dx = z[0] - oi_obs[imin].cx;
        double dy = z[1] - oi_obs[imin].cy;
        grad_z_lse[0] = 2.0 * dx;
        grad_z_lse[1] = 2.0 * dy;
        grad_z_lse[2] = 0.0;
    } else {
        *h_lse = -(log(sumexp) + zmax) / kappa;
        for (int m = 0; m < 3; m++) grad_z_lse[m] /= sumexp;
    }

    grad_x[0] = grad_x[1] = grad_x[2] = 0.0;
    for (int j = 0; j < 3; j++) {
        for (int m = 0; m < 3; m++) grad_x[j] += grad_z_lse[m] * Phi[m][j];
    }
}

static inline void oi_gu_from_grad(double grad_x[3], double theta, double gu[3])
{
    double G[3][3];
    oi_input_matrix(theta, G);
    gu[0] = gu[1] = gu[2] = 0.0;
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) gu[j] += grad_x[k] * G[k][j];
    }
}

static void oi_build_scalar_constraints(double x, double y, double theta,
                                        const double u_state[3],
                                        const double uhat_nom[3], const double uhat_b[3],
                                        double a_mu[OI_N_CONS], double b_mu[OI_N_CONS])
{
    /*
     * Match the ROS2 OI-CF-CD implementation:
     *   - propagate the augmented backup state z_B with uhat_B = 0,
     *   - propagate only the two OI directional sensitivities q_p and q_b,
     *     rather than the full Phi_B matrix,
     *   - build scalar constraints a_i + b_i mu >= 0 for every rollout
     *     obstacle term, every terminal backup-set term, and current-time
     *     input-limit CBF rows.
     */
    double z0[6] = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double z[6]  = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double q_p[6], q_b[6];
    double dt = oi_T_backup / (double)(OI_N_SAMPLES - 1);
    int row = 0;

    // Initialize all rows to inactive-safe rows. This keeps the fixed-size
    // arrays valid even if fewer than OI_N_OBS obstacles are active.
    for (int r = 0; r < OI_N_CONS; r++) {
        a_mu[r] = 1.0;
        b_mu[r] = 0.0;
    }

    // Initial directional vector fields at z0:
    //   q_p(0) = fhat(z0) + ghat(z0) uhat_nom,
    //   q_b(0) = fhat(z0) + ghat(z0) uhat_b.
    oi_fhat_with_uhat_aug(z0, uhat_nom, q_p);
    oi_fhat_with_uhat_aug(z0, uhat_b,   q_b);

    for (int k = 0; k < OI_N_SAMPLES; k++) {
        for (int i = 0; i < OI_N_OBS; i++) {
            if (!oi_obs[i].active) continue;

            double h = oi_h_obs(i, z[0], z[1]);
            double dx = z[0] - oi_obs[i].cx;
            double dy = z[1] - oi_obs[i].cy;
            double grad_local[6] = {2.0 * dx, 2.0 * dy, 0.0, 0.0, 0.0, 0.0};

            double hdot_p = 0.0;
            double hdot_b = 0.0;
            for (int j = 0; j < 6; j++) {
                hdot_p += grad_local[j] * q_p[j];
                hdot_b += grad_local[j] * q_b[j];
            }

            // a + b*mu >= 0, where uhat(mu)=uhat_nom+mu*(uhat_b-uhat_nom)
            a_mu[row] = hdot_p + oi_obs[i].alpha * h;
            b_mu[row] = hdot_b - hdot_p;
            row++;
        }

        if (k < OI_N_SAMPLES - 1) {
            oi_rk4_step_backup_state_and_directionals(z, q_p, q_b, dt);
        }
    }

    // Terminal backup-set rows, one row per active obstacle, matching the
    // real robot OI-CF-CD implementation. No terminal LSE compression here.
    for (int i = 0; i < OI_N_OBS; i++) {
        if (!oi_obs[i].active) continue;

        double hB = oi_h_backup_obs(i, z[0], z[1]);
        double dx = z[0] - oi_obs[i].cx;
        double dy = z[1] - oi_obs[i].cy;
        double grad_local[6] = {2.0 * dx, 2.0 * dy, 0.0, 0.0, 0.0, 0.0};

        double hdot_p = 0.0;
        double hdot_b = 0.0;
        for (int j = 0; j < 6; j++) {
            hdot_p += grad_local[j] * q_p[j];
            hdot_b += grad_local[j] * q_b[j];
        }

        a_mu[row] = hdot_p + oi_alpha_B * hB;
        b_mu[row] = hdot_b - hdot_p;
        row++;
    }

    // Current-time command-state input-limit CBF rows.
    double a_in[6], b_in[6];
    oi_add_input_limit_mu_rows(z0, uhat_nom, uhat_b, a_in, b_in);
    for (int r = 0; r < 6; r++) {
        a_mu[row] = a_in[r];
        b_mu[row] = b_in[r];
        row++;
    }
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
    double uhat_nom[3], uhat_b[3];
    oi_backup_controller(x, y, theta, uhat_b);

    if (!oi_cd_init) {
        for (int j = 0; j < 3; j++) {
            oi_cd_u[j] = oi_clip(ud[j], oi_u_min[j], oi_u_max[j]);
            oi_cd_ud_prev[j] = ud[j];
        }
        oi_cd_init = 1;
    }
    oi_nominal_surrogate(ud, uhat_nom);

    double a_mu[OI_N_CONS];
    double b_mu[OI_N_CONS];
    oi_build_scalar_constraints(x, y, theta, oi_cd_u, uhat_nom, uhat_b, a_mu, b_mu);

    double qp_t0 = cbf_now_sec();
    double mu_low = 0.0;
    double mu_high = 1.0;
    const double eps = 1e-9;
    int infeasible = 0;

    for (int i = 0; i < OI_N_CONS; i++) {
        if (b_mu[i] > eps) {
            double lb = -a_mu[i] / b_mu[i];
            if (lb > mu_low) mu_low = lb;
        } else if (b_mu[i] < -eps) {
            double ub = -a_mu[i] / b_mu[i];
            if (ub < mu_high) mu_high = ub;
        } else {
            if (a_mu[i] < 0.0) infeasible = 1;
        }
    }

    double mu;
    if (infeasible || mu_low > mu_high) {
        mu = 1.0;
    } else {
        mu = oi_clip(mu_low, 0.0, 1.0);
    }
    cbf_log_qp_time_sec = cbf_now_sec() - qp_t0;

    double uhat[3];
    for (int j = 0; j < 3; j++) {
        uhat[j] = uhat_nom[j] + mu * (uhat_b[j] - uhat_nom[j]);
    }
    oi_update_command_state(uhat);

    *vx_cmd = oi_cd_u[0];
    *vy_cmd = oi_cd_u[1];
    *wz_cmd = oi_cd_u[2];

    // Logging
    for (int j = 0; j < 3; j++) {
        cbf_log_u_nom[j]    = ud[j];
        cbf_log_u_backup[j] = uhat_b[j];
    }
    cbf_log_u_safe[0] = *vx_cmd;
    cbf_log_u_safe[1] = *vy_cmd;
    cbf_log_u_safe[2] = *wz_cmd;
    cbf_log_mu = mu;
    cbf_log_status = (infeasible || mu_low > mu_high) ? 0 : 1;

    for (int j = 0; j < 3; j++) oi_cd_ud_prev[j] = ud[j];
    cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
    return (mu > 1e-6);
}
