#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <time.h>
#include <string.h>
#include "osqp.h"


// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure
static double cbf_log_qp_time_sec = 0.0;
static double cbf_log_filter_time_sec = 0.0;  // wall time for the whole safety filter call


// =============================================================================
// OSQP-based Optimal Interpolation safety filter, v2 zero-action backup for the high-level quadruped model.
// Same public function signature as cbf_circle_obstacles_qp_lim.c.
//
// It uses the same backup-horizon constraints as bCBF, but restricts the applied
// command to
//      u(mu) = u_nom + mu (u_b - u_nom),  mu in [0,1].
// Then each constraint becomes scalar-affine in mu:
//      a_i + b_i mu >= 0.
// This file solves the scalar 1D OI-QP with OSQP. Use the oi_cf_* files for the paper-style closed-form solution.
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
#define OI_N_CONS (OI_N_OBS * OI_N_SAMPLES + 1 + 6)

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
    double z0[6] = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double z[6]  = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double Phi[6][6] = {{0}};
    for (int i = 0; i < 6; i++) Phi[i][i] = 1.0;
    double dt = oi_T_backup / (double)(OI_N_SAMPLES - 1);
    int row = 0;

    for (int k = 0; k < OI_N_SAMPLES; k++) {
        for (int i = 0; i < OI_N_OBS; i++) {
            double dx = z[0] - oi_obs[i].cx;
            double dy = z[1] - oi_obs[i].cy;
            double grad_z[6] = {2.0 * dx, 2.0 * dy, 0.0, 0.0, 0.0, 0.0};
            double grad0[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
            for (int j = 0; j < 6; j++) for (int m = 0; m < 6; m++) grad0[j] += grad_z[m] * Phi[m][j];
            double h = oi_h_obs(i, z[0], z[1]);
            oi_aug_row_to_mu(grad0, z0, h, oi_obs[i].alpha, uhat_nom, uhat_b,
                             &a_mu[row], &b_mu[row]);
            row++;
        }
        if (k < OI_N_SAMPLES - 1) {
            double fz[6], F[6][6], Phidot[6][6];
            oi_fhat_drift_aug(z, fz);
            oi_backup_jacobian_aug(z, F);
            oi_matmul6(F, Phi, Phidot);
            for (int m = 0; m < 6; m++) z[m] += dt * fz[m];
            z[2] = oi_wrap_pi(z[2]);
            for (int ii = 0; ii < 6; ii++) for (int jj = 0; jj < 6; jj++) Phi[ii][jj] += dt * Phidot[ii][jj];
        }
    }

    // Terminal backup set as one LSE row.
    const double kappa = oi_kappa_lse_B;
    double hvals[OI_N_OBS];
    double zmax = -1.0e300;
    for (int i = 0; i < OI_N_OBS; i++) {
        hvals[i] = oi_h_backup_obs(i, z[0], z[1]);
        double zi = -kappa * hvals[i];
        if (zi > zmax) zmax = zi;
    }
    double sumexp = 0.0;
    double grad_z_lse[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for (int i = 0; i < OI_N_OBS; i++) {
        double wi = exp((-kappa * hvals[i]) - zmax);
        sumexp += wi;
        grad_z_lse[0] += wi * 2.0 * (z[0] - oi_obs[i].cx);
        grad_z_lse[1] += wi * 2.0 * (z[1] - oi_obs[i].cy);
    }
    double hB_lse;
    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int i = 1; i < OI_N_OBS; i++) if (hvals[i] < hvals[imin]) imin = i;
        hB_lse = hvals[imin];
        memset(grad_z_lse, 0, sizeof(grad_z_lse));
        grad_z_lse[0] = 2.0 * (z[0] - oi_obs[imin].cx);
        grad_z_lse[1] = 2.0 * (z[1] - oi_obs[imin].cy);
    } else {
        hB_lse = -(log(sumexp) + zmax) / kappa;
        for (int j = 0; j < 6; j++) grad_z_lse[j] /= sumexp;
    }
    double grad0B[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for (int j = 0; j < 6; j++) for (int m = 0; m < 6; m++) grad0B[j] += grad_z_lse[m] * Phi[m][j];
    oi_aug_row_to_mu(grad0B, z0, hB_lse, oi_alpha_B, uhat_nom, uhat_b,
                     &a_mu[row], &b_mu[row]);
    row++;

    double a_in[6], b_in[6];
    oi_add_input_limit_mu_rows(z0, uhat_nom, uhat_b, a_in, b_in);
    for (int r = 0; r < 6; r++) {
        a_mu[row] = a_in[r];
        b_mu[row] = b_in[r];
        row++;
    }
}

// -------------------------- OSQP scalar OI-QP state --------------------------// -------------------------- OSQP scalar OI-QP state --------------------------
// QP solved online:
//      min_mu  0.5 * mu^2
//      s.t.    a_i + b_i * mu >= 0,  i = 1,...,OI_N_CONS
//              0 <= mu <= 1
// The objective is equivalent to minimum intervention along the interpolation
// direction because ||u(mu)-u_nom||^2 = mu^2 ||u_b-u_nom||^2. The coefficient
// does not change the optimizer as long as it is positive.
#define OI_QP_N_ROWS (OI_N_CONS + 1)

static OSQPSolver   *oi_qp_solver = NULL;
static OSQPSettings  oi_qp_settings;

static OSQPFloat oi_qp_P_x[1] = {1.0};
static OSQPInt   oi_qp_P_i[1] = {0};
static OSQPInt   oi_qp_P_p[2] = {0, 1};
static OSQPCscMatrix oi_qp_P_mat;

static OSQPFloat oi_qp_A_x[OI_QP_N_ROWS];
static OSQPInt   oi_qp_A_i[OI_QP_N_ROWS];
static OSQPInt   oi_qp_A_p[2] = {0, OI_QP_N_ROWS};
static OSQPCscMatrix oi_qp_A_mat;

static OSQPFloat oi_qp_q[1] = {0.0};
static OSQPFloat oi_qp_l[OI_QP_N_ROWS];
static OSQPFloat oi_qp_u[OI_QP_N_ROWS];

static int oi_qp_init_solver(void)
{
    oi_qp_P_mat.m = 1;
    oi_qp_P_mat.n = 1;
    oi_qp_P_mat.nz = -1;
    oi_qp_P_mat.nzmax = 1;
    oi_qp_P_mat.x = oi_qp_P_x;
    oi_qp_P_mat.i = oi_qp_P_i;
    oi_qp_P_mat.p = oi_qp_P_p;

    for (int r = 0; r < OI_QP_N_ROWS; r++) {
        oi_qp_A_i[r] = r;
        oi_qp_A_x[r] = 0.0;
        oi_qp_l[r] = -OSQP_INFTY;
        oi_qp_u[r] = OSQP_INFTY;
    }

    // Last row is the simple bound row: 0 <= 1*mu <= 1.
    oi_qp_A_x[OI_N_CONS] = 1.0;
    oi_qp_l[OI_N_CONS] = 0.0;
    oi_qp_u[OI_N_CONS] = 1.0;

    oi_qp_A_mat.m = OI_QP_N_ROWS;
    oi_qp_A_mat.n = 1;
    oi_qp_A_mat.nz = -1;
    oi_qp_A_mat.nzmax = OI_QP_N_ROWS;
    oi_qp_A_mat.x = oi_qp_A_x;
    oi_qp_A_mat.i = oi_qp_A_i;
    oi_qp_A_mat.p = oi_qp_A_p;

    osqp_set_default_settings(&oi_qp_settings);
    oi_qp_settings.verbose = 0;
    oi_qp_settings.warm_starting = 1;
    oi_qp_settings.eps_abs = 1e-4;
    oi_qp_settings.eps_rel = 1e-4;
    oi_qp_settings.max_iter = 10000;

    OSQPInt flag = osqp_setup(&oi_qp_solver, &oi_qp_P_mat, oi_qp_q,
                              &oi_qp_A_mat, oi_qp_l, oi_qp_u,
                              OI_QP_N_ROWS, 1, &oi_qp_settings);
    if (flag != 0) {
        fprintf(stderr, "[OI-QP] OSQP setup failed: %d\n", (int)flag);
        return 0;
    }
    return 1;
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

    if (oi_qp_solver == NULL) {
        if (!oi_qp_init_solver()) {
            cbf_log_qp_time_sec = 0.0;
            double mu = 1.0;
            double uhat_fb[3];
            for (int j = 0; j < 3; j++) uhat_fb[j] = uhat_b[j];
            oi_update_command_state(uhat_fb);
            *vx_cmd = oi_cd_u[0];
            *vy_cmd = oi_cd_u[1];
            *wz_cmd = oi_cd_u[2];
            for (int j = 0; j < 3; j++) { cbf_log_u_nom[j] = ud[j]; cbf_log_u_backup[j] = uhat_b[j]; }
            cbf_log_u_safe[0] = *vx_cmd; cbf_log_u_safe[1] = *vy_cmd; cbf_log_u_safe[2] = *wz_cmd;
            cbf_log_mu = mu; cbf_log_status = 0;
            for (int j = 0; j < 3; j++) oi_cd_ud_prev[j] = ud[j];
            cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
            return 1;
        }
    }

    for (int i = 0; i < OI_N_CONS; i++) {
        oi_qp_A_x[i] = (OSQPFloat)b_mu[i];
        oi_qp_l[i] = (OSQPFloat)(-a_mu[i]);
        oi_qp_u[i] = (OSQPFloat)OSQP_INFTY;
    }
    oi_qp_A_x[OI_N_CONS] = 1.0;
    oi_qp_l[OI_N_CONS] = 0.0;
    oi_qp_u[OI_N_CONS] = 1.0;
    oi_qp_q[0] = 0.0;

    osqp_update_data_vec(oi_qp_solver, oi_qp_q, oi_qp_l, oi_qp_u);
    osqp_update_data_mat(oi_qp_solver, NULL, NULL, 0, oi_qp_A_x, NULL, oi_qp_A_mat.nzmax);

    double qp_t0 = cbf_now_sec();
    osqp_solve(oi_qp_solver);
    cbf_log_qp_time_sec = cbf_now_sec() - qp_t0;

    int status_val = oi_qp_solver->info->status_val;
    int status_ok = (status_val == OSQP_SOLVED || status_val == OSQP_SOLVED_INACCURATE);
    double mu;
    if (status_ok && oi_qp_solver->solution != NULL) {
        mu = oi_clip((double)oi_qp_solver->solution->x[0], 0.0, 1.0);
    } else {
        fprintf(stderr, "[OI-QP-CD] OSQP failed (%s). Falling back to backup surrogate.\n", oi_qp_solver->info->status);
        mu = 1.0;
    }
    int infeasible = !status_ok;
    double mu_low = 0.0, mu_high = 1.0;

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
    cbf_log_status = status_ok ? 1 : 0;

    for (int j = 0; j < 3; j++) oi_cd_ud_prev[j] = ud[j];
    cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
    return (mu > 1e-6);
}
