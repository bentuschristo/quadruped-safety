#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>


// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure


// =============================================================================
// Optimal Interpolation safety filter for the high-level quadruped model.
// Same public function signature as cbf_circle_obstacles_qp_lim.c.
//
// It uses the same backup-horizon constraints as bCBF, but restricts the applied
// command to
//      u(mu) = u_nom + mu (u_b - u_nom),  mu in [0,1].
// Then each constraint becomes scalar-affine in mu:
//      a_i + b_i mu >= 0.
// The code solves the 1D OI-QP in closed form by intersecting interval bounds.
// =============================================================================

typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} OIObstacle;

#define OI_N_OBS 4
static OIObstacle oi_obs[OI_N_OBS] = {
    { -2.0,  1.0, 0.45, 1.0, 1 },
    { -4.0, -1.0, 0.45, 1.0, 1 },
    { -2.0, -1.0, 0.45, 1.0, 1 },
    { -4.0,  1.0, 0.45, 1.0, 1 }
};

// static OIObstacle oi_obs[OI_N_OBS] = {
//     { -2.0,  1.0, 0.42, 1.0, 1 },
//     { -4.0, -1.0, 0.42, 1.0, 1 },
//     { -2.0, -1.0, 0.42, 1.0, 1 },
//     { -4.0,  1.0, 0.42, 1.0, 1 },
//     { -2.8,  0.0, 0.30, 1.0, 1 },
//     { -3.3,  0.0, 0.30, 1.0, 1 }
// };

static double oi_u_max[3] = {  1.0,  0.3,  1.0 };
static double oi_u_min[3] = { -1.0, -0.3, -1.0 };

static const double oi_backup_y = 1.5;
static const double oi_alpha_B  = 0.8;
static const double oi_ky_backup = 0.8;
static const double oi_vup_max   = 0.55;
static const double oi_kw_backup = 0.0;

#define OI_N_SAMPLES 20
static const double oi_T_backup = 4.0;
#define OI_N_CONS (OI_N_OBS * OI_N_SAMPLES + 1)

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

static inline void oi_backup_controller(double x, double y, double th, double ub[3])
{
    (void)x;

    double v_up = oi_ky_backup * (oi_backup_y - y);
    if (v_up < 0.0) v_up = 0.0;
    v_up = oi_clip(v_up, 0.0, oi_vup_max);

    // Convert desired world velocity [0, v_up] to body-frame command.
    ub[0] = oi_clip(v_up * sin(th), oi_u_min[0], oi_u_max[0]);
    ub[1] = oi_clip(v_up * cos(th), oi_u_min[1], oi_u_max[1]);
    ub[2] = oi_clip(oi_kw_backup * 0.0, oi_u_min[2], oi_u_max[2]);
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

static inline double oi_h_backup(double x, double y)
{
    (void)x;
    return y - oi_backup_y; // backup set y >= backup_y iff h_b >= 0
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
                                        const double u_nom[3], const double u_b[3],
                                        double a_mu[OI_N_CONS], double b_mu[OI_N_CONS])
{
    double z[3] = {x, y, theta};
    double Phi[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double dt = oi_T_backup / (double)(OI_N_SAMPLES - 1);
    int row = 0;

    double du[3] = {
        u_b[0] - u_nom[0],
        u_b[1] - u_nom[1],
        u_b[2] - u_nom[2]
    };

    for (int k = 0; k < OI_N_SAMPLES; k++) {
        for (int i = 0; i < OI_N_OBS; i++) {
            double dx = z[0] - oi_obs[i].cx;
            double dy = z[1] - oi_obs[i].cy;
            double grad_z[3] = {2.0 * dx, 2.0 * dy, 0.0};
            double grad_x[3] = {0.0, 0.0, 0.0};

            // grad_x = grad_z * Phi
            for (int j = 0; j < 3; j++) {
                for (int m = 0; m < 3; m++) grad_x[j] += grad_z[m] * Phi[m][j];
            }

            double gu[3];
            oi_gu_from_grad(grad_x, theta, gu);

            // Constraint: gu*u(mu) + alpha*h >= 0.
            // With u(mu) = u_nom + mu*du:
            //      a_mu + b_mu * mu >= 0.
            double h = oi_h_obs(i, z[0], z[1]);
            a_mu[row] = gu[0] * u_nom[0] + gu[1] * u_nom[1] + gu[2] * u_nom[2]
                      + oi_obs[i].alpha * h;
            b_mu[row] = gu[0] * du[0]    + gu[1] * du[1]    + gu[2] * du[2];
            row++;
        }

        if (k < OI_N_SAMPLES - 1) {
            double fz[3];
            double J[3][3];
            double Phidot[3][3];
            oi_closed_loop(z, fz);
            oi_cl_jacobian_fd(z, J);
            oi_matmul3(J, Phi, Phidot);

            for (int m = 0; m < 3; m++) z[m] += dt * fz[m];
            z[2] = oi_wrap_pi(z[2]);
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) Phi[ii][jj] += dt * Phidot[ii][jj];
            }
        }
    }

    // Terminal backup-set constraint h_b(phi_b(T,x)) >= 0.
    // h_b = y - backup_y, grad_z h_b = [0,1,0].
    double grad_zB[3] = {0.0, 1.0, 0.0};
    double grad_xB[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 3; j++) {
        for (int m = 0; m < 3; m++) grad_xB[j] += grad_zB[m] * Phi[m][j];
    }

    double guB[3];
    oi_gu_from_grad(grad_xB, theta, guB);
    double hb = oi_h_backup(z[0], z[1]);

    a_mu[row] = guB[0] * u_nom[0] + guB[1] * u_nom[1] + guB[2] * u_nom[2]
              + oi_alpha_B * hb;
    b_mu[row] = guB[0] * du[0]    + guB[1] * du[1]    + guB[2] * du[2];
}

static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    (void)wx; (void)wy; (void)ww;

    double u_nom[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    for (int j = 0; j < 3; j++) u_nom[j] = oi_clip(u_nom[j], oi_u_min[j], oi_u_max[j]);

    double u_b[3];
    oi_backup_controller(x, y, theta, u_b);

    double a_mu[OI_N_CONS];
    double b_mu[OI_N_CONS];
    oi_build_scalar_constraints(x, y, theta, u_nom, u_b, a_mu, b_mu);

    // Solve min 0.5*mu^2 subject to a_i + b_i*mu >= 0 and 0 <= mu <= 1.
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

    int oi_feasible = !(infeasible || mu_low > mu_high);

    double mu;
    if (!oi_feasible) {
        // No interpolation parameter satisfies all constraints. Fall back to
        // full backup command. This is the structural failure mode we want to
        // expose when OI is too restrictive.
        mu = 1.0;
    } else {
        // Objective is minimized by the smallest feasible nonnegative mu.
        mu = oi_clip(mu_low, 0.0, 1.0);
    }

    double u[3];
    for (int j = 0; j < 3; j++) {
        u[j] = u_nom[j] + mu * (u_b[j] - u_nom[j]);
        u[j] = oi_clip(u[j], oi_u_min[j], oi_u_max[j]);
    }

    *vx_cmd = u[0];
    *vy_cmd = u[1];
    *wz_cmd = u[2];

    // Logging
    for (int j = 0; j < 3; j++) {
        cbf_log_u_nom[j]    = u_nom[j];
        cbf_log_u_backup[j] = u_b[j];
    }

    cbf_log_u_safe[0] = *vx_cmd;
    cbf_log_u_safe[1] = *vy_cmd;
    cbf_log_u_safe[2] = *wz_cmd;

    cbf_log_mu = mu;
    cbf_log_status = oi_feasible ? 1 : 0;

    return (mu > 1e-6);
}
