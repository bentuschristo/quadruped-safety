#pragma once
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <string.h>
#include "osqp.h"
#include <time.h>

static int qp_success_count = 0;
static int qp_fail_count = 0;

// Logging the control inputs
static double cbf_log_u_nom[3]    = {0.0, 0.0, 0.0};
static double cbf_log_u_backup[3] = {0.0, 0.0, 0.0};
static double cbf_log_u_safe[3]   = {0.0, 0.0, 0.0};

static double cbf_log_mu     = -1.0;  // -1 for bCBF, actual mu for blending/OI
static int    cbf_log_status = 0;     // 1 = success/active, 0 = fallback/inactive/failure
static double cbf_log_qp_time_sec = 0.0;
static double cbf_log_filter_time_sec = 0.0;  // wall time for the whole safety filter call
// =============================================================================
// LSE-aggregated Backup-CBF filter, v2 zero-action backup, for the high-level quadruped model
// State:   X = [x, y, theta]
// Input:   U = [vx, vy, omega]  body-frame command
// Dynamics:
//      xdot     = cos(theta) vx - sin(theta) vy
//      ydot     = sin(theta) vx + cos(theta) vy
//      thetadot = omega
//
// Safe-set convention used here:
//      C_S = { X : h_i(X) >= 0 for every obstacle i }
//      h_i(X) = ||p - p_obs,i||^2 - r_i^2
//
// Backup-set convention, v2:
//      k_b(X) = 0  (zero high-level velocity command)
//      C_B = { X : h_{b,i}(X) >= 0 for every obstacle i }
//      h_{b,i}(X) = ||p - p_obs,i||^2 - (r_i + d_B)^2
//
// In words: the backup action is to stop. The backup set is the set of
// positions with an extra clearance margin from every circular obstacle.
// Since the high-level model is kinematic, zero command gives phi_b(tau,x)=x.
//
// The bCBF-QP uses LSE/soft-min to aggregate all horizon obstacle constraints
// and the terminal enlarged-clearance backup constraints into ONE smooth QP constraint.
// =============================================================================

// --------------------------- Scenario parameters -----------------------------
typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int active;
} BCBFObstacle;

// Original four obstacles plus two center obstacles to create a lateral-escape
// scenario where scalar interpolation methods can fail.
#define BCBF_N_OBS 4
static BCBFObstacle bcbf_obs[BCBF_N_OBS] = {
    { -2.0,  1.1, 0.45, 0.5, 1 }, 
    { -4.0, -1.1, 0.45, 0.5, 1 },
    { -2.0, -1.1, 0.45, 0.5, 1 },
    { -4.0,  1.1, 0.45, 0.5, 1 }
};

// Input limits. Keep these consistent with the rest of the high-level controller.
static double bcbf_u_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_min[3] = { -1.0, -0.3, -1.0 };

// Backup set, v2: enlarged-obstacle clearance set.
// C_B = { h_{b,i}(x,y) >= 0 for all i }, where
// h_{b,i} = ||p - p_i||^2 - (r_i + bcbf_backup_margin)^2.
static const double bcbf_backup_margin = 0.1;
static const double bcbf_alpha_B = 0.8;
static const double bcbf_kappa_lse_B = 30.0;  // larger -> closer to hard min over backup-set margins
static const double bcbf_kappa_lse = 30.0;    // larger -> closer to hard min over all horizon/terminal terms
static const double bcbf_alpha_lse = 0.8;     // class-K gain for the single LSE aggregate barrier

// Backup rollout discretization. Longer horizon helps the backup controller
// certify recovery to the right-side backup set from states left of the origin.
#define BCBF_N_SAMPLES 20
static const double bcbf_T_backup = 4.0;

// Total constraints = obstacle constraints at each backup sample + terminal hb + boxes
// LSE terms = obstacle margins at all backup samples + terminal enlarged-obstacle backup margins.
// LSE aggregation reduces all bCBF horizon/terminal constraints to ONE QP row.
#define BCBF_N_BACKUP_LSE_TERMS (BCBF_N_OBS * BCBF_N_SAMPLES + BCBF_N_OBS)
#define BCBF_N_INPUT_LSE_TERMS 6
#define BCBF_N_LSE_TERMS (BCBF_N_BACKUP_LSE_TERMS + BCBF_N_INPUT_LSE_TERMS)
#define BCBF_N_CBF_ROWS  1
#define BCBF_N_ROWS      1


static inline double cbf_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

// ---------------------------- Utility functions ------------------------------
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

static inline void bcbf_input_matrix(double th, double G[3][3])
{
    double c = cos(th);
    double s = sin(th);
    G[0][0] =  c; G[0][1] = -s; G[0][2] = 0.0;
    G[1][0] =  s; G[1][1] =  c; G[1][2] = 0.0;
    G[2][0] = 0.0; G[2][1] = 0.0; G[2][2] = 1.0;
}

// ----------------------- Control-dynamics extension --------------------------
// The QP decision variable is now the surrogate command uhat.  The actually
// applied command u is kept as an internal state and evolves as
//      u_dot = -A u + A uhat.
// Input limits are therefore handled as CBF state constraints on u, not as
// static box constraints on the QP decision variable.
static double bcbf_cd_a[3]     = {8.0, 8.0, 8.0};
static double bcbf_sigma0[3]   = {8.0, 8.0, 8.0};
static double bcbf_alpha_u     = 2.0;
static double bcbf_cd_dt       = 0.001;
static int    bcbf_cd_init     = 0;
static double bcbf_cd_u[3]     = {0.0, 0.0, 0.0};
static double bcbf_cd_ud_prev[3] = {0.0, 0.0, 0.0};

static inline void bcbf_fhat_drift_aug(const double z[6], double fhat[6])
{
    double c = cos(z[2]);
    double s = sin(z[2]);
    fhat[0] = c * z[3] - s * z[4];
    fhat[1] = s * z[3] + c * z[4];
    fhat[2] = z[5];
    fhat[3] = -bcbf_cd_a[0] * z[3];
    fhat[4] = -bcbf_cd_a[1] * z[4];
    fhat[5] = -bcbf_cd_a[2] * z[5];
}

static inline void bcbf_backup_jacobian_aug(const double z[6], double F[6][6])
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
    F[3][3] = -bcbf_cd_a[0];
    F[4][4] = -bcbf_cd_a[1];
    F[5][5] = -bcbf_cd_a[2];
}

static inline void bcbf_matmul6(double A[6][6], double B[6][6], double C[6][6])
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

static inline void bcbf_push_aug_row(const double grad0[6], const double z0[6],
                                     double hval, double alpha,
                                     double Arow[3], double *bval)
{
    double fhat0[6];
    bcbf_fhat_drift_aug(z0, fhat0);
    double Lf = 0.0;
    for (int j = 0; j < 6; j++) Lf += grad0[j] * fhat0[j];

    double Lg[3] = {
        grad0[3] * bcbf_cd_a[0],
        grad0[4] * bcbf_cd_a[1],
        grad0[5] * bcbf_cd_a[2]
    };

    Arow[0] = -Lg[0];
    Arow[1] = -Lg[1];
    Arow[2] = -Lg[2];
    *bval = Lf + alpha * hval;
}

static inline void bcbf_push_input_limit_rows(const double z0[6],
                                              double Arow[6][3], double bval[6])
{
    for (int j = 0; j < 3; j++) {
        // lower: h = u_j - u_min_j
        for (int k = 0; k < 3; k++) Arow[2*j][k] = 0.0;
        Arow[2*j][j] = -bcbf_cd_a[j];
        bval[2*j] = -bcbf_cd_a[j] * z0[3+j]
                    + bcbf_alpha_u * (z0[3+j] - bcbf_u_min[j]);

        // upper: h = u_max_j - u_j
        for (int k = 0; k < 3; k++) Arow[2*j+1][k] = 0.0;
        Arow[2*j+1][j] = bcbf_cd_a[j];
        bval[2*j+1] = bcbf_cd_a[j] * z0[3+j]
                      + bcbf_alpha_u * (bcbf_u_max[j] - z0[3+j]);
    }
}

static inline void bcbf_update_command_state(const double uhat[3])
{
    for (int j = 0; j < 3; j++) {
        double rho = exp(-bcbf_cd_a[j] * bcbf_cd_dt);
        bcbf_cd_u[j] = rho * bcbf_cd_u[j] + (1.0 - rho) * uhat[j];
    }
}

static inline void bcbf_nominal_surrogate(const double ud[3], double uhat_d[3])
{
    for (int j = 0; j < 3; j++) {
        double ud_dot = (ud[j] - bcbf_cd_ud_prev[j]) / fmax(bcbf_cd_dt, 1.0e-9);
        uhat_d[j] = bcbf_cd_u[j]
                    + (ud_dot / bcbf_cd_a[j])
                    + (bcbf_sigma0[j] / bcbf_cd_a[j]) * (ud[j] - bcbf_cd_u[j]);
    }
}

static inline void bcbf_backup_controller(double x, double y, double th, double ub[3])
{
    // Zero-action backup: commanded stop in the high-level kinematic model.
    // Under this controller, phi_b(tau,x)=x and the sensitivity is identity.
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
    return dx * dx + dy * dy - bcbf_obs[i].r * bcbf_obs[i].r; // safe iff >= 0
}

static inline double bcbf_h_backup_obs(int i, double x, double y)
{
    double dx = x - bcbf_obs[i].cx;
    double dy = y - bcbf_obs[i].cy;
    double rb = bcbf_obs[i].r + bcbf_backup_margin;
    return dx * dx + dy * dy - rb * rb; // backup set iff >= 0 for every obstacle
}


static inline void bcbf_h_backup_lse_and_gradx(double z[3], double Phi[3][3],
                                               double *h_lse, double grad_x[3])
{
    // Soft-min over terminal backup-set margins:
    // h_B,LSE = -(1/kappa) log(sum_i exp(-kappa h_{b,i})).
    // grad h_B,LSE = sum_i w_i grad h_{b,i}.
    const double kappa = bcbf_kappa_lse_B;
    double hvals[BCBF_N_OBS];
    double zmax = -1.0e300;

    for (int i = 0; i < BCBF_N_OBS; i++) {
        hvals[i] = bcbf_h_backup_obs(i, z[0], z[1]);
        double zi = -kappa * hvals[i];
        if (zi > zmax) zmax = zi;
    }

    double sumexp = 0.0;
    double grad_z_lse[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < BCBF_N_OBS; i++) {
        double wi = exp((-kappa * hvals[i]) - zmax);
        sumexp += wi;
        double dx = z[0] - bcbf_obs[i].cx;
        double dy = z[1] - bcbf_obs[i].cy;
        grad_z_lse[0] += wi * 2.0 * dx;
        grad_z_lse[1] += wi * 2.0 * dy;
    }

    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int i = 1; i < BCBF_N_OBS; i++) if (hvals[i] < hvals[imin]) imin = i;
        *h_lse = hvals[imin];
        double dx = z[0] - bcbf_obs[imin].cx;
        double dy = z[1] - bcbf_obs[imin].cy;
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

static inline void bcbf_push_constraint_row(double grad_x[3], double theta,
                                            double hval, double alpha,
                                            double Arow[3], double *bval)
{
    double G[3][3];
    bcbf_input_matrix(theta, G);

    // gu = grad_x * G
    double gu[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) gu[j] += grad_x[k] * G[k][j];
    }

    // CBF: gu*u + alpha*h >= 0  ->  -gu*u <= alpha*h
    Arow[0] = -gu[0];
    Arow[1] = -gu[1];
    Arow[2] = -gu[2];
    *bval = alpha * hval;
}

static void bcbf_build_constraints(double x, double y, double theta,
                                   const double u_state[3],
                                   double A_cbf[BCBF_N_CBF_ROWS][3],
                                   double b_cbf[BCBF_N_CBF_ROWS])
{
    double z0[6] = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double z[6]  = {x, y, theta, u_state[0], u_state[1], u_state[2]};
    double Phi[6][6] = {{0}};
    for (int i = 0; i < 6; i++) Phi[i][i] = 1.0;
    double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);

    double h_terms[BCBF_N_LSE_TERMS];
    double grad_terms[BCBF_N_LSE_TERMS][6];
    int term = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        for (int i = 0; i < BCBF_N_OBS; i++) {
            double dx = z[0] - bcbf_obs[i].cx;
            double dy = z[1] - bcbf_obs[i].cy;
            double grad_z[6] = {2.0 * dx, 2.0 * dy, 0.0, 0.0, 0.0, 0.0};
            h_terms[term] = bcbf_h_obs(i, z[0], z[1]);
            for (int j = 0; j < 6; j++) {
                grad_terms[term][j] = 0.0;
                for (int m = 0; m < 6; m++) grad_terms[term][j] += grad_z[m] * Phi[m][j];
            }
            term++;
        }
        if (k < BCBF_N_SAMPLES - 1) {
            double fz[6], F[6][6], Phidot[6][6];
            bcbf_fhat_drift_aug(z, fz);
            bcbf_backup_jacobian_aug(z, F);
            bcbf_matmul6(F, Phi, Phidot);
            for (int i = 0; i < 6; i++) z[i] += dt * fz[i];
            z[2] = bcbf_wrap_pi(z[2]);
            for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) Phi[i][j] += dt * Phidot[i][j];
        }
    }

    for (int i = 0; i < BCBF_N_OBS; i++) {
        double dx = z[0] - bcbf_obs[i].cx;
        double dy = z[1] - bcbf_obs[i].cy;
        double grad_z[6] = {2.0 * dx, 2.0 * dy, 0.0, 0.0, 0.0, 0.0};
        h_terms[term] = bcbf_h_backup_obs(i, z[0], z[1]);
        for (int j = 0; j < 6; j++) {
            grad_terms[term][j] = 0.0;
            for (int m = 0; m < 6; m++) grad_terms[term][j] += grad_z[m] * Phi[m][j];
        }
        term++;
    }

    // Current-time input-limit terms are composed into the same LSE barrier.
    for (int j = 0; j < 3; j++) {
        h_terms[term] = z0[3+j] - bcbf_u_min[j];
        for (int k = 0; k < 6; k++) grad_terms[term][k] = 0.0;
        grad_terms[term][3+j] = 1.0;
        term++;

        h_terms[term] = bcbf_u_max[j] - z0[3+j];
        for (int k = 0; k < 6; k++) grad_terms[term][k] = 0.0;
        grad_terms[term][3+j] = -1.0;
        term++;
    }

    const double kappa = bcbf_kappa_lse;
    double zmax = -1.0e300;
    for (int q = 0; q < term; q++) {
        double zq = -kappa * h_terms[q];
        if (zq > zmax) zmax = zq;
    }

    double sumexp = 0.0;
    double grad_lse[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for (int q = 0; q < term; q++) {
        double w = exp((-kappa * h_terms[q]) - zmax);
        sumexp += w;
        for (int j = 0; j < 6; j++) grad_lse[j] += w * grad_terms[q][j];
    }

    double H_lse;
    if (sumexp <= 0.0 || !isfinite(sumexp)) {
        int imin = 0;
        for (int q = 1; q < term; q++) if (h_terms[q] < h_terms[imin]) imin = q;
        H_lse = h_terms[imin];
        for (int j = 0; j < 6; j++) grad_lse[j] = grad_terms[imin][j];
    } else {
        H_lse = -(log(sumexp) + zmax) / kappa;
        for (int j = 0; j < 6; j++) grad_lse[j] /= sumexp;
    }

    bcbf_push_aug_row(grad_lse, z0, H_lse, bcbf_alpha_lse,
                      A_cbf[0], &b_cbf[0]);
}

// ------------------------------ OSQP state ----------------------------------// ------------------------------ OSQP state ----------------------------------
static OSQPSolver   *bcbf_solver = NULL;
static OSQPSettings  bcbf_settings;

static OSQPFloat bcbf_P_x[3];
static OSQPInt   bcbf_P_i[3] = {0,1,2};
static OSQPInt   bcbf_P_p[4] = {0,1,2,3};
static OSQPCscMatrix bcbf_P_mat;

static OSQPFloat bcbf_A_x[BCBF_N_ROWS * 3];
static OSQPInt   bcbf_A_i[BCBF_N_ROWS * 3];
static OSQPInt   bcbf_A_p[4];
static OSQPCscMatrix bcbf_A_mat;

static OSQPFloat bcbf_q[3];
static OSQPFloat bcbf_l[BCBF_N_ROWS];
static OSQPFloat bcbf_u[BCBF_N_ROWS];

static int bcbf_init_solver(double wx, double wy, double ww)
{
    bcbf_P_x[0] = (OSQPFloat)wx;
    bcbf_P_x[1] = (OSQPFloat)wy;
    bcbf_P_x[2] = (OSQPFloat)ww;

    bcbf_P_mat.m = 3; bcbf_P_mat.n = 3; bcbf_P_mat.nz = -1;
    bcbf_P_mat.nzmax = 3; bcbf_P_mat.x = bcbf_P_x;
    bcbf_P_mat.i = bcbf_P_i; bcbf_P_mat.p = bcbf_P_p;

    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        bcbf_A_p[col] = nnz;
        for (int row = 0; row < BCBF_N_ROWS; row++) {
            bcbf_A_x[nnz] = 0.0;
            bcbf_A_i[nnz] = row;
            nnz++;
        }
    }
    bcbf_A_p[3] = nnz;

    bcbf_A_mat.m = BCBF_N_ROWS; bcbf_A_mat.n = 3; bcbf_A_mat.nz = -1;
    bcbf_A_mat.nzmax = nnz; bcbf_A_mat.x = bcbf_A_x;
    bcbf_A_mat.i = bcbf_A_i; bcbf_A_mat.p = bcbf_A_p;

    for (int i = 0; i < BCBF_N_ROWS; i++) {
        bcbf_l[i] = (OSQPFloat)(-OSQP_INFTY);
        bcbf_u[i] = 0.0;
    }

    osqp_set_default_settings(&bcbf_settings);
    bcbf_settings.verbose = 0;
    bcbf_settings.warm_starting = 1;
    bcbf_settings.eps_abs = 1e-4;
    bcbf_settings.eps_rel = 1e-4;
    bcbf_settings.max_iter = 10000;

    OSQPInt flag = osqp_setup(&bcbf_solver, &bcbf_P_mat, bcbf_q,
                              &bcbf_A_mat, bcbf_l, bcbf_u,
                              BCBF_N_ROWS, 3, &bcbf_settings);
    if (flag != 0) {
        fprintf(stderr, "[bCBF] OSQP setup failed: %d\n", (int)flag);
        return 0;
    }
    return 1;
}

// =============================================================================
// Public filter function. Same signature as the original CBF filter.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    const double cbf_filter_t0 = cbf_now_sec();
    cbf_log_filter_time_sec = 0.0;
    double ud[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    double uhat_d[3];
    double A_cbf[BCBF_N_CBF_ROWS][3];
    double b_cbf[BCBF_N_CBF_ROWS];

    cbf_log_u_nom[0] = ud[0];
    cbf_log_u_nom[1] = ud[1];
    cbf_log_u_nom[2] = ud[2];

    double ub_log[3];
    bcbf_backup_controller(x, y, theta, ub_log);
    cbf_log_u_backup[0] = ub_log[0];
    cbf_log_u_backup[1] = ub_log[1];
    cbf_log_u_backup[2] = ub_log[2];

    if (!bcbf_cd_init) {
        for (int j = 0; j < 3; j++) {
            bcbf_cd_u[j] = bcbf_clip(ud[j], bcbf_u_min[j], bcbf_u_max[j]);
            bcbf_cd_ud_prev[j] = ud[j];
        }
        bcbf_cd_init = 1;
    }

    cbf_log_mu = -1.0;
    cbf_log_status = 0;

    bcbf_nominal_surrogate(ud, uhat_d);
    bcbf_build_constraints(x, y, theta, bcbf_cd_u, A_cbf, b_cbf);

    if (bcbf_solver == NULL) {
        if (!bcbf_init_solver(wx, wy, ww)) {
            cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
            return 0;
        }
    }

    bcbf_q[0] = (OSQPFloat)(-wx * uhat_d[0]);
    bcbf_q[1] = (OSQPFloat)(-wy * uhat_d[1]);
    bcbf_q[2] = (OSQPFloat)(-ww * uhat_d[2]);

    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < BCBF_N_ROWS; row++) {
            bcbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        }
    }

    for (int i = 0; i < BCBF_N_ROWS; i++) bcbf_u[i] = (OSQPFloat)b_cbf[i];

    osqp_update_data_vec(bcbf_solver, bcbf_q, bcbf_l, bcbf_u);
    osqp_update_data_mat(bcbf_solver, NULL, NULL, 0, bcbf_A_x, NULL, bcbf_A_mat.nzmax);
    double qp_t0 = cbf_now_sec();
    osqp_solve(bcbf_solver);
    cbf_log_qp_time_sec = cbf_now_sec() - qp_t0;

    int status_val = bcbf_solver->info->status_val;
    int status_ok = (status_val == OSQP_SOLVED || status_val == OSQP_SOLVED_INACCURATE);
    if (status_ok) qp_success_count++; else qp_fail_count++;
    if ((qp_success_count + qp_fail_count) % 100 == 0) {
        printf("[bCBF-CD] success=%d fail=%d\n", qp_success_count, qp_fail_count);
    }

    double uhat[3];
    if (status_ok) {
        uhat[0] = (double)bcbf_solver->solution->x[0];
        uhat[1] = (double)bcbf_solver->solution->x[1];
        uhat[2] = (double)bcbf_solver->solution->x[2];
        cbf_log_status = 1;
    } else {
        fprintf(stderr, "[bCBF-CD] OSQP failed (%s). Falling back to backup surrogate.\n",
                bcbf_solver->info->status);
        uhat[0] = ub_log[0];
        uhat[1] = ub_log[1];
        uhat[2] = ub_log[2];
        cbf_log_status = 0;
    }

    bcbf_update_command_state(uhat);

    *vx_cmd = bcbf_cd_u[0];
    *vy_cmd = bcbf_cd_u[1];
    *wz_cmd = bcbf_cd_u[2];

    cbf_log_u_safe[0] = *vx_cmd;
    cbf_log_u_safe[1] = *vy_cmd;
    cbf_log_u_safe[2] = *wz_cmd;
    cbf_log_mu = -1.0;

    for (int j = 0; j < 3; j++) bcbf_cd_ud_prev[j] = ud[j];
    cbf_log_filter_time_sec = cbf_now_sec() - cbf_filter_t0;
    return 1;
}
