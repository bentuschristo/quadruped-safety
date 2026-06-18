#pragma once
#include <math.h>
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

// =============================================================================
// Original Backup-CBF filter for the high-level quadruped model
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
// Backup-set convention:
//      C_B = { X : h_b(X) >= 0 }
//      h_b(X) = y - y_B, with y_B = 1.5
//
// In words: the backup set is the upper half-plane y >= 2.
// The backup controller drives the robot upward in world coordinates.
//
// The bCBF-QP enforces one constraint for every obstacle at every sampled point
// along the backup rollout, plus one terminal backup-set constraint.
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
    { -2.0,  1.0, 0.45, 1.0, 1 },
    { -4.0, -1.0, 0.45, 1.0, 1 },
    { -2.0, -1.0, 0.45, 1.0, 1 },
    { -4.0,  1.0, 0.45, 1.0, 1 }
};

// static BCBFObstacle bcbf_obs[BCBF_N_OBS] = {
//     { -2.0,  1.0, 0.42, 1.0, 1 },
//     { -4.0, -1.0, 0.42, 1.0, 1 },
//     { -2.0, -1.0, 0.42, 1.0, 1 },
//     { -4.0,  1.0, 0.42, 1.0, 1 },
//     { -2.8,  0.0, 0.30, 1.0, 1 },
//     { -3.3,  0.0, 0.30, 1.0, 1 }
// };

// Input limits. Keep these consistent with the rest of the high-level controller.
static double bcbf_u_max[3] = {  1.0,  0.3,  1.0 };
static double bcbf_u_min[3] = { -1.0, -0.3, -1.0 };

// Backup set: upper refuge half-plane y >= BACKUP_Y.
// This avoids forcing the backup controller to drive through the obstacle
// cluster to reach a point. Instead, the backup controller escapes upward
// above the obstacle corridor.
static const double bcbf_backup_y = 1.5;
static const double bcbf_alpha_B = 0.8;

// Backup controller: world-frame upward escape command.
// Desired world velocity is [0, v_up]. This is converted to the body-frame
// command [vx, vy] using R(theta)^T.
static const double bcbf_ky_backup = 0.8;
static const double bcbf_vup_max   = 0.55;
static const double bcbf_kw_backup = 0.0;

// Backup rollout discretization. Longer horizon helps the backup controller
// certify recovery to the right-side backup set from states left of the origin.
#define BCBF_N_SAMPLES 20
static const double bcbf_T_backup = 4.0;

// Total constraints = obstacle constraints at each backup sample + terminal hb + boxes
#define BCBF_N_CBF_ROWS (BCBF_N_OBS * BCBF_N_SAMPLES + 1)
#define BCBF_N_ROWS     (BCBF_N_CBF_ROWS + 3)

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

static inline void bcbf_backup_controller(double x, double y, double th, double ub[3])
{
    (void)x;

    // Drive upward in world coordinates toward the backup half-plane y >= 2.
    // If already above the backup line, do not intentionally drive downward.
    double v_up = bcbf_ky_backup * (bcbf_backup_y - y);
    if (v_up < 0.0) v_up = 0.0;
    v_up = bcbf_clip(v_up, 0.0, bcbf_vup_max);

    // Convert world velocity [0, v_up] into body-frame [vx, vy]:
    // [vx; vy] = R(theta)^T [0; v_up] = [sin(theta) v_up; cos(theta) v_up].
    ub[0] = bcbf_clip(v_up * sin(th), bcbf_u_min[0], bcbf_u_max[0]);
    ub[1] = bcbf_clip(v_up * cos(th), bcbf_u_min[1], bcbf_u_max[1]);

    // Keep heading unchanged by default. You can add heading stabilization here
    // later if needed.
    ub[2] = bcbf_clip(bcbf_kw_backup * 0.0, bcbf_u_min[2], bcbf_u_max[2]);
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

static inline double bcbf_h_backup(double x, double y)
{
    (void)x;
    return y - bcbf_backup_y; // backup set y >= backup_y iff h_b >= 0
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
                                   double A_cbf[BCBF_N_CBF_ROWS][3],
                                   double b_cbf[BCBF_N_CBF_ROWS])
{
    double z[3] = {x, y, theta};
    double Phi[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double dt = bcbf_T_backup / (double)(BCBF_N_SAMPLES - 1);
    int row = 0;

    for (int k = 0; k < BCBF_N_SAMPLES; k++) {
        // Obstacle constraints at backup-flow sample z = phi_b(tau_k, x)
        for (int i = 0; i < BCBF_N_OBS; i++) {
            double dx = z[0] - bcbf_obs[i].cx;
            double dy = z[1] - bcbf_obs[i].cy;
            double grad_z[3] = {2.0 * dx, 2.0 * dy, 0.0};
            double grad_x[3] = {0.0, 0.0, 0.0};

            // grad_x = grad_z * Phi
            for (int j = 0; j < 3; j++) {
                for (int m = 0; m < 3; m++) grad_x[j] += grad_z[m] * Phi[m][j];
            }

            double hval = bcbf_h_obs(i, z[0], z[1]);
            bcbf_push_constraint_row(grad_x, theta, hval, bcbf_obs[i].alpha,
                                     A_cbf[row], &b_cbf[row]);
            row++;
        }

        // Propagate backup flow and sensitivity, except after last sample.
        if (k < BCBF_N_SAMPLES - 1) {
            double fz[3];
            double J[3][3];
            double Phidot[3][3];
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

    // Terminal backup-set constraint h_b(phi_b(T,x)) >= 0.
    // For C_B = {y >= backup_y}, h_b = y - backup_y and grad_z h_b = [0, 1, 0].
    double grad_zB[3] = {0.0, 1.0, 0.0};
    double grad_xB[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 3; j++) {
        for (int m = 0; m < 3; m++) grad_xB[j] += grad_zB[m] * Phi[m][j];
    }
    double hb = bcbf_h_backup(z[0], z[1]);
    bcbf_push_constraint_row(grad_xB, theta, hb, bcbf_alpha_B,
                             A_cbf[row], &b_cbf[row]);
}

// ------------------------------ OSQP state ----------------------------------
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
        for (int row = 0; row < BCBF_N_CBF_ROWS; row++) {
            bcbf_A_x[nnz] = 0.0;
            bcbf_A_i[nnz] = row;
            nnz++;
        }
        bcbf_A_x[nnz] = 1.0;
        bcbf_A_i[nnz] = BCBF_N_CBF_ROWS + col;
        nnz++;
    }
    bcbf_A_p[3] = nnz;

    bcbf_A_mat.m = BCBF_N_ROWS; bcbf_A_mat.n = 3; bcbf_A_mat.nz = -1;
    bcbf_A_mat.nzmax = nnz; bcbf_A_mat.x = bcbf_A_x;
    bcbf_A_mat.i = bcbf_A_i; bcbf_A_mat.p = bcbf_A_p;

    for (int i = 0; i < BCBF_N_CBF_ROWS; i++) {
        bcbf_l[i] = (OSQPFloat)(-OSQP_INFTY);
        bcbf_u[i] = 0.0;
    }
    for (int j = 0; j < 3; j++) {
        bcbf_l[BCBF_N_CBF_ROWS + j] = (OSQPFloat)bcbf_u_min[j];
        bcbf_u[BCBF_N_CBF_ROWS + j] = (OSQPFloat)bcbf_u_max[j];
    }

    osqp_set_default_settings(&bcbf_settings);
    bcbf_settings.verbose = 0;
    bcbf_settings.warm_starting = 1;
    bcbf_settings.eps_abs = 1e-4;
    bcbf_settings.eps_rel = 1e-4;
    bcbf_settings.max_iter = 1000;

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
    double u0[3] = {*vx_cmd, *vy_cmd, *wz_cmd};
    double A_cbf[BCBF_N_CBF_ROWS][3];
    double b_cbf[BCBF_N_CBF_ROWS];

    // Log nominal command before the QP modifies it
    cbf_log_u_nom[0] = u0[0];
    cbf_log_u_nom[1] = u0[1];
    cbf_log_u_nom[2] = u0[2];

    // Log backup command
    double ub_log[3];
    bcbf_backup_controller(x, y, theta, ub_log);

    cbf_log_u_backup[0] = ub_log[0];
    cbf_log_u_backup[1] = ub_log[1];
    cbf_log_u_backup[2] = ub_log[2];

    cbf_log_mu = -1.0;
    cbf_log_status = 0;

    bcbf_build_constraints(x, y, theta, A_cbf, b_cbf);

    if (bcbf_solver == NULL) {
        if (!bcbf_init_solver(wx, wy, ww)) return 0;
    }

    bcbf_q[0] = (OSQPFloat)(-wx * u0[0]);
    bcbf_q[1] = (OSQPFloat)(-wy * u0[1]);
    bcbf_q[2] = (OSQPFloat)(-ww * u0[2]);

    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < BCBF_N_CBF_ROWS; row++) {
            bcbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        }
        nnz++; // fixed identity element
    }

    for (int i = 0; i < BCBF_N_CBF_ROWS; i++) bcbf_u[i] = (OSQPFloat)b_cbf[i];

    osqp_update_data_vec(bcbf_solver, bcbf_q, bcbf_l, bcbf_u);
    osqp_update_data_mat(bcbf_solver, NULL, NULL, 0, bcbf_A_x, NULL, bcbf_A_mat.nzmax);
    
    double qp_t0 = cbf_now_sec();
    osqp_solve(bcbf_solver);
    cbf_log_qp_time_sec = cbf_now_sec() - qp_t0;

    int status_val = bcbf_solver->info->status_val;

    int status_ok =
        (status_val == OSQP_SOLVED ||
        status_val == OSQP_SOLVED_INACCURATE);

    if (status_ok) {
        qp_success_count++;
    } else {
        qp_fail_count++;
    }

    if ((qp_success_count + qp_fail_count) % 100 == 0) {
        printf("[bCBF] success=%d fail=%d\n", qp_success_count, qp_fail_count);
    }

    int modified = 0;

    if (status_ok) {
        *vx_cmd = (double)bcbf_solver->solution->x[0];
        *vy_cmd = (double)bcbf_solver->solution->x[1];
        *wz_cmd = (double)bcbf_solver->solution->x[2];
        modified = 1;

        // Log the safe command after the QP modifies it
        cbf_log_u_safe[0] = *vx_cmd;
        cbf_log_u_safe[1] = *vy_cmd;
        cbf_log_u_safe[2] = *wz_cmd;
        cbf_log_status = 1;
        cbf_log_mu = -1.0;
    } else {
        fprintf(stderr, "[bCBF] OSQP failed (%s). Falling back to backup control.\n",
                bcbf_solver->info->status);

        double ub[3];
        bcbf_backup_controller(x, y, theta, ub);

        *vx_cmd = bcbf_clip(ub[0], bcbf_u_min[0], bcbf_u_max[0]);
        *vy_cmd = bcbf_clip(ub[1], bcbf_u_min[1], bcbf_u_max[1]);
        *wz_cmd = bcbf_clip(ub[2], bcbf_u_min[2], bcbf_u_max[2]);

        // Log the fallback command as the safe/applied command.
        cbf_log_u_safe[0] = *vx_cmd;
        cbf_log_u_safe[1] = *vy_cmd;
        cbf_log_u_safe[2] = *wz_cmd;
        cbf_log_status = 0;
        cbf_log_mu = -1.0;

        // The command was still modified, even though the QP failed.
        modified = 1;
    }

    return modified;
}
