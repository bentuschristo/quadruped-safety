#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  Policy Control Barrier Function (PCBF) — Circle Obstacle Avoidance (v4_lse)
//
//  Uses plain Log-Sum-Exp (LSE) as a smooth approximation to max over the
//  discrete rollout values {y_0, ..., y_H}, with an explicit closed-form
//  gradient via softmax weights and the sensitivity matrix.
//
//  LSE Value (original, no numerical stabilization):
//    V_LSE = (1/beta) * log( sum_k exp(beta * y_k) )
//
//  LSE Gradient (explicit, closed-form):
//    dV/dx0 = sum_k  w_k * Phi_k^T * grad_x h(x_k)
//    where:
//      w_k  = exp(beta * y_k) / sum_j exp(beta * y_j)   (softmax weights)
//      Phi_k = d(x_k)/d(x0)                              (sensitivity matrix)
//      grad_x h = [-2(x_k-cx), -2(y_k-cy), 0]^T
//
//  NOTE: Plain LSE without y_max subtraction will overflow for large beta
//  or large |y_k|. This is the raw original form for study purposes.
//
//  Based on:
//    Knoedler et al., "Safety on the Fly: Constructing Robust Safety Filters
//    via Policy Control Barrier Functions at Runtime", IEEE RA-L 2025.
// =============================================================================


// ---- Obstacle definition ----------------------------------------------------
typedef struct {
    double cx;
    double cy;
    double r;
    double alpha;
    int    active;
} CircleCBF;

// ---- Edit obstacles here ----------------------------------------------------
#define N_OBS 4
const double buffer = 0.1;

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, 0.5, 1 },
    { -4.0, -1.0, 0.35, 0.5, 1 },
    { -2.0, -1.0, 0.35, 0.5, 1 },
    { -4.0,  1.0, 0.35, 0.5, 1 },
};
// -----------------------------------------------------------------------------


// ---- Rollout parameters -----------------------------------------------------
#define PCBF_H      20
#define PCBF_DT     0.1

#define POLICY_ZERO_VELOCITY  0
#define POLICY_NOMINAL        1
#define POLICY_PUSH_VX        2

static int    pcbf_policy = POLICY_PUSH_VX;
static double vx_push     = 0.5;
// -----------------------------------------------------------------------------


// ---- LSE parameter ----------------------------------------------------------
//  beta: sharpness of the LSE approximation.
//    beta -> inf  : approaches exact discrete max  (may overflow)
//    beta small   : smooth average of all y_k      (underestimates max)
static double lse_beta = 20.0;
// -----------------------------------------------------------------------------


// ---- Input limits -----------------------------------------------------------
static double u_max[3] = {  1.0,  0.3,  1.0 };
static double u_min[3] = { -1.0, -0.3, -1.0 };

// ---- Slack variable penalty -------------------------------------------------
static double p_slack = 1e4;
// -----------------------------------------------------------------------------


// ---- OSQP persistent state --------------------------------------------------
static OSQPSolver    *pcbf_solver   = NULL;
static OSQPSettings   pcbf_settings;

#define N_VARS  4
#define N_ROWS  (N_OBS + 3 + 1)

static OSQPFloat pcbf_P_x[N_VARS];
static OSQPInt   pcbf_P_i[N_VARS];
static OSQPInt   pcbf_P_p[N_VARS+1];
static OSQPCscMatrix pcbf_P_mat;

static OSQPFloat pcbf_A_x[N_ROWS * N_VARS];
static OSQPInt   pcbf_A_i[N_ROWS * N_VARS];
static OSQPInt   pcbf_A_p[N_VARS+1];
static OSQPCscMatrix pcbf_A_mat;

static OSQPFloat pcbf_q[N_VARS];
static OSQPFloat pcbf_l[N_ROWS];
static OSQPFloat pcbf_u[N_ROWS];
// -----------------------------------------------------------------------------


// =============================================================================
//  lse_compute()
//
//  Plain Log-Sum-Exp — no numerical stabilization.
//
//  Value:
//    V_LSE = (1/beta) * log( sum_k exp(beta * y[k]) )
//
//  Softmax weights (= dV_LSE/dy_k, used directly by chain rule):
//    w_k = exp(beta * y[k]) / sum_j exp(beta * y[j])
//
//  Arguments:
//    y[0..n-1]    : discrete h values along rollout
//    n            : number of points
//    beta         : LSE sharpness (>0)
//    w_out[0..n-1]: output softmax weights (sum to 1)
//
//  Returns: V_LSE value
// =============================================================================
static double lse_compute(const double *y, int n, double beta, double *w_out)
{
    // Step 1: compute exp(beta * y_k) for each k and sum
    double sum_exp = 0.0;
    for (int k = 0; k < n; k++) {
        w_out[k] = exp(beta * y[k]);
        sum_exp += w_out[k];
    }

    // Step 2: normalize -> softmax weights  w_k = exp(beta*y_k) / sum_exp
    for (int k = 0; k < n; k++)
        w_out[k] /= sum_exp;

    // Step 3: LSE value = (1/beta) * log(sum_exp)
    return log(sum_exp) / beta;
}


// =============================================================================
//  sensitivity_step()
//
//  Propagates the 3x3 sensitivity matrix one Euler step forward:
//    Phi_{k+1} = (I + dt * J_k) * Phi_k
//
//  Unicycle Jacobian J_k = df/dx at step k:
//    row 0: [0,  0,  -sin(th)*ux - cos(th)*uy]
//    row 1: [0,  0,   cos(th)*ux - sin(th)*uy]
//    row 2: [0,  0,   0                       ]
//
//  Phi stored row-major: Phi[row*3 + col]
// =============================================================================
static void sensitivity_step(
    const double *Phi_in,
    double       *Phi_out,
    double        theta_k,
    double        ux,
    double        uy)
{
    double J02 = -sin(theta_k)*ux - cos(theta_k)*uy;
    double J12 =  cos(theta_k)*ux - sin(theta_k)*uy;

    for (int c = 0; c < 3; c++) {
        Phi_out[0*3 + c] = Phi_in[0*3 + c] + PCBF_DT * J02 * Phi_in[2*3 + c];
        Phi_out[1*3 + c] = Phi_in[1*3 + c] + PCBF_DT * J12 * Phi_in[2*3 + c];
        Phi_out[2*3 + c] = Phi_in[2*3 + c];
    }
}


// =============================================================================
//  pcbf_rollout_and_value()
//
//  Rolls out the design policy for H steps. For each obstacle computes:
//
//  V_LSE = (1/beta) * log( sum_k exp(beta * h_i(x_k)) )
//
//  dV/dx0 = sum_k  w_k * Phi_k^T * grad_x h_i(x_k)
//
//  where w_k are the softmax weights from lse_compute(), Phi_k is the
//  sensitivity matrix d(x_k)/d(x0), and grad_x h = [-2(x-cx), -2(y-cy), 0].
// =============================================================================
static void pcbf_rollout_and_value(
    double x0, double y0, double theta0,
    double vx_nom, double vy_nom, double wz_nom,
    double V_val[N_OBS],
    double dV_dx0[N_OBS][3])
{
    // ---- Design policy input ------------------------------------------------
    double ux, uy, uw;
    if (pcbf_policy == POLICY_ZERO_VELOCITY) {
        ux = 0.0;  uy = 0.0;  uw = 0.0;
    } else if (pcbf_policy == POLICY_NOMINAL) {
        ux = vx_nom;  uy = vy_nom;  uw = wz_nom;
    } else {
        ux = vx_push;  uy = 0.0;  uw = 0.0;
    }

    // ---- Storage ------------------------------------------------------------
    double xs[PCBF_H+1], ys[PCBF_H+1], ths[PCBF_H+1];
    double h_traj[N_OBS][PCBF_H+1];
    double Phi[PCBF_H+1][9];
    double w[PCBF_H+1];

    // ---- Initialize ---------------------------------------------------------
    xs[0]  = x0;
    ys[0]  = y0;
    ths[0] = theta0;

    memset(Phi[0], 0, sizeof(Phi[0]));
    Phi[0][0*3+0] = 1.0;
    Phi[0][1*3+1] = 1.0;
    Phi[0][2*3+2] = 1.0;

    for (int i = 0; i < N_OBS; i++) {
        if (!obstacles[i].active) { h_traj[i][0] = -1e9; continue; }
        double dx = x0 - obstacles[i].cx;
        double dy = y0 - obstacles[i].cy;
        h_traj[i][0] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
    }

    // ---- Forward rollout ----------------------------------------------------
    for (int k = 0; k < PCBF_H; k++)
    {
        double ck = cos(ths[k]);
        double sk = sin(ths[k]);

        xs[k+1]  = xs[k]  + PCBF_DT * (ck*ux - sk*uy);
        ys[k+1]  = ys[k]  + PCBF_DT * (sk*ux + ck*uy);
        ths[k+1] = ths[k] + PCBF_DT * uw;

        sensitivity_step(Phi[k], Phi[k+1], ths[k], ux, uy);

        for (int i = 0; i < N_OBS; i++) {
            if (!obstacles[i].active) { h_traj[i][k+1] = -1e9; continue; }
            double dx = xs[k+1] - obstacles[i].cx;
            double dy = ys[k+1] - obstacles[i].cy;
            h_traj[i][k+1] = obstacles[i].r*obstacles[i].r - dx*dx - dy*dy;
        }
    }

    int n_pts = PCBF_H + 1;

    // ---- Per obstacle: LSE value + explicit gradient ------------------------
    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) {
            V_val[i]     = -1e9;
            dV_dx0[i][0] = 0.0;
            dV_dx0[i][1] = 0.0;
            dV_dx0[i][2] = 0.0;
            continue;
        }

        // V_LSE and softmax weights
        V_val[i] = lse_compute(h_traj[i], n_pts, lse_beta, w);

        // Explicit gradient: dV/dx0 = sum_k  w_k * Phi_k^T * grad_x h(x_k)
        double dV[3] = {0.0, 0.0, 0.0};

        for (int k = 0; k < n_pts; k++)
        {
            double gh_x  = -2.0 * (xs[k] - obstacles[i].cx);
            double gh_y  = -2.0 * (ys[k] - obstacles[i].cy);
            double gh_th = 0.0;

            for (int c = 0; c < 3; c++) {
                dV[c] += w[k] * (
                    Phi[k][0*3 + c] * gh_x  +
                    Phi[k][1*3 + c] * gh_y  +
                    Phi[k][2*3 + c] * gh_th
                );
            }
        }

        dV_dx0[i][0] = dV[0];
        dV_dx0[i][1] = dV[1];
        dV_dx0[i][2] = dV[2];
    }
}


// ---- OSQP setup (called once on first use) ----------------------------------
static int pcbf_init_solver(double wx, double wy, double ww, int nc)
{
    int nr = nc + 3 + 1;

    pcbf_P_x[0] = (OSQPFloat)wx;
    pcbf_P_x[1] = (OSQPFloat)wy;
    pcbf_P_x[2] = (OSQPFloat)ww;
    pcbf_P_x[3] = (OSQPFloat)p_slack;
    for (int j = 0; j < N_VARS; j++) { pcbf_P_i[j] = j; pcbf_P_p[j] = j; }
    pcbf_P_p[N_VARS] = N_VARS;

    pcbf_P_mat.m = N_VARS; pcbf_P_mat.n = N_VARS;
    pcbf_P_mat.nz = -1;    pcbf_P_mat.nzmax = N_VARS;
    pcbf_P_mat.x = pcbf_P_x; pcbf_P_mat.i = pcbf_P_i; pcbf_P_mat.p = pcbf_P_p;

    int nnz = 0;

    pcbf_A_p[0] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 0; nnz++;

    pcbf_A_p[1] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 1; nnz++;

    pcbf_A_p[2] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = 0.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = 1.0f; pcbf_A_i[nnz] = nc + 2; nnz++;

    pcbf_A_p[3] = nnz;
    for (int row = 0; row < nc; row++) {
        pcbf_A_x[nnz] = -1.0f; pcbf_A_i[nnz] = row; nnz++;
    }
    pcbf_A_x[nnz] = -1.0f; pcbf_A_i[nnz] = nc + 3; nnz++;

    pcbf_A_p[4] = nnz;

    pcbf_A_mat.m = nr; pcbf_A_mat.n = N_VARS;
    pcbf_A_mat.nz = -1; pcbf_A_mat.nzmax = nnz;
    pcbf_A_mat.x = pcbf_A_x; pcbf_A_mat.i = pcbf_A_i; pcbf_A_mat.p = pcbf_A_p;

    for (int j = 0; j < N_VARS; j++) pcbf_q[j] = 0.0f;
    for (int i = 0; i < nc; i++) {
        pcbf_l[i] = (OSQPFloat)(-OSQP_INFTY);
        pcbf_u[i] = 0.0f;
    }
    for (int j = 0; j < 3; j++) {
        pcbf_l[nc + j] = (OSQPFloat)u_min[j];
        pcbf_u[nc + j] = (OSQPFloat)u_max[j];
    }
    pcbf_l[nc + 3] = (OSQPFloat)(-OSQP_INFTY);
    pcbf_u[nc + 3] = 0.0f;

    osqp_set_default_settings(&pcbf_settings);
    pcbf_settings.verbose       = 0;
    pcbf_settings.warm_starting = 1;
    pcbf_settings.eps_abs       = 1e-3;
    pcbf_settings.eps_rel       = 1e-3;
    pcbf_settings.max_iter      = 1000;

    OSQPInt flag = osqp_setup(&pcbf_solver, &pcbf_P_mat, pcbf_q,
                               &pcbf_A_mat, pcbf_l, pcbf_u, nr, N_VARS, &pcbf_settings);
    if (flag != 0) {
        fprintf(stderr, "[PCBF_LSE] OSQP setup failed (code %d)\n", (int)flag);
        return 0;
    }
    return 1;
}


// =============================================================================
//  cbf_circle_obstacles_filter() — drop-in replacement for v4 (spline)
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    // STEP 1: rollout + LSE value + explicit gradient
    double V_val[N_OBS];
    double dV_dx0[N_OBS][3];

    pcbf_rollout_and_value(x, y, theta,
                           *vx_cmd, *vy_cmd, *wz_cmd,
                           V_val, dV_dx0);

    // STEP 2: build QP rows
    //   dV/dx * x_dot + alpha * V <= 0
    //   x_dot = [c*vx - s*vy,  s*vx + c*vy,  wz]
    const double c = cos(theta);
    const double s = sin(theta);

    double A_cbf[N_OBS][3];
    double b_cbf[N_OBS];
    int nc = 0;

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) continue;

        A_cbf[nc][0] =  dV_dx0[i][0]*c + dV_dx0[i][1]*s;
        A_cbf[nc][1] = -dV_dx0[i][0]*s + dV_dx0[i][1]*c;
        A_cbf[nc][2] =  dV_dx0[i][2];
        b_cbf[nc]    = -obstacles[i].alpha * V_val[i] - buffer;

        nc++;
    }

    if (nc == 0) return 0;

    // STEP 3: init OSQP on first call
    if (pcbf_solver == NULL) {
        if (!pcbf_init_solver(wx, wy, ww, nc)) return 0;
    }

    // STEP 4: update OSQP data
    pcbf_q[0] = (OSQPFloat)(-wx * (*vx_cmd));
    pcbf_q[1] = (OSQPFloat)(-wy * (*vy_cmd));
    pcbf_q[2] = (OSQPFloat)(-ww * (*wz_cmd));
    pcbf_q[3] = 0.0f;

    int nnz = 0;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < nc; row++)
            pcbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        nnz++;
    }

    for (int i = 0; i < nc; i++)
        pcbf_u[i] = (OSQPFloat)b_cbf[i];

    osqp_update_data_vec(pcbf_solver, pcbf_q, pcbf_l, pcbf_u);
    osqp_update_data_mat(pcbf_solver, NULL, NULL, 0, pcbf_A_x, NULL, pcbf_A_mat.nzmax);

    // STEP 5: solve and extract
    osqp_solve(pcbf_solver);

    int modified = 0;

    if (pcbf_solver->info->status_val == OSQP_SOLVED ||
        pcbf_solver->info->status_val == OSQP_SOLVED_INACCURATE)
    {
        *vx_cmd = (double)pcbf_solver->solution->x[0];
        *vy_cmd = (double)pcbf_solver->solution->x[1];
        *wz_cmd = (double)pcbf_solver->solution->x[2];
        double delta = (double)pcbf_solver->solution->x[3];
        if (delta > 1e-3)
            fprintf(stderr, "[PCBF_LSE] slack active: delta=%.4f (CBF softened)\n", delta);
        modified = 1;
    }
    else
    {
        fprintf(stderr, "[PCBF_LSE] OSQP did not converge (%s), clamping to limits\n",
                pcbf_solver->info->status);
        *vx_cmd = fmax(u_min[0], fmin(u_max[0], *vx_cmd));
        *vy_cmd = fmax(u_min[1], fmin(u_max[1], *vy_cmd));
        *wz_cmd = fmax(u_min[2], fmin(u_max[2], *wz_cmd));
    }

    return modified;
}
