#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "osqp.h"

// =============================================================================
//  CBF Circular Obstacle Avoidance  —  OSQP version
//
//  KEY DESIGN: the OSQP solver is initialized ONCE (first call) as a static
//  variable. Every subsequent call only updates the constraint data via
//  osqp_update_data_vectors() — no memory allocation, no re-factorization.
//  This makes it safe to call at 1kHz inside the control loop.
//
//  THE QP solved each step:
//    minimize    0.5*(u - u0)^T * W * (u - u0)    (stay close to nominal)
//    subject to  A * u <= b                        (CBF safety constraints)
//
//  where u = [vx, vy, wz],  W = diag(wx, wy, ww),  one row in A per obstacle.
// =============================================================================


// ---- Obstacle definition ----------------------------------------------------
typedef struct {
    double cx;     // obstacle center x (world frame)
    double cy;     // obstacle center y (world frame)
    double r;      // safety radius
    double alpha;  // CBF gain: larger = robot reacts sooner
    int    active; // 1 = on,  0 = off
} CircleCBF;

// ---- Edit obstacle positions and sizes here ---------------------------------
#define N_OBS 4

static CircleCBF obstacles[N_OBS] = {
    { -2.0,  1.0, 0.35, 1.0, 1 },
    { -4.0, -1.0, 0.35, 1.0, 1 },
    { -2.0, -1.0, 0.35, 1.0, 1 },
    { -4.0,  1.0, 0.35, 1.0, 1 },
};

// static CircleCBF obstacles[N_OBS] = {
//     { -2.0,  1.0, 0.35, 0.5, 1 },
//     { -4.0, -1.0, 0.35, 0.5, 1 }
// };

// ---- Input limits (keep in sync with control_params.h) ---------------------
// Added as box constraints inside the QP so the solver always respects them,
// even during obstacle avoidance corrections.
// [ vx,   vy,   wz  ]
static double u_max[3] = {  1.0,  0.3,  1.0 };
static double u_min[3] = { -1.0, -0.3, -1.0 };
// -----------------------------------------------------------------------------


// ---- Static solver state (persists across calls) ----------------------------
static OSQPSolver   *cbf_solver   = NULL;
static OSQPSettings  cbf_settings;

// Total rows = N_OBS (CBF) + 3 (box constraints on vx, vy, wz)
#define N_ROWS (N_OBS + 3)

// Fixed-size storage for P and A matrices
static OSQPFloat cbf_P_x[3];
static OSQPInt   cbf_P_i[3] = {0,1,2};
static OSQPInt   cbf_P_p[4] = {0,1,2,3};
static OSQPCscMatrix cbf_P_mat;

static OSQPFloat cbf_A_x[N_ROWS * 3]; // constraint matrix values (col-major)
static OSQPInt   cbf_A_i[N_ROWS * 3]; // row indices
static OSQPInt   cbf_A_p[4];           // column pointers
static OSQPCscMatrix cbf_A_mat;

static OSQPFloat cbf_q[3];             // linear cost: -W * u0
static OSQPFloat cbf_l[N_ROWS];        // lower bounds
static OSQPFloat cbf_u[N_ROWS];        // upper bounds
// -----------------------------------------------------------------------------


// ---- Initialize solver (called once on first use) ---------------------------
static int cbf_init_solver(double wx, double wy, double ww, int nc)
{
    int nr = nc + 3; // total rows: nc CBF rows + 3 box constraint rows

    // P matrix: diagonal cost W = diag(wx, wy, ww)
    cbf_P_x[0] = (OSQPFloat)wx;
    cbf_P_x[1] = (OSQPFloat)wy;
    cbf_P_x[2] = (OSQPFloat)ww;

    cbf_P_mat.m     = 3;
    cbf_P_mat.n     = 3;
    cbf_P_mat.nz    = -1;
    cbf_P_mat.nzmax = 3;
    cbf_P_mat.x     = cbf_P_x;
    cbf_P_mat.i     = cbf_P_i;
    cbf_P_mat.p     = cbf_P_p;

    // A matrix: (nc+3) x 3, column-major CSC
    // Rows 0..nc-1   : CBF constraints (A_cbf rows, filled with zeros for init)
    // Rows nc..nc+2  : identity block  (box constraints on vx, vy, wz)
    //
    // Column j of A contains:
    //   rows 0..nc-1  : A_cbf[row][j]   (CBF part)
    //   row  nc+j     : 1.0             (identity part — only one entry per col)
    int nnz = 0;
    for (int col = 0; col < 3; col++)
    {
        cbf_A_p[col] = nnz;

        // CBF rows for this column (zeros at init, updated each step)
        for (int row = 0; row < nc; row++)
        {
            cbf_A_x[nnz] = 0.0f;
            cbf_A_i[nnz] = row;
            nnz++;
        }

        // Identity row for this column (box constraint on u[col])
        cbf_A_x[nnz] = 1.0f;
        cbf_A_i[nnz] = nc + col;
        nnz++;
    }
    cbf_A_p[3] = nnz;

    cbf_A_mat.m     = nr;
    cbf_A_mat.n     = 3;
    cbf_A_mat.nz    = -1;
    cbf_A_mat.nzmax = nnz;
    cbf_A_mat.x     = cbf_A_x;
    cbf_A_mat.i     = cbf_A_i;
    cbf_A_mat.p     = cbf_A_p;

    // Initial bounds
    for (int i = 0; i < 3;  i++) cbf_q[i] = 0.0f;

    // CBF rows: l = -inf, u = 0 (will be overwritten each step)
    for (int i = 0; i < nc; i++)
    {
        cbf_l[i] = (OSQPFloat)(-OSQP_INFTY);
        cbf_u[i] = 0.0f;
    }

    // Box constraint rows: l = u_min, u = u_max (fixed, never change)
    for (int j = 0; j < 3; j++)
    {
        cbf_l[nc + j] = (OSQPFloat)u_min[j];
        cbf_u[nc + j] = (OSQPFloat)u_max[j];
    }

    // OSQP settings: tuned for real-time use
    osqp_set_default_settings(&cbf_settings);
    cbf_settings.verbose       = 0;
    cbf_settings.warm_starting = 1;
    cbf_settings.eps_abs       = 1e-4;
    cbf_settings.eps_rel       = 1e-4;
    cbf_settings.max_iter      = 1000;

    OSQPInt flag = osqp_setup(&cbf_solver,
                              &cbf_P_mat, cbf_q,
                              &cbf_A_mat, cbf_l, cbf_u,
                              nr, 3, &cbf_settings);

    if (flag != 0)
    {
        fprintf(stderr, "[CBF] OSQP setup failed (code %d)\n", (int)flag);
        return 0;
    }

    return 1;
}


// =============================================================================
//  cbf_circle_obstacles_filter()
//
//  Call after computing your nominal (vx, vy, wz) commands.
//  Modifies commands in-place if any obstacle constraint would be violated.
//
//  Inputs:
//    x, y, theta         — robot world position and yaw
//    vx_cmd, vy_cmd, wz_cmd — nominal body-frame commands (modified in-place)
//    wx, wy, ww          — QP cost weights (how much to penalise each change)
//
//  Returns 1 if commands were modified, 0 otherwise.
// =============================================================================
static inline int cbf_circle_obstacles_filter(
    double  x,      double  y,      double theta,
    double *vx_cmd, double *vy_cmd, double *wz_cmd,
    double  wx,     double  wy,     double  ww)
{
    const double c = cos(theta);
    const double s = sin(theta);

    // =========================================================================
    // STEP 1: Compute CBF constraint rows  A_cbf * u <= b_cbf
    // =========================================================================

    double A_cbf[N_OBS][3];
    double b_cbf[N_OBS];
    int nc = 0;

    for (int i = 0; i < N_OBS; i++)
    {
        if (!obstacles[i].active) continue;

        double dx = x - obstacles[i].cx;
        double dy = y - obstacles[i].cy;
        double h  = dx*dx + dy*dy - obstacles[i].r * obstacles[i].r;

        A_cbf[nc][0] = -(2.0*dx*c + 2.0*dy*s);
        A_cbf[nc][1] = -(2.0*dy*c - 2.0*dx*s);
        A_cbf[nc][2] =  0.0;
        b_cbf[nc]    =  obstacles[i].alpha * h - 0.1;

        nc++;
    }

    if (nc == 0) return 0;

    // =========================================================================
    // STEP 2: Initialize solver on first call
    // =========================================================================

    if (cbf_solver == NULL)
    {
        if (!cbf_init_solver(wx, wy, ww, nc)) return 0;
    }

    // =========================================================================
    // STEP 3: Update OSQP data for this timestep (fast — no re-factorization)
    //
    //  q = -W * u0   (centers the cost at the nominal command)
    //  A values      (updated CBF constraint rows)sssss
    //  u_bound       (updated CBF right-hand side)
    // =========================================================================

    // Update q vector
    cbf_q[0] = (OSQPFloat)(-wx * (*vx_cmd));
    cbf_q[1] = (OSQPFloat)(-wy * (*vy_cmd));
    cbf_q[2] = (OSQPFloat)(-ww * (*wz_cmd));

    // Update CBF rows of A matrix (column-major, only the first nc rows per col)
    int nnz = 0;
    for (int col = 0; col < 3; col++)
    {
        for (int row = 0; row < nc; row++)
            cbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];
        nnz++; // skip the identity entry (stays 1.0, never changes)
    }

    // Update CBF upper bounds only (rows 0..nc-1)
    // Box bounds (rows nc..nc+2) are fixed at u_min/u_max — no update needed
    for (int i = 0; i < nc; i++)
        cbf_u[i] = (OSQPFloat)b_cbf[i];

    // Push updates to solver
    osqp_update_data_vec(cbf_solver, cbf_q, cbf_l, cbf_u);
    osqp_update_data_mat(cbf_solver, NULL, NULL, 0, cbf_A_x, NULL, cbf_A_mat.nzmax);
    // =========================================================================
    // STEP 4: Solve and extract result
    // =========================================================================

    osqp_solve(cbf_solver);

    int modified = 0;

    if (cbf_solver->info->status_val == OSQP_SOLVED ||
        cbf_solver->info->status_val == OSQP_SOLVED_INACCURATE)
    {
        *vx_cmd = (double)cbf_solver->solution->x[0];
        *vy_cmd = (double)cbf_solver->solution->x[1];
        *wz_cmd = (double)cbf_solver->solution->x[2];
        modified = 1;
    }
    else
    {
        fprintf(stderr, "[CBF] OSQP did not converge (%s)\n",
                cbf_solver->info->status);
        // fallback: clamp to input limits even without CBF
        *vx_cmd = fmax(u_min[0], fmin(u_max[0], *vx_cmd));
        *vy_cmd = fmax(u_min[1], fmin(u_max[1], *vy_cmd));
        *wz_cmd = fmax(u_min[2], fmin(u_max[2], *wz_cmd));
    }

    return modified;
}
