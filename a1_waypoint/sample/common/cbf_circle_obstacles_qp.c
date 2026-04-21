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
    { -2.0,  1.0, 0.35, 1, 1 },
    { -4.0, -1.0, 0.35, 1, 1 },
    { -2.0, -1.0, 0.35, 1, 1 },
    { -4.0,  1.0, 0.35, 1, 1 },
};

// static CircleCBF obstacles[N_OBS] = {
//     { -2.0,  1.0, 0.35, 0.5, 1 },   // obstacle 1
//     { -4.0, -1.0, 0.35, 0.5, 1 },   // obstacle 2
// };
// -----------------------------------------------------------------------------


// ---- Static solver state (persists across calls) ----------------------------
static OSQPSolver   *cbf_solver   = NULL;
static OSQPSettings  cbf_settings;

// Fixed-size storage for P and A matrices (upper bound: N_OBS constraints)
// These arrays are reused every call — only their values change, not structure.
static OSQPFloat cbf_P_x[3];          // P diagonal: wx, wy, ww
static OSQPInt   cbf_P_i[3] = {0,1,2};
static OSQPInt   cbf_P_p[4] = {0,1,2,3};
static OSQPCscMatrix cbf_P_mat;

static OSQPFloat cbf_A_x[N_OBS * 3]; // constraint matrix values (col-major)
static OSQPInt   cbf_A_i[N_OBS * 3]; // row indices
static OSQPInt   cbf_A_p[4];          // column pointers
static OSQPCscMatrix cbf_A_mat;

static OSQPFloat cbf_q[3];            // linear cost: -W * u0
static OSQPFloat cbf_l[N_OBS];        // lower bounds (-infinity)
static OSQPFloat cbf_u[N_OBS];        // upper bounds (b_cbf)
// -----------------------------------------------------------------------------


// ---- Initialize solver (called once on first use) ---------------------------
static int cbf_init_solver(double wx, double wy, double ww, int nc)
{
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

    // A matrix: nc x 3, column-major CSC, filled with zeros for init
    int nnz = 0;
    for (int col = 0; col < 3; col++)
    {
        cbf_A_p[col] = nnz;
        for (int row = 0; row < nc; row++)
        {
            cbf_A_x[nnz] = 0.0f;
            cbf_A_i[nnz] = row;
            nnz++;
        }
    }
    cbf_A_p[3] = nnz;

    cbf_A_mat.m     = nc;
    cbf_A_mat.n     = 3;
    cbf_A_mat.nz    = -1;
    cbf_A_mat.nzmax = nnz;
    cbf_A_mat.x     = cbf_A_x;
    cbf_A_mat.i     = cbf_A_i;
    cbf_A_mat.p     = cbf_A_p;

    // Initial q and bounds (will be overwritten each call)
    for (int i = 0; i < 3;  i++) cbf_q[i] = 0.0f;
    for (int i = 0; i < nc; i++) { cbf_l[i] = (OSQPFloat)(-OSQP_INFTY); cbf_u[i] = 0.0f; }

    // OSQP settings: tuned for real-time use
    osqp_set_default_settings(&cbf_settings);
    cbf_settings.verbose       = 0;
    cbf_settings.warm_starting = 1;   // reuse previous solution each step
    cbf_settings.eps_abs       = 1e-4;
    cbf_settings.eps_rel       = 1e-4;
    cbf_settings.max_iter      = 200; // low cap keeps latency bounded

    OSQPInt flag = osqp_setup(&cbf_solver,
                              &cbf_P_mat, cbf_q,
                              &cbf_A_mat, cbf_l, cbf_u,
                              nc, 3, &cbf_settings);

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
    //  A values      (updated CBF constraint rows)
    //  u_bound       (updated CBF right-hand side)
    // =========================================================================

    // Update q vector
    cbf_q[0] = (OSQPFloat)(-wx * (*vx_cmd));
    cbf_q[1] = (OSQPFloat)(-wy * (*vy_cmd));
    cbf_q[2] = (OSQPFloat)(-ww * (*wz_cmd));

    // Update A matrix values (column-major order)
    int nnz = 0;
    for (int col = 0; col < 3; col++)
        for (int row = 0; row < nc; row++)
            cbf_A_x[nnz++] = (OSQPFloat)A_cbf[row][col];

    // Update upper bounds
    for (int i = 0; i < nc; i++)
        cbf_u[i] = (OSQPFloat)b_cbf[i];

    // Push updates to solver (only vectors, sparsity pattern stays fixed)
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
    }

    return modified;
}
