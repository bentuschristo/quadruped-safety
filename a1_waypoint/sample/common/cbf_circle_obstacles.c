#pragma once
#include <math.h>
#include <stdio.h>

// -----------------------------------------------------------------------
//  Circular/cylindrical obstacle avoidance via Control Barrier Functions
//
//  CBF:  h(x,y) = (x-cx)^2 + (y-cy)^2 - r^2  >= 0  (stay outside circle)
//
//  Safety condition (1st-order CBF):
//    hdot + alpha*h >= 0
//    hdot = 2*(x-cx)*xdot + 2*(y-cy)*ydot
//
//  Body-to-world velocity mapping:
//    xdot =  cos(theta)*vx - sin(theta)*vy
//    ydot =  sin(theta)*vx + cos(theta)*vy
//
//  Rewritten as a linear inequality on [vx, vy, wz] (A*u <= b form):
//    a0*vx + a1*vy + 0*wz <= alpha*h
//  where:
//    a0 = -2*(x-cx)*c - 2*(y-cy)*s      (c = cos theta, s = sin theta)
//    a1 =  2*(x-cx)*s - 2*(y-cy)*c
//
//  With two obstacles this gives a QP with 2 linear inequality constraints.
//  Solved analytically via active-set enumeration (no external solver).
// -----------------------------------------------------------------------

// Obstacle parameters
typedef struct {
    double cx;     // center x  (world frame)
    double cy;     // center y  (world frame)
    double r;      // radius (safety margin already included)
    double alpha;  // CBF gain — larger = more aggressive near boundary
    int    active; // 1 = enabled, 0 = disabled
} CircleCBF;

// Velocity command triple
typedef struct {
    double vx, vy, wz;
} QP3;

// ---- Obstacle definitions (edit as needed) ----
static CircleCBF CBF_CIRC1 = { -2.0,  1, 0.35, 0.5, 1 };
static CircleCBF CBF_CIRC2 = { -4.0, -1, 0.35, 0.5, 1 };
// -----------------------------------------------


// ---- Internal QP helpers ----------------------

static double Au_eval(int i, const double A[2][3], const QP3 *u)
{
    return A[i][0]*u->vx + A[i][1]*u->vy + A[i][2]*u->wz;
}

static int qp_feasible(const double A[2][3], const double b[2], const QP3 *u)
{
    return (Au_eval(0, A, u) <= b[0] + 1e-12) &&
           (Au_eval(1, A, u) <= b[1] + 1e-12);
}

static double obj_cost_diagW(const QP3 *u, const QP3 *u0,
                              double wx, double wy, double ww)
{
    const double dx = u->vx - u0->vx;
    const double dy = u->vy - u0->vy;
    const double dw = u->wz - u0->wz;
    return 0.5 * (wx*dx*dx + wy*dy*dy + ww*dw*dw);
}

// Solve   min  0.5*(u-u0)^T W (u-u0)   s.t.  A u <= b
// W = diag(wx, wy, ww),  2 constraints,  3 unknowns.
// Active-set enumeration: 0, 1, or 2 active constraints.
static QP3 cbf_qp_solve_3var_2ineq(
    const QP3   *u0,
    const double A[2][3],
    const double b[2],
    double wx, double wy, double ww)
{
    const double inv_wx = 1.0 / wx;
    const double inv_wy = 1.0 / wy;
    const double inv_ww = 1.0 / ww;

    QP3 best = *u0;
    double best_cost = 1e300;

    // ---- 0 active constraints (nominal is feasible) ----
    if (qp_feasible(A, b, u0)) return *u0;

    // ---- 1 active constraint (project onto constraint i) ----
    for (int i = 0; i < 2; i++)
    {
        const double a0 = A[i][0], a1 = A[i][1], a2 = A[i][2];

        // M = a W^{-1} a^T  (scalar)
        const double denom = a0*a0*inv_wx + a1*a1*inv_wy + a2*a2*inv_ww;
        if (denom < 1e-12) continue;

        // KKT multiplier  lambda = (A u0 - b) / M
        const double lhs0 = Au_eval(i, A, u0);
        const double lambda = (lhs0 - b[i]) / denom;
        if (lambda < 0.0) continue; // constraint not actually active

        QP3 u;
        u.vx = u0->vx - inv_wx * a0 * lambda;
        u.vy = u0->vy - inv_wy * a1 * lambda;
        u.wz = u0->wz - inv_ww * a2 * lambda;

        if (!qp_feasible(A, b, &u)) continue;

        const double J = obj_cost_diagW(&u, u0, wx, wy, ww);
        if (J < best_cost) { best_cost = J; best = u; }
    }

    // ---- 2 active constraints (both equalities) ----
    {
        // M = A W^{-1} A^T  (2x2 symmetric)
        double M00 = 0, M01 = 0, M11 = 0;
        for (int k = 0; k < 3; k++)
        {
            const double inv_wk = (k == 0) ? inv_wx : (k == 1) ? inv_wy : inv_ww;
            M00 += A[0][k] * inv_wk * A[0][k];
            M01 += A[0][k] * inv_wk * A[1][k];
            M11 += A[1][k] * inv_wk * A[1][k];
        }

        const double det = M00*M11 - M01*M01;
        if (fabs(det) >= 1e-12)
        {
            const double rhs0 = Au_eval(0, A, u0) - b[0];
            const double rhs1 = Au_eval(1, A, u0) - b[1];

            // lambda = M^{-1} * rhs
            const double lam0 = ( M11*rhs0 - M01*rhs1) / det;
            const double lam1 = (-M01*rhs0 + M00*rhs1) / det;

            if (lam0 >= 0.0 && lam1 >= 0.0)
            {
                QP3 u = *u0;
                u.vx -= inv_wx * (A[0][0]*lam0 + A[1][0]*lam1);
                u.vy -= inv_wy * (A[0][1]*lam0 + A[1][1]*lam1);
                u.wz -= inv_ww * (A[0][2]*lam0 + A[1][2]*lam1);

                if (qp_feasible(A, b, &u))
                {
                    const double J = obj_cost_diagW(&u, u0, wx, wy, ww);
                    if (J < best_cost) { best_cost = J; best = u; }
                }
            }
        }
    }

    // Last resort: no feasible candidate found, return nominal unchanged
    if (best_cost >= 1e250) return *u0;

    return best;
}


// -----------------------------------------------------------------------
//  Public interface
//
//  cbf_circle_obstacles_filter()
//
//  Inputs:
//    x, y, theta   : world-frame robot position and yaw
//    vx_cmd etc.   : nominal body-frame velocity commands (modified in-place)
//    wx, wy, ww    : diagonal QP cost weights (how much to penalise each deviation)
//
//  Returns 1 if any command was modified, 0 otherwise.
// -----------------------------------------------------------------------
static inline int cbf_circle_obstacles_filter(
    double  x,
    double  y,
    double  theta,
    double *vx_cmd,
    double *vy_cmd,
    double *wz_cmd,
    double  wx,
    double  wy,
    double  ww)
{
    const double c = cos(theta);
    const double s = sin(theta);

    // Build constraint matrix A (2x3) and rhs b (2)
    // Row i corresponds to obstacle i.
    double A[2][3] = {{0}};
    double b[2]    = {0};
    int    n_active = 0; // count actually-binding obstacles

    const CircleCBF *obs[2] = { &CBF_CIRC1, &CBF_CIRC2 };

    for (int i = 0; i < 2; i++)
    {
        if (!obs[i]->active) continue;

        const double dx = x - obs[i]->cx;
        const double dy = y - obs[i]->cy;
        const double h  = dx*dx + dy*dy - obs[i]->r * obs[i]->r;

        // hdot constraint: a*u >= -alpha*h  =>  -a*u <= alpha*h
        // a[vx] = 2*dx*c + 2*dy*s
        // a[vy] = -2*dx*s + 2*dy*c      (note sign: ydot = s*vx + c*vy)
        // a[wz] = 0
        A[i][0] = -(2.0*dx*c + 2.0*dy*s);
        A[i][1] = -(2.0*dy*c - 2.0*dx*s);
        A[i][2] =  0.0;
        b[i]    =  obs[i]->alpha * h;

        n_active++;
    }

    if (n_active == 0) return 0;

    QP3 u0 = { *vx_cmd, *vy_cmd, *wz_cmd };

    // If already feasible, nothing to do
    if (qp_feasible(A, b, &u0)) return 0;

    QP3 u_safe = cbf_qp_solve_3var_2ineq(&u0, (const double (*)[3])A, b, wx, wy, ww);

    const int modified = (u_safe.vx != u0.vx ||
                          u_safe.vy != u0.vy ||
                          u_safe.wz != u0.wz);

    *vx_cmd = u_safe.vx;
    *vy_cmd = u_safe.vy;
    *wz_cmd = u_safe.wz;

    return modified;
}
