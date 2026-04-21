#pragma once
#include <math.h>
#include <stdio.h>

// Enforce x <= x_max using a 1st-order CBF on h = x_max - x
// Modifies vx_cmd in-place (optionally using vy_cmd for coupling).
//
// Inputs:
//  x, theta: world position x and yaw
//  vy_cmd: body-frame lateral velocity command (used for coupling term)
// Params:
//  x_max: boundary
//  alpha: CBF gain (bigger = more aggressive near boundary)
// Returns:
//  1 if it modified vx_cmd, 0 otherwise
static inline int cbf_geofence_xmax_filter(
    double x,
    double theta,
    double vy_cmd,
    double x_max,
    double alpha,
    double *vx_cmd
){
    const double c = cos(theta);
    const double s = sin(theta);

    // h = x_max - x
    const double h = x_max - x;

    // Desired inequality: xdot <= alpha*h
    // xdot ≈ c*vx - s*vy
    // => c*vx <= alpha*h + s*vy
    // => vx <= (alpha*h + s*vy)/c   if c>0
    // => vx >= (alpha*h + s*vy)/c   if c<0 (inequality flips)
    const double eps = 1e-6;

    // If facing nearly sideways, vx has almost no effect on world-x.
    // In that case, don't do anything (or you can clamp vy instead).
    if (fabs(c) < eps) return 0;

    const double rhs = (alpha * h + s * vy_cmd) / c;

    int modified = 0;

    if (c > 0.0) {
        if (*vx_cmd > rhs) { *vx_cmd = rhs; modified = 1; }
    } else {
        // c < 0: inequality flips because dividing by negative
        if (*vx_cmd < rhs) { *vx_cmd = rhs; modified = 1; }
    }

    return modified;
}
