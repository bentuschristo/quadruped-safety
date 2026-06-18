extern int g_traj_type;
extern double a;

/* --------------------------------------------------------------------------
 * Editable waypoint path for g_traj_type == 4
 *
 * Add, remove, or modify rows below. Each row is {x, y}.
 * The path uses a C1-continuous cubic Hermite spline parameterized by
 * cumulative chord length. End-point tangents are one-sided; interior
 * tangents are centered.
 * -------------------------------------------------------------------------- */
static const double waypoint_path[][2] = {
    { 0.0, 0.0}, // start: should match the robot's initial position
    { 0.0, 1.0},
    { 0.2, 2.0},
    { 0.0, 3.0},
    {-0.5, 4.0},
    {-0.7, 5.0},
    {-0.3, 6.0},
    { 0.0, 7.0},
    { 0.0, 8.0} // Endpoint
};

#define WAYPOINT_COUNT ((int)(sizeof(waypoint_path) / sizeof(waypoint_path[0])))


static double wrap_to_pi_local(double angle)
{
    while (angle > pi)  angle -= 2.0*pi;
    while (angle < -pi) angle += 2.0*pi;
    return angle;
}

/* Return the equivalent of angle that is closest to reference. */
static double unwrap_near_local(double angle, double reference)
{
    return reference + wrap_to_pi_local(angle - reference);
}

static double clamp01(double v)
{
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

static void eval_waypoint_spline(double t, double duration,
                                 double *x, double *y,
                                 double *xdot, double *ydot,
                                 double *xddot, double *yddot)
{
    double s[WAYPOINT_COUNT];
    double tx[WAYPOINT_COUNT];
    double ty[WAYPOINT_COUNT];

    s[0] = 0.0;
    for (int i = 1; i < WAYPOINT_COUNT; ++i) {
        double dx = waypoint_path[i][0] - waypoint_path[i-1][0];
        double dy = waypoint_path[i][1] - waypoint_path[i-1][1];
        double ds = sqrt(dx*dx + dy*dy);
        if (ds < 1e-9) ds = 1e-9;
        s[i] = s[i-1] + ds;
    }

    for (int i = 0; i < WAYPOINT_COUNT; ++i) {
        if (i == 0) {
            double ds = s[1] - s[0];
            tx[i] = (waypoint_path[1][0] - waypoint_path[0][0]) / ds;
            ty[i] = (waypoint_path[1][1] - waypoint_path[0][1]) / ds;
        } else if (i == WAYPOINT_COUNT - 1) {
            double ds = s[i] - s[i-1];
            tx[i] = (waypoint_path[i][0] - waypoint_path[i-1][0]) / ds;
            ty[i] = (waypoint_path[i][1] - waypoint_path[i-1][1]) / ds;
        } else {
            double ds = s[i+1] - s[i-1];
            tx[i] = (waypoint_path[i+1][0] - waypoint_path[i-1][0]) / ds;
            ty[i] = (waypoint_path[i+1][1] - waypoint_path[i-1][1]) / ds;
        }
    }

    double rho = clamp01(t / duration);
    double sigma = 10.0*rho*rho*rho - 15.0*rho*rho*rho*rho + 6.0*rho*rho*rho*rho*rho;
    double sigma_dot = (30.0*rho*rho - 60.0*rho*rho*rho + 30.0*rho*rho*rho*rho) / duration;
    double sigma_ddot = (60.0*rho - 180.0*rho*rho + 120.0*rho*rho*rho) / (duration*duration);

    double s_total = s[WAYPOINT_COUNT - 1];
    double s_now = sigma * s_total;
    double s_dot = sigma_dot * s_total;
    double s_ddot = sigma_ddot * s_total;

    int seg = WAYPOINT_COUNT - 2;
    for (int i = 0; i < WAYPOINT_COUNT - 1; ++i) {
        if (s_now <= s[i+1]) { seg = i; break; }
    }

    double ds = s[seg+1] - s[seg];
    double u = (s_now - s[seg]) / ds;
    u = clamp01(u);

    double h00 =  2*u*u*u - 3*u*u + 1;
    double h10 =      u*u*u - 2*u*u + u;
    double h01 = -2*u*u*u + 3*u*u;
    double h11 =      u*u*u -   u*u;

    double dh00 =  6*u*u - 6*u;
    double dh10 =  3*u*u - 4*u + 1;
    double dh01 = -6*u*u + 6*u;
    double dh11 =  3*u*u - 2*u;

    double d2h00 = 12*u - 6;
    double d2h10 =  6*u - 4;
    double d2h01 = -12*u + 6;
    double d2h11 =  6*u - 2;

    double x0 = waypoint_path[seg][0], x1 = waypoint_path[seg+1][0];
    double y0 = waypoint_path[seg][1], y1 = waypoint_path[seg+1][1];

    *x = h00*x0 + h10*ds*tx[seg] + h01*x1 + h11*ds*tx[seg+1];
    *y = h00*y0 + h10*ds*ty[seg] + h01*y1 + h11*ds*ty[seg+1];

    double dx_du = dh00*x0 + dh10*ds*tx[seg] + dh01*x1 + dh11*ds*tx[seg+1];
    double dy_du = dh00*y0 + dh10*ds*ty[seg] + dh01*y1 + dh11*ds*ty[seg+1];
    double d2x_du2 = d2h00*x0 + d2h10*ds*tx[seg] + d2h01*x1 + d2h11*ds*tx[seg+1];
    double d2y_du2 = d2h00*y0 + d2h10*ds*ty[seg] + d2h01*y1 + d2h11*ds*ty[seg+1];

    double du_dt = s_dot / ds;
    double d2u_dt2 = s_ddot / ds;

    *xdot = dx_du * du_dt;
    *ydot = dy_du * du_dt;
    *xddot = d2x_du2 * du_dt * du_dt + dx_du * d2u_dt2;
    *yddot = d2y_du2 * du_dt * du_dt + dy_du * d2u_dt2;
}

void get_trajectory(double ttime,double ts,double tend,
                    double *xref,double *yref,
                    double *xdotref,double *ydotref,
                    double *thetaref,double *thetadotref,
                    double a,double x_center,double y_center)
{
  if (g_traj_type == 1) {
      // lemniscate code
      double tau = 2*pi*(ttime-ts)/tend;
      double  den = 1+sin(tau)*sin(tau);
      double den2 = den*den;
      double b = 2*pi/tend;
      *xref = x_center + a*cos(tau)/den;
      *yref = y_center + a*cos(tau)*sin(tau)/den;
      *xdotref = (a*b*sin(tau)*(sin(tau)*sin(tau) - 3))/den2;
      *ydotref = -(a*b*(3*sin(tau)*sin(tau) - 1))/den2;

      double xdot = *xdotref;
      double ydot = *ydotref;
      double epsilon = 1e-6;
      double xddot = a*(b*b)*cos(tau)*1.0/pow(pow(cos(tau),2.0)-2.0,3.0)*(pow(cos(tau),2.0)*1.0E+1+pow(cos(tau),4.0)-8.0);
      double yddot = a*(b*b)*cos(tau)*sin(tau)*1.0/pow(pow(sin(tau),2.0)+1.0,3.0)*(pow(cos(tau),2.0)*3.0+2.0)*-2.0;

      double theta = atan2(ydot,xdot);
      double sec_theta = 1/cos(theta);
      double num3 = xdot*yddot - ydot*xddot;
      double den3 = (xdot*xdot+epsilon);
      double thetadot;

      //constant reference
        if (flag_constant)
          thetadot = 0;
        else
          thetadot = (1/(sec_theta*sec_theta))*(num3/den3);

        //changing reference
      *thetadotref = thetadot;
      theta_ref0 += t_stance*thetadot;
      *thetaref = theta_ref0;
  }
  else if (g_traj_type == 2) {
      // circle code
      // a = radius of circle
      // tend = period (time to complete one lap)

      double omega = 2.0 * pi / tend;
      double tau   = omega * (ttime - ts);

      /* Position (circle) */
      *xref = x_center + a * cos(tau);
      *yref = y_center + a * sin(tau);

      /* Velocity */
      *xdotref = -a * omega * sin(tau);
      *ydotref =  a * omega * cos(tau);

      /* Acceleration */
      double xddot = -a * omega * omega * cos(tau);
      double yddot = -a * omega * omega * sin(tau);

      /* Heading = tangent direction */
      double xdot = *xdotref;
      double ydot = *ydotref;

      double v2 = xdot*xdot + ydot*ydot + 1e-9;

      double thetadot = (xdot*yddot - ydot*xddot) / v2;

      if (flag_constant)
          thetadot = 0.0;

      *thetadotref = thetadot;

      /* Integrate yaw like original code */
      theta_ref0 += t_stance * thetadot;
      *thetaref = theta_ref0;
  }
  else if (g_traj_type == 3) {
      // straight-line code
      // Straight line toward -Y direction.
      // Reuse:
      //   a    = speed (m/s)
      //   tend = duration (s)
      //   (x_center, y_center) = start position at t=ts

      double t = ttime - ts;
      if (t < 0.0) t = 0.0;
      if (t > tend) t = tend;

      // Position: x constant, y decreases linearly
      *xref = x_center;
      *yref = y_center - a * t;

      // Velocity: constant along -Y
      *xdotref = 0.0;
      *ydotref = -a;

      // Straight line => zero curvature => thetadot = 0
      double thetadot = 0.0;

      if (flag_constant) {
          thetadot = 0.0;
      }

      *thetadotref = thetadot;

      // Keep your original yaw integration behavior
      theta_ref0 += t_stance * thetadot;
      *thetaref = theta_ref0;
  }

  else if (g_traj_type == 4) {
      // Editable waypoint-based forest trajectory.
      // tend = total traversal duration.
      double t = ttime - ts;
      if (t < 0.0) t = 0.0;
      if (t > tend) t = tend;

      double xddot = 0.0, yddot = 0.0;
      eval_waypoint_spline(t, tend,
                           xref, yref,
                           xdotref, ydotref,
                           &xddot, &yddot);

      /*
       * Heading option:
       *   flag_constant != 0 : hold the existing constant heading theta_ref0.
       *   flag_constant == 0 : align heading with the path tangent so that
       *                        the reference motion is primarily forward (v_x).
       */
      if (flag_constant) {
          *thetadotref = 0.0;
          *thetaref = theta_ref0;
      } else {
          double v2 = (*xdotref)*(*xdotref) + (*ydotref)*(*ydotref);

          if (v2 > 1e-9) {
              double theta_tangent = atan2(*ydotref, *xdotref);

              /* Prevent artificial 2*pi jumps when the tangent crosses +/-pi. */
              theta_tangent = unwrap_near_local(theta_tangent, theta_ref0);

              *thetaref = theta_tangent;
              *thetadotref = ((*xdotref)*yddot - (*ydotref)*xddot) / v2;
              theta_ref0 = theta_tangent;
          } else {
              /* At the start/end of the fifth-order timing law, speed is zero.
               * Hold the most recent tangent heading until the path starts moving.
               */
              *thetaref = theta_ref0;
              *thetadotref = 0.0;
          }
      }
  }
  else {
      // fallback default: lemniscate
      double tau = 2*pi*(ttime-ts)/tend;
      double  den = 1+sin(tau)*sin(tau);
      double den2 = den*den;
      double b = 2*pi/tend;
      *xref = x_center + a*cos(tau)/den;
      *yref = y_center + a*cos(tau)*sin(tau)/den;
      *xdotref = (a*b*sin(tau)*(sin(tau)*sin(tau) - 3))/den2;
      *ydotref = -(a*b*(3*sin(tau)*sin(tau) - 1))/den2;

      double xdot = *xdotref;
      double ydot = *ydotref;
      double epsilon = 1e-6;
      double xddot = a*(b*b)*cos(tau)*1.0/pow(pow(cos(tau),2.0)-2.0,3.0)*(pow(cos(tau),2.0)*1.0E+1+pow(cos(tau),4.0)-8.0);
      double yddot = a*(b*b)*cos(tau)*sin(tau)*1.0/pow(pow(sin(tau),2.0)+1.0,3.0)*(pow(cos(tau),2.0)*3.0+2.0)*-2.0;

      double theta = atan2(ydot,xdot);
      double sec_theta = 1/cos(theta);
      double num3 = xdot*yddot - ydot*xddot;
      double den3 = (xdot*xdot+epsilon);
      double thetadot;

      //constant reference
        if (flag_constant)
          thetadot = 0;
        else
          thetadot = (1/(sec_theta*sec_theta))*(num3/den3);

        //changing reference
      *thetadotref = thetadot;
      theta_ref0 += t_stance*thetadot;
      *thetaref = theta_ref0;
  }


  
  
}
