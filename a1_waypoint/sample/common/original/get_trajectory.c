extern int g_traj_type;
extern double a;

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
