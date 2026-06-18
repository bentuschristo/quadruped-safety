double setpt[nact]={0};
int step_no = 0;
int prev_step = 0;
double theta_ref0 = 0;

double Kp_ik = 2; //original = 2.5

int flag_constant = 0; //1 for constant yaw and 0 for varying yaw

int cbf_on = 1;

double px = 0.0;
double py = 0;
#define t_curve 40
double tend = t_curve;
double a = 3;

#include <sys/time.h>

#include "utility_matrix.c"
#include "set_model_state.c"
#include "poly_trajectory.c"
#include "get_Jac.c"
#include "get_kinematics2.c"
#include "inverse_kinematics_3D_analytic.c"
#include "control_params.h"
#include "get_trq_stance.c"
#include "get_trq_swing.c"
#include "set_limit.c"
#include "leg_state_machine.c"
#include "get_stance_forces.c"
#include "get_trajectory.c"

// #include "cbf_circle_obstacles_qp.c"
// #include "cbf_circle_obstacles_qp_lim.c"
// #include "cbf_circle_obstacles_cf.c"

// ==== Baselines ====
// #include "bcbf_circle_obstacles.c"
// #include "bcbf_circle_obstacles_v2.c"

// #include "bcbf_lse_circle_obstacles.c"
// #include "bcbf_lse_circle_obstacles_v2.c"

// #include "oi_qp_circle_obstacles.c"
// #include "oi_qp_circle_obstacles_v2.c"
// #include "oi_cf_circle_obstacles.c"
// #include "oi_cf_circle_obstacles_v2.c"

// #include "blended_circle_obstacles.c"
// #include "blended_circle_obstacles_v2.c"
// ==================

// #include "bcbf_lse_control_dynamics_closedform.c"
// #include "bcbf_lse_control_dynamics_closedform_ybackup.c"

// #include "bcbf_cf_zero.c"
// #include "bcbf_cf_y.c"

// ===Proposed Method ===
// #include "bcbf_cf_zero_v2.c"
// #include "bcbf_cf_y_v2.c"
// ==================

// ===== Control Dynamics Version =====
#include "bcbf_circle_obstacles_v2_cd.c"
// #include "bcbf_lse_circle_obstacles_v2_cd.c"
// #include "oi_qp_circle_obstacles_v2_cd.c"
// #include "oi_cf_circle_obstacles_v2_cd.c"
// #include "blended_circle_obstacles_v2_cd.c"
// =================


// ************ Change this as per the problem being solved *********//
#define DATA_PTS 3 //Set this based on columns in data file
float data[5000][DATA_PTS]; //Structure that will store the data.
int STEPS;
int COUNTER=0;
// **************************************************************** //

void read_data()
{
    int i, j;

	//fid = fopen(DATA_FILE,"r");
  fidC = fopen(pathC,"r");

	/* Read file */
    i = 0;
    while(!feof(fidC))
    {
        for (j=0;j<DATA_PTS;j++)
            fscanf(fidC, "%f", &data[i][j]);
        i = i+1;
    }
    fclose(fidC);
    STEPS = i-1;

	/* For Display only */
    // for (i=0;i<STEPS;i++)
    // {
    //     for (j=0;j<DATA_PTS;j++)
    //         printf("%f ",data[i][j]);
    //     printf("\n");
    // }
}



int fsm[4] = {0};

void init_controller(int motiontime, double quat_act[4], double omega_act[3],
                     double q_act[nact], double u_act[nact])
{
    int i,j;

    sdinit(); sdprinterr(stderr);


    for (j=0;j<4;j++)
      fsm[j] = fsm_init;

    fidmat = fopen(pathmat,"w");
    fprintf(fidmat,"all_data = [ \n");
    flag_openmat = 1;

  fidC = fopen(pathC,"w");
  flag_C = 1;

  double ang1,ang2, theta;
  double dc[3][3];
  sdquat2dc(quat_act[1],quat_act[2],quat_act[3],quat_act[0],dc); sdprinterr(stderr);
  sddc2ang(dc,&ang1,&ang2,&theta); sdprinterr(stderr);
  theta_ref0 = theta; //initialize to initial position

  // Keep the control-dynamics filter integration step synchronized with the actual controller call period.
  // This matters for filters that define static double cbf_cd_dt inside the included CBF file.
  // cbf_cd_dt = dtime;

}

void my_controller(int motiontime,double quat_act[4], double omega_act[3],
                   double q_act[nact], double u_act[nact], double tau[nact],
                  double Kp[nact],double Kd[nact],
                  double q_ref[nact], double u_ref[nact],double tauEst[nact],
                  uint16_t footForce[4]) //, double tracking_error[2])
{
  int i,j;
  double tau_r[3]={0}, kp_r[3]={0}, kd_r[3] = {0};
  double F[4][3]={0};

  double ttime = motiontime*dtime;
  double ts = 5; // > t_free + t_trans + t_stance

  // -------------------------------------------------------------------------
  // Pose is computed EVERY low-level call. The safety filter should use the
  // newest pose, not only the pose at the last gait-step update.
  // -------------------------------------------------------------------------
  double tmp_dc[3][3]={0};
  double ang1,ang2,theta;
  sdquat2dc(quat_act[1],quat_act[2],quat_act[3],quat_act[0],tmp_dc); sdprinterr(stderr);
  sddc2ang(tmp_dc,&ang1,&ang2,&theta); sdprinterr(stderr);

  double c = cos(theta);
  double s = sin(theta);

  // Robot body/world coordinates.
  double x_c = bodyWorldPos[0];
  double y_c = bodyWorldPos[1];

  // Position of point P used by the tracking controller.
  double x_p = x_c + c*px - s*py;
  double y_p = y_c + s*px + c*py;

  // -------------------------------------------------------------------------
  // Nominal high-level command memory.
  // The trajectory/IK command is updated only when step_no increments, but the
  // safety filter below is called EVERY my_controller call using this held
  // nominal command and the latest pose.
  // -------------------------------------------------------------------------
  static double vx_nom = 0.0;
  static double vy_nom = 0.0;
  static double omega_nom = 0.0;

  static double x_ref_log = 0.0;
  static double y_ref_log = 0.0;
  static double theta_ref_log = 0.0;
  static double err_x_log = 0.0;
  static double err_y_log = 0.0;
  static double err_theta_log = 0.0;

  // Low-level/gait computation still runs at the usual controller rate.
  get_stance_forces(quat_act,omega_act,q_act,u_act,F);

  for (j=0;j<4;j++)
  {
    double qr_ref[3]={0}, ur_ref[3]={0};
    double qr_act[3]={0}, ur_act[3]={0};
    for (i=0;i<3;i++)
    {
      qr_act[i]=q_act[index_leg_act[j][i]];
      ur_act[i]=u_act[index_leg_act[j][i]];
    }

    double force[3] = {0};
    for (i=0;i<3;i++)
      force[i] = F[j][i];

    fsm[j] = leg_state_machine(j,fsm[j],motiontime,quat_act,omega_act,qr_act,ur_act,qr_ref,ur_ref,
                                 tau_r,kp_r,kd_r,footForce[j],force);

    for (i=0;i<3;i++)
    {
      q_ref[index_leg_act[j][i]] = qr_ref[i];
      u_ref[index_leg_act[j][i]] = ur_ref[i];
      tau[index_leg_act[j][i]] = tau_r[i];
      Kp[index_leg_act[j][i]] = kp_r[i];
      Kd[index_leg_act[j][i]] = kd_r[i];
    }
  }

  // -------------------------------------------------------------------------
  // Update the NOMINAL trajectory command only once per gait step.
  // This keeps your original waypoint/trajectory update structure.
  // -------------------------------------------------------------------------
  if (ttime >= ts && ttime <= ts+tend)
  {
    if (prev_step < step_no)
    {
      prev_step = step_no;

      double x_ref, y_ref;
      double xdot_ref, ydot_ref;
      double theta_ref, thetadot_ref;

      double x_center = -a;
      double y_center = 0;
      get_trajectory(ttime,ts,tend,&x_ref,&y_ref,
                      &xdot_ref,&ydot_ref,
                      &theta_ref,&thetadot_ref,
                      a,x_center,y_center);

      double dc[3][3]={0}, dcT[3][3]={0}, dc_ref[3][3]={0}, dc_diff[3][3]={0};
      double euler_ref[3]={0,0,theta_ref};
      double euler[3] ={0,0,theta};
      sdang2dc(euler[0],euler[1],euler[2],dc); sdprinterr(stderr);
      sdang2dc(euler_ref[0],euler_ref[1],euler_ref[2],dc_ref); sdprinterr(stderr);
      ram_transpose(&dc[0][0],3,3,&dcT[0][0]);
      ram_multMatMat(&dc_ref[0][0],3,3,&dcT[0][0],3,3,&dc_diff[0][0]);
      double error_euler[3]={0};
      sddc2ang(dc_diff,&error_euler[0],&error_euler[1],&error_euler[2]); sdprinterr(stderr);
      double error_theta = error_euler[2];

      double xdot_p, ydot_p, thetadot_p;
      xdot_p = xdot_ref + Kp_ik*(x_ref - x_p);
      ydot_p = ydot_ref + Kp_ik*(y_ref - y_p);
      thetadot_p = thetadot_ref + Kp_ik*error_theta;

      // Xdot = B*U where U = [vx; vy; omega]. Here this is inv(B).
      double a11 = c; double a12 = s; double a13 = py;
      double a21 = -s; double a22 =  c; double a23 = -px;
      double a31 = 0;  double a32 = 0; double a33 = 1;

      vx_nom    = a11*xdot_p + a12*ydot_p + a13*thetadot_p;
      vy_nom    = a21*xdot_p + a22*ydot_p + a23*thetadot_p;
      omega_nom = a31*xdot_p + a32*ydot_p + a33*thetadot_p;

      x_ref_log = x_ref;
      y_ref_log = y_ref;
      theta_ref_log = theta_ref;
      err_x_log = x_ref - x_p;
      err_y_log = y_ref - y_p;
      err_theta_log = error_theta;
    }

    // -----------------------------------------------------------------------
    // Run safety filter EVERY low-level call, not only every gait step.
    // The nominal command is held between gait-step updates, but the safety
    // correction uses the latest x_c, y_c, theta each call.
    // -----------------------------------------------------------------------
    vx = vx_nom;
    vy = vy_nom;
    omega = omega_nom;

    if (cbf_on == 1){
      cbf_circle_obstacles_filter(x_c, y_c, theta,
                                  &vx, &vy, &omega,
                                  1.0, 1.0, 1.0);
    }

    // Avoid writing huge logs unless you want full 1 kHz logs.
    // Change 10 to 1 for every-call logging, or 100 for lighter logs.
    if (motiontime % 10 == 0) {
      fprintf(fidC,"%d %f %f %f %f %f %f %f %f %f %f %d %.9e\n",
              step_no,
              vx, vy, omega,
              cbf_log_u_nom[0], cbf_log_u_nom[1], cbf_log_u_nom[2],
              cbf_log_u_backup[0], cbf_log_u_backup[1], cbf_log_u_backup[2],
              cbf_log_mu, cbf_log_status, cbf_log_qp_time_sec);

      fprintf(fidmat,"%d %f %f %f %f %f %f %f %f %f %f %f %f; \n",
              step_no, x_ref_log, y_ref_log, theta_ref_log,
              x_p, y_p, theta_ref_log - err_theta_log,
              err_x_log, err_y_log, err_theta_log,
              vx, vy, omega);
    }
  }
  else if (ttime > ts+tend)
  {
    if (flag_openmat==1)
    {
      fprintf(fidmat,"]; ");
      fclose(fidmat);
      flag_openmat = 0;
    }

    if (flag_C==1)
    {
      fclose(fidC);
      flag_C = 0;
    }

    vx = 0;
    omega = 0;
    vy = 0;
  }
}
