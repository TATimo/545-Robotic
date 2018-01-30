/*============================================================================
==============================================================================
                      
                              min_jerk_task.cpp
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"

// defines

// local variables
static double      start_time = 0.0;
static int const   N_TARGETS = 4;
static int         target = 0;
static SL_DJstate  targets[N_TARGETS][N_DOFS+1];
static double      delta_t = 0.01;
static double      duration = 1.0;
static double      time_to_go;


// global functions 
extern "C" void
add_min_jerk_task( void );

// local functions
static int  init_min_jerk_task(void);
static int  run_min_jerk_task(void);
static int  change_min_jerk_task(void);

static int 
min_jerk_next_step (double x,double xd, double xdd, double t, double td, double tdd,
		    double t_togo, double dt,
		    double *x_next, double *xd_next, double *xdd_next);


// configuration variables for different homework questions
static bool SQUARE = false;
static bool USE_SPLINE = false;


/*****************************************************************************
******************************************************************************
Function Name	: add_min_jerk_task
Date		: Feb 1999
Remarks:

adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_min_jerk_task( void )
{
  int i, j;
  
  addTask("Min Jerk Task", init_min_jerk_task, 
	  run_min_jerk_task, change_min_jerk_task);

}    

/*****************************************************************************
******************************************************************************
  Function Name	: init_min_jerk_task
  Date		: Dec. 1997

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_min_jerk_task(void)
{
  int j, i;
  int ans;
  static int firsttime = TRUE;
  
  if (firsttime){
    firsttime = FALSE;
  }

  // fill all 4 targets with the default position
  for (j=0; j<N_TARGETS; j++)
    for (i=1; i<=N_DOFS; i++)
      targets[j][i] = joint_default_state[i];

  if (SQUARE) {
    // modify the targets to get a square using L_EB and L_SAA
    targets[0][L_SFE].th += 0.2;
    targets[1][L_SFE].th += 0.2;
    targets[1][L_SAA].th += 0.2;
    targets[2][L_SAA].th += 0.2;
  }
  else {
    // for reaching movement, only use the first target
    targets[0][R_SFE].th += 0.4;
    targets[0][R_SAA].th -= 0.4;
    targets[0][R_EB].th  -= 0.5;
  }

  // go to the default position target using inverse dynamics (ID)
  if (!go_target_wait_ID(targets[3])) 
    return FALSE;

  // this variable will increment each time we reach the current target
  target = 0;

  // ready to go
  ans = 999;
  while (ans == 999) {
    if (!get_int("Enter 1 to start or anthing else to abort ...",ans,&ans))
      return FALSE;
  }
  
  // only go when user really types the right thing
  if (ans != 1) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
	 start_time, task_servo_time);

  // start data collection
  scd();

  // time to go
  time_to_go = duration;

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: run_min_jerk_task
  Date		: Dec. 1997

  Remarks:

  run the task from the task servo: REAL TIME requirements!

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
run_min_jerk_task(void)
{
  int j, i;

  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  // compute the update for the desired states
  for (i=1; i<=N_DOFS; ++i) {
    min_jerk_next_step(joint_des_state[i].th,
		       joint_des_state[i].thd,
		       joint_des_state[i].thdd,
		       targets[target][i].th,
		       targets[target][i].thd,
		       targets[target][i].thdd,
		       time_to_go,
		       delta_t,
		       &(joint_des_state[i].th),
		       &(joint_des_state[i].thd),
		       &(joint_des_state[i].thdd));
  }

  // compute inverse dynamics torques
  SL_InvDynNE(joint_state,joint_des_state,endeff,&base_state,&base_orient);

  // decrement time to go
  time_to_go -= delta_t;
  if (time_to_go <= 0) {
    target++;
    if (!SQUARE || target == N_TARGETS) {
      freeze();
    }
    else {
      time_to_go = duration;
    }
  }

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_min_jerk_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_min_jerk_task(void)
{
  int    ivar;
  double dvar;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}


static void solve_min_jerk_poly(
	double x0, double dx0, double ddx0, double xf, double dxf, double ddxf, double T, double poly[6])
{
	double T2 = T * T;
	double T3 = T2 * T;
	double T4 = T3 * T;
	double T5 = T4 * T;
	poly[0] = x0;
	poly[1] = dx0;
	poly[2] = ddx0 / 2;
	poly[3] = (-12*dx0*T - 8*dxf*T - 3*ddx0*T2 + ddxf*T2 - 20*x0 + 20*xf)/(2*T3);
	poly[4] = (16*dx0*T + 14*dxf*T + 3*ddx0*T2 - 2*ddxf*T2 + 30*x0 - 30*xf)/(2*T4);
	poly[5] = (-6*dx0*T - 6*dxf*T - ddx0*T2 + ddxf*T2 - 12*x0 + 12*xf)/(2*T5);
}

// evaluate a polynomial using horner's rule, coefs are ascending
static double polyval(double const *poly, int deg, double t)
{
	if (deg == 0) return poly[0];
	return poly[0] + t * polyval(poly+1, deg-1, t);
}

// compute derivative of a polynomial
static void polyder(double const *p, int deg, double *dp)
{
	for (int i = 1; i <= deg; ++i) {
		dp[i-1] = i * p[i];
	}
}

/*!*****************************************************************************
 *******************************************************************************
\note  min_jerk_next_step
\date  April 2014
   
\remarks 

Given the time to go, the current state is updated to the next state
using min jerk splines

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]          x,xd,xdd : the current state, vel, acceleration
 \param[in]          t,td,tdd : the target state, vel, acceleration
 \param[in]          t_togo   : time to go until target is reached
 \param[in]          dt       : time increment
 \param[in]          x_next,xd_next,xdd_next : the next state after dt

 ******************************************************************************/

static int 
min_jerk_next_step(
	double x, double xd, double xdd, double t, double td, double tdd,
	double t_togo, double dt,
	double *x_next, double *xd_next, double *xdd_next)

{
	if (USE_SPLINE) {
		double p[6], dp[5], ddp[4];
		solve_min_jerk_poly(x, xd, xdd, t, td, tdd, t_togo, p);
		polyder(p, 5, dp);
		polyder(dp, 4, ddp);

		*x_next = polyval(p, 5, dt);
		*xd_next = polyval(dp, 4, dt);
		*xdd_next = polyval(ddp, 3, dt);
	}
	else {
		// dynamical system version
		double alpha = 25.0;
		double beta = 6.0;
		double tau = 1.0;
		*xdd_next = (alpha / tau) * (beta * (t - x) - xd);
		*xd_next = xd + dt * (*xdd_next);
		*x_next = x + dt * (*xd_next);
	}

	return TRUE;
}
