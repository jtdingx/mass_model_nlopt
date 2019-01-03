/**
MPCClass.h

Description:	Header file of MPCClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	Thu 02 Aug 2018 11:53:47 AM CEST
@Update:	Thu 02 Aug 2018 11:53:41 AM CEST
*/
#pragma once

#include "QP/QPBaseClass.h"
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>

#include <vector> 

#include "RobotPara/RobotParaClass.h"

using namespace Eigen;
using namespace std;


class MPCClass : public QPBaseClass
{
public:
	MPCClass();
	virtual ~MPCClass() {};
	
	//void FootStepNumberInputs(int footstepsnumber);
	void FootStepInputs(int footstepsnumber, double stepwidth, double steplength, double stepheight);
	
	void Initialize(double dt_mpc);
	
	
	
	void CoM_foot_trajection_generation_local(int i, Eigen::Matrix<double,18,1> estimated_state, Eigen::Vector3d _Rfoot_location_feedback, Eigen::Vector3d _Lfoot_location_feedback,bool _ISwalking);
	
	
	void solve_reactive_step_body_inclination_CoMz(); 
	void solve_reactive_step_body_inclination(); 
	void solve_reactive_step(); 
	
	
//  	int Indexfind(double goalvari, Eigen::VectorXd goalarray, int xyz);
	int Indexfind(double goalvari, int xyz);

	Eigen::MatrixXd Matrix_ps(Eigen::MatrixXd a, int nh, Eigen::MatrixXd cxps);
	Eigen::MatrixXd Matrix_pu(Eigen::MatrixXd a, Eigen::MatrixXd b, int nh, Eigen::MatrixXd cxpu);		
	void Solve();

	void Foot_trajectory_solve(int j_index, bool _ISwalking);
	

		
	
	// current state based on the past one and two actual sampling time;
	Vector3d XGetSolution_CoM_position(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);
	Vector3d XGetSolution_Foot_positionR(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);
	Vector3d XGetSolution_Foot_positionL(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);	
	Vector3d XGetSolution_body_inclination(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3);
	
	
	int Get_maximal_number(double dtx);
	
	int Get_maximal_number_reference();


	int is_sigular(int num);
	
	void File_wl();	
	
	std::string _robot_name;
	double _robot_mass;
	double _lift_height;
	double _tstep;
	
	int _method_flag;

	int _n_end_walking;

	//////////////////////////////// substitute the 
	int _j_period;

	Eigen::VectorXd _tx;		
	
protected: 


private:          
        //parameters declaration  
        int _footstepsnumber;
	

	
        Eigen::VectorXd _steplength, _stepwidth, _stepheight,_lift_height_ref;	
        Eigen::VectorXd _footx_ref, _footy_ref, _footz_ref;	
	Eigen::VectorXd _ts, _td;	
	
        double _dt; 
	Eigen::VectorXd _t;
	int _nsum;
	int _nT;	
	
	Eigen::RowVectorXd _zmpx_real, _zmpy_real;
	Eigen::RowVectorXd _comx, _comvx, _comax;
	Eigen::RowVectorXd _comy, _comvy, _comay;
	Eigen::RowVectorXd _comz, _comvz, _comaz;	
	Eigen::RowVectorXd _thetax, _thetavx, _thetaax;
	Eigen::RowVectorXd _thetay, _thetavy, _thetaay;	
	Eigen::RowVectorXd _thetaz, _thetavz, _thetaaz;	
	
	Eigen::RowVectorXd _torquex_real, _torquey_real;	
	
	// CoM+angular momentum state and contro input
	Eigen::MatrixXd _xk,_yk,_zk,_thetaxk,_thetayk;
	Eigen::RowVectorXd _x_vacc_k,_y_vacc_k,_z_vacc_k,_thetax_vacc_k,_thetay_vacc_k;
	
	std::vector<Eigen::Vector3d> _footsteps;
	
	// initial parameters for MPC
	double _hcom;
	Eigen::MatrixXd _Hcom1;
	double _g;
	Eigen::VectorXd _ggg;
	int    _nh;
	Eigen::MatrixXd _a,_b,_c,_cp,_cv,_ca;
	
	
	//vertical height constraints
	Eigen::VectorXd _z_max;
	Eigen::VectorXd _z_min;
	
	
	//footz
	Eigen::MatrixXd _Zsc;
	
	//predictive model
	Eigen::MatrixXd _pps,_ppu,_pvs,_pvu,_pas,_pau;
	
	
	Eigen::VectorXd _footx_real, _footy_real, _footz_real;
	Eigen::VectorXd _footx_real_next, _footy_real_next, _footz_real_next;
	Eigen::VectorXd _footx_real_next1, _footy_real_next1, _footz_real_next1;
	

	Eigen::MatrixXd _footxyz_real;
	Eigen::MatrixXd _Lfootxyz, _Rfootxyz;
	
	
	Eigen::VectorXd _footx_max, _footx_min,_footy_max,_footy_min;
	
	double _mass,  _rad,  _j_ini;
	
	// zmp-constraints
	Eigen::VectorXd _zmpx_ub, _zmpx_lb, _zmpy_ub, _zmpy_lb;
	
	// com-support range
	Eigen::VectorXd _comx_max,  _comx_min,  _comy_max, _comy_min;
	
	// angle range
	Eigen::VectorXd _thetax_max,  _thetax_min,  _thetay_max, _thetay_min;
	
	
	// torque range
	Eigen::VectorXd _torquex_max, _torquex_min, _torquey_max, _torquey_min; 
	
	// swing foot velocity constraints
	Eigen::VectorXd _footx_vmax, _footx_vmin,_footy_vmax,_footy_vmin;


	// solution preparation
	Eigen::MatrixXd _V_optimal;
	Eigen::VectorXd _flag;
	Eigen::VectorXd _flag_global;
	
	
	double _Rx,     _Ry,     _Rz;
	double _alphax, _alphay, _alphaz;
	double _beltax, _beltay, _beltaz;
	double _gamax,  _gamay,  _gamaz;
	double _Rthetax, _Rthetay;
	double _alphathetax, _alphathetay;
	double _beltathetax, _beltathetay;
	
	// time cost consumption
	Eigen::RowVectorXd _tcpu;
	Eigen::RowVectorXd _tcpu_iterative;
	Eigen::RowVectorXd _tcpu_prepara;
	Eigen::RowVectorXd _tcpu_prepara2;
	Eigen::RowVectorXd _tcpu_qp;	
	
	
// predictive model control_tracking with time_varying height		
	int _bjxx_1;
	int _bjxx;
	Eigen::VectorXd _t_f;
	
	int _bjx1, _bjx2, _mx;
	Eigen::VectorXd  _tnx;
	
	Eigen::MatrixXd _v_i, _VV_i;
	
	int _n_vis;
	
        Eigen::VectorXd _Lx_ref, _Ly_ref,_Lz_ref;
	
	int _Nt, _nstep;
	Eigen::MatrixXd _V_ini;
	
	int xxx, xxx1,xxx2;

        Eigen::VectorXd _fx, _fy;
	Eigen::VectorXd _fxx_global, _fyy_global;
	
	Eigen::MatrixXd _comx_center_ref, _comy_center_ref,_comz_center_ref;
	
	Eigen::MatrixXd _thetax_center_ref, _thetay_center_ref;      
	
	int _loop;
	
	// sqp model
	
	Eigen::MatrixXd _ppu_2, _pvu_2;
	
	Eigen::MatrixXd _WX, _WY, _WZ, _WthetaX, _WthetaY, _PHIX, _PHIY, _PHIZ, _Q_goal, _q_goal, _Q_goal1, _q_goal1;
	
	Eigen::MatrixXd A_unit, C_unit;	
	Eigen::MatrixXd _Sjx,_Sjy,_Sjz,_Sjthetax,_Sjthetay,_Sfx,_Sfy,_Sfz;
	
	// inequality constraints
	
	Eigen::MatrixXd _H_q_upx,_F_zmp_upx,_H_q_lowx,_F_zmp_lowx,_H_q_upy,_F_zmp_upy,_H_q_lowy,_F_zmp_lowy;
	
	Eigen::MatrixXd _H_h_upz,_F_h_upz,_H_h_lowz,_F_h_lowz, _delta_footz_up, _delta_footz_low;
	
	Eigen::MatrixXd _H_hacc_lowz,_F_hacc_lowz, _delta_footzacc_up;
	
	Eigen::MatrixXd _q_upx,_qq_upx,_q_lowx,_qq_lowx,_q_upy,_qq_upy,_q_lowy,_qq_lowy,_qq1_upx,_qq1_lowx,_qq1_upy,_qq1_lowy;
	
	Eigen::MatrixXd _t_upx,_tt_upx,_t_lowx,_tt_lowx,_t_upy,_tt_upy,_t_lowy,_tt_lowy,_tt1_upx,_tt1_lowx,_tt1_upy,_tt1_lowy;
	
	Eigen::MatrixXd _H_q_footx_up,_F_foot_upx,_H_q_footx_low,_F_foot_lowx,_H_q_footy_up,_F_foot_upy,_H_q_footy_low,_F_foot_lowy;
	
	Eigen::MatrixXd _Footvx_max,_Footvx_min,_Footvy_max,_Footvy_min, _footubxv,_footlbxv,_footubyv,_footlbyv;
	
	// equality constraints
	Eigen::MatrixXd _H_q_footz,_F_footz;
	
	Eigen::MatrixXd _h_h, _hhhx;
	Eigen::MatrixXd _a_hx, _a_hxx, _a_hy, _a_hyy;
	
	
	
	
	
	// zmp constraints
	Eigen::MatrixXd _Si;
	Eigen::MatrixXd _phi_i_x_up,_p_i_x_t_up,_del_i_x_up,_phi_i_x_low,_p_i_x_t_low,_del_i_x_low;
	Eigen::MatrixXd _phi_i_y_up,_p_i_y_t_up,_del_i_y_up,_phi_i_y_low,_p_i_y_t_low,_del_i_y_low;
	
	//foot location constraints
	Eigen::MatrixXd _Sfoot;
	
	//com support leg constraints
	Eigen::MatrixXd _S1;	
	
	Eigen::MatrixXd _V_inix;
	
// 	int xxx_vector=30;
	vector <Eigen::MatrixXd> ZMPx_constraints_offfline;
	vector <Eigen::MatrixXd> ZMPy_constraints_offfline;
	
	vector <Eigen::MatrixXd> ZMPx_constraints_half;
	vector <Eigen::MatrixXd> ZMPy_constraints_half;
	
	int xyz1;  //flag for find function 
	int xyz2;
	
	

	
	Eigen::RowVectorXd _Lfootx, _Lfooty,_Lfootz, _Lfootvx, _Lfootvy,_Lfootvz, _Lfootax, _Lfootay,_Lfootaz;
	Eigen::RowVectorXd _Rfootx, _Rfooty,_Rfootz, _Rfootvx, _Rfootvy,_Rfootvz, _Rfootax, _Rfootay,_Rfootaz;
        double _ry_left_right;
	
	
	////result CoM_foot_trajection_generation
	Eigen::MatrixXd _CoM_position_optimal;
	Eigen::MatrixXd _torso_angle_optimal;
	Eigen::MatrixXd _L_foot_optition_optimal;
	Eigen::MatrixXd _R_foot_optition_optimal;
	Eigen::MatrixXd _foot_location_optimal;
	
	

	
	double _Footx_global_relative,_Footy_global_relative;
	
	
	Eigen::MatrixXd _ZMPx_constraints_half2, _ZMPy_constraints_half2, _phi_i_x_up1, _phi_i_y_up1;
	vector <Eigen::MatrixXd> _phi_i_x_up_est, _phi_i_y_up_est;
	
	Eigen::MatrixXd CoMMM_ZMP_foot;
	
// 	Eigen::Matrix4d AAAaaa1;
	double _ZMP_ratio;
	
	double _t_end_walking;
// 	int _n_end_walking;
	
	Eigen::MatrixXd  _comy_matrix_inv;

// 	protected:
// 	void Rfooty_plan(int arg1);
};
