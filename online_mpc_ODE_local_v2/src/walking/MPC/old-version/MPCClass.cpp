/*****************************************************************************
MPCClass.cpp

Description:    source file of MPCClass

@Version:   1.0
@Author:    Chengxu Zhou (zhouchengxu@gmail.com)
@Release:   Thu 02 Aug 2018 12:33:23 PM CEST
@Update:    Thu 02 Aug 2018 12:33:19 PM CEST
*****************************************************************************/
#include "MPC/MPCClass.h"
#include "MPC/spline.h"
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>



using namespace Eigen;
using namespace std;


MPCClass::MPCClass()                    ///declaration function
	: QPBaseClass()
	, _robot_name("")
	, _robot_mass(0.0)
	, _lift_height(0.0)
	, _tstep(0.0)
	, _method_flag(0)
{
  
}

////////////////step parameters input============================================================
void MPCClass::FootStepInputs(int footstepsnumber, double stepwidth, double steplength, double stepheight)
{	
        _footstepsnumber = footstepsnumber;
	_steplength.setConstant(_footstepsnumber,steplength);
	_steplength(0) = 0;
	
	_stepwidth.setConstant(_footstepsnumber,stepwidth);
	_stepwidth(0) = _stepwidth(0)/2;
	
        _stepwidth1	= _stepwidth;
	
	_stepheight.setConstant(_footstepsnumber,stepheight);
//         
        _steplength(footstepsnumber-1) = 0;
        _steplength(footstepsnumber-2) = 0;
        _steplength(footstepsnumber-3) = 0;
        _steplength(footstepsnumber-4) = 0;	
//         _steplength(footstepsnumber-5) = 0;
	
	_lift_height_ref.setConstant(_footstepsnumber,_lift_height);

        _lift_height_ref(footstepsnumber-1) = 0;
        _lift_height_ref(footstepsnumber-2) = 0;	
        _lift_height_ref(footstepsnumber-3) = 0; 
/*        _lift_height_ref(footstepsnumber-4) = 0; */	
}

/////////////////////// initialize ============================================================
void MPCClass::Initialize(double dt_mpc)
{
/////////////===============================================================================	

	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
///////////////////////////////////////////////////////////////////////	
///////////////////////////////////////////////////////////////////////	
// 	_footsteps
 	// ==step loctions setup==
        _footx_ref.setZero(_footstepsnumber);
	_footy_ref.setZero(_footstepsnumber);
	_footz_ref.setZero(_footstepsnumber);
		
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _footx_ref(i) = _footx_ref(i-1) + _steplength(i-1);
	  _footy_ref(i) = _footy_ref(i-1) + (int)pow(-1,i-1)*_stepwidth(i-1);   
	  _footz_ref(i) = _footz_ref(i-1) + _stepheight(i-1);
	}
	
	_footy_ref1 =_footy_ref;
        
//       sampling time & step cycle

	
	_dt =dt_mpc;
	_ts.setConstant(_footstepsnumber,_tstep);
	_td = 0.2*_ts; 
	
	_tx.setZero(_footstepsnumber);
  	for (int i = 1; i < _footstepsnumber; i++) {
 	  _tx(i) = _tx(i-1) + _ts(i-1);
	  _tx(i) = round(_tx(i)/_dt)*_dt -0.00001;	  
	}	
// 	cout << _tx<<endl;
	      
	_t.setLinSpaced(round(_tx(_footstepsnumber-1)/_dt),_dt,_tx(_footstepsnumber-1));
	_nsum = _t.size();
	_nT = round(_ts(0)/_dt);
	

	
//////////////==============================state variable============================================
	
	_zmpx_real.setZero(_nsum); _zmpy_real.setZero(_nsum);
	_comx.setZero(_nsum); _comvx.setZero(_nsum); _comax.setZero(_nsum);
	_comy.setZero(_nsum); _comvy.setZero(_nsum); _comay.setZero(_nsum);
	_comz.setZero(_nsum); _comvz.setZero(_nsum); _comaz.setZero(_nsum);	
	_thetax.setZero(_nsum); _thetavx.setZero(_nsum); _thetaax.setZero(_nsum);
	_thetay.setZero(_nsum); _thetavy.setZero(_nsum); _thetaay.setZero(_nsum);
	_thetaz.setZero(_nsum); _thetavz.setZero(_nsum); _thetaaz.setZero(_nsum);

	
	_torquex_real.setZero(_nsum); _torquey_real.setZero(_nsum);


	 _CoM_position_optimal.setZero(3,_nsum);
	 _torso_angle_optimal.setZero(3,_nsum);
	 _L_foot_optition_optimal.setZero(3,_nsum);
	 _R_foot_optition_optimal.setZero(3,_nsum);
	 _foot_location_optimal.setZero(3,_footstepsnumber);
	 
	
	_xk.setZero(3,_nsum); _yk.setZero(3,_nsum); _zk.setZero(3,_nsum);
	_thetaxk.setZero(3,_nsum); _thetayk.setZero(3,_nsum);
	_x_vacc_k.setZero(_nsum); _y_vacc_k.setZero(_nsum); _z_vacc_k.setZero(_nsum); 
	_thetax_vacc_k.setZero(_nsum); _thetay_vacc_k.setZero(_nsum); 
	
///////================================================================================================	
        // ==initial parameters & matrix for MPC==
        _hcom = RobotParaClass::Z_C();

	_g = RobotParaClass::G();
	_ggg.setConstant(1, RobotParaClass::G());
	_nh = round(RobotParaClass::PreviewT()/_dt);	//0.02-0.09m stairs
	
	_Hcom1.setConstant(_nh,1,_hcom);
        	
	_a.setZero(3,3);
	_a << 1, _dt, pow(_dt,2)/2,    
	      0,   1,            _dt,
	      0,   0,              1;
	_b.setZero(3,1);
	_b << pow(_dt,3)/6,    
	      pow(_dt,2)/2,
	               _dt;	
	_c.setZero(1,3);
	_c << 1,0,(-1)*_hcom /_g;	
	
	_cp.setZero(1,3);
	_cp(0,0) = 1;
	_cv.setZero(1,3);
	_cv(0,1) = 1;
	_ca.setZero(1,3);
	_ca(0,2) = 1;	
		

	
	
        //footz refer: height of step
	_Zsc.setZero(_nsum,1);		

	xyz1 = 0;  //flag for find function 
	xyz2 = 1;
	

	_j_period = 0;
		

  	for (int i = 0; i < _nsum-1; i++) {	  
          Indexfind(_t(i),xyz1);
	  
	  _Zsc(i,0) = _footz_ref(_j_period);   
	  _j_period = 0; 
	}		

        _yk.topRows(1).setConstant(_footy_ref(0)); 
	_zk.topRows(1).setConstant(_hcom);
	
	
	_v_i.setZero(_nh,1);
	_VV_i.setZero(_nh, _nstep);	
	

	//predictive model
	_pps.setZero(_nh,3); _ppu.setZero(_nh,_nh);
	_pvs.setZero(_nh,3); _pvu.setZero(_nh,_nh);
	_pas.setZero(_nh,3); _pau.setZero(_nh,_nh);

        _pps = Matrix_ps(_a,_nh,_cp);
	_pvs = Matrix_ps(_a,_nh,_cv);
	_pas = Matrix_ps(_a,_nh,_ca);

	_ppu = Matrix_pu(_a,_b,_nh,_cp);
	_pvu = Matrix_pu(_a,_b,_nh,_cv);
	_pau = Matrix_pu(_a,_b,_nh,_ca);

	
	_footx_real.setZero(_footstepsnumber);  _footy_real.setZero(_footstepsnumber); _footz_real.setZero(_footstepsnumber);	
	_footxyz_real.setZero(3,_footstepsnumber);


	
	_footx_real_next.setZero(_nsum);  _footy_real_next.setZero(_nsum); _footz_real_next.setZero(_nsum);
	_footx_real_next1.setZero(_nsum);  _footy_real_next1.setZero(_nsum); _footz_real_next1.setZero(_nsum);
	
	
	_Lfootx.setZero(_nsum); _Lfooty.setConstant(_nsum,_stepwidth(0));_Lfootz.setZero(_nsum); _Lfootvx.setZero(_nsum); _Lfootvy.setZero(_nsum);_Lfootvz.setZero(_nsum); 
	_Lfootax.setZero(_nsum); _Lfootay.setZero(_nsum);_Lfootaz.setZero(_nsum);
	_Rfootx.setZero(_nsum); _Rfooty.setConstant(_nsum,-_stepwidth(0));_Rfootz.setZero(_nsum); _Rfootvx.setZero(_nsum); _Rfootvy.setZero(_nsum);_Rfootvz.setZero(_nsum); 
	_Rfootax.setZero(_nsum); _Rfootay.setZero(_nsum);_Rfootaz.setZero(_nsum);
	_ry_left_right = 0;
	
	
/////=========================================constraints========================================
	_ZMP_ratio = 0.8;
	
	if(_robot_name == "coman"){
	  //vertical height constraints
	  _z_max.setConstant(_nsum,0.1);
	  _z_min.setConstant(_nsum,-0.1);
	  
	  _rad = 0.1; 
	  
	  _footx_max.setConstant(1, 0.3);
	  _footx_min.setConstant(1, -0.2);	  
	  
	  /// zmp-constraints	
	  _zmpx_ub.setConstant(_nsum,0.07);  _zmpx_lb.setConstant(_nsum,-0.03);
	  _zmpy_ub.setConstant(_nsum,0.05); _zmpy_lb.setConstant(_nsum,-0.05);		  	  
	}
	else if (_robot_name == "bigman")
	{
	  //vertical height constraints
	  _z_max.setConstant(_nsum,0.1);
	  _z_min.setConstant(_nsum,-0.1);
	  
	  _rad = 0.2; 
	  
	  _footx_max.setConstant(1, 0.5);
	  _footx_min.setConstant(1, -0.2);	  
	  
	  /// zmp-constraints	
	  _zmpx_ub.setConstant(_nsum,(RobotParaClass::FOOT_LENGTH()/2+RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio);  
	  _zmpx_lb.setConstant(_nsum,-(RobotParaClass::FOOT_LENGTH()/2-RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio);
	  _zmpy_ub.setConstant(_nsum,RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio); 
	  _zmpy_lb.setConstant(_nsum,-RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio);	
	
	}
	else if (_robot_name == "cogimon")
        {
	  //vertical height constraints
	  _z_max.setConstant(_nsum,0.1);
	  _z_min.setConstant(_nsum,-0.1);
	  
	  _rad = 0.2; 
	  
	  _footx_max.setConstant(1, 0.4);
	  _footx_min.setConstant(1, -0.2);	  
	  
	  /// zmp-constraints	
	  _zmpx_ub.setConstant(_nsum,(RobotParaClass::FOOT_LENGTH()/2+RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio);  
	  _zmpx_lb.setConstant(_nsum,-(RobotParaClass::FOOT_LENGTH()/2-RobotParaClass::HIP_TO_ANKLE_X_OFFSET())*_ZMP_ratio);
	  _zmpy_ub.setConstant(_nsum,RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio); 
	  _zmpy_lb.setConstant(_nsum,-RobotParaClass::FOOT_WIDTH()/2*_ZMP_ratio);	
        } 
        else {
	  DPRINTF("Errorrrrrrrr for IK\n");}
	
	
	//vertical height constraints
		
	_mass = _robot_mass; 
		
	_j_ini = _mass* pow(_rad,2);
	
	
// 	/// zmp-constraints	

	_footy_max.setConstant(1, 2*RobotParaClass::HALF_HIP_WIDTH() + 0.2); 
	_footy_min.setConstant(1, RobotParaClass::HALF_HIP_WIDTH() - 0.03);
	
	_fx.setZero(1);
	_fy.setZero(1);

	_fxx_global.setZero(1);
	_fyy_global.setZero(1);	
	
	
	// com-support range=== not using here
	_comx_max.setConstant(1,0.06);
	_comx_min.setConstant(1,-0.04);  
	_comy_max.setConstant(1,0.6);  
	_comy_min.setConstant(1,0.02);
	
	// angle range
	_thetax_max.setConstant(1,10*M_PI/180);  
	_thetax_min.setConstant(1,-5*M_PI/180);
	_thetay_max.setConstant(1,10*M_PI/180);  
	_thetay_min.setConstant(1,-10*M_PI/180);
	
	// torque range
	_torquex_max.setConstant(1,80/_j_ini); 
	_torquex_min.setConstant(1,-60/_j_ini);
	_torquey_max.setConstant(1,80/_j_ini);  
	_torquey_min.setConstant(1,-80/_j_ini);	

	// swing foot velocity constraints	
	_footx_vmax.setConstant(1,3);
	_footx_vmin.setConstant(1,-2);
	_footy_vmax.setConstant(1,2); 
	_footy_vmin.setConstant(1,-2);	
	
	
	
	
////////////===================================================================	
///===========initiallize: preparation for MPC solution
	_nstep = 2;
	_Nt = 5*_nh + 3*_nstep;	
	
	// sulotion preparation		
	_V_optimal.setZero(_Nt, _nsum);	
	_Lx_ref.setZero(_nstep); _Ly_ref.setZero(_nstep); _Lz_ref.setZero(_nstep);
	_V_ini.setZero(_Nt,1);
	_comx_center_ref.setZero(_nh,1);
	_comy_center_ref.setZero(_nh,1);
	_comz_center_ref.setZero(_nh,1);
	
	_thetax_center_ref.setZero(_nh,1); 
	_thetay_center_ref.setZero(_nh,1);	
	
        
// 	 store n_vis
	_flag.setZero(_nsum);
	
	_flag_global.setZero(_nsum);
	


	if(_robot_name == "coman"){
// // 	parameters for objective function======================	
// /////////////////////////////////////local coordinate generation///////////////
// /////////////////////////////////////
         //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline
	 _Rx = 1;           _Ry = 1;            _Rz =1;
	_alphax = 1;       _alphay = 1;        _alphaz = 100; 
	_beltax = 5000;   _beltay = 10;        _beltaz = 20000000;
	_gamax =  10000000; _gamay = 10000000;  _gamaz = 200;
	_Rthetax = 1; _Rthetay = 1;
	_alphathetax =1; _alphathetay = 1;
	_beltathetax = 10; _beltathetay = 10;
	  
	}
	else if(_robot_name  == "bigman"){
// // 	parameters for objective function======================	
// /////////////////////////////////////local coordinate generation///////////////
         //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline
	 _Rx = 1;           _Ry = 1;            _Rz =1;
	_alphax = 1;       _alphay = 10;        _alphaz = 100; 
	_beltax = 100;   _beltay = 1000;        _beltaz = 20000000;
	_gamax =  50000000; _gamay = 1000000000;  _gamaz = 200;
	_Rthetax = 1; _Rthetay = 1;
	_alphathetax =1; _alphathetay = 1;
	_beltathetax = 1000; _beltathetay = 1000;
// 
// 

	}
	else if (_robot_name == "cogimon"){
	  //////// for methx ==2:reactive step + body inclination + height variance: for flat ground walking and up-down stairs: offline
	  _Rx = 1;           _Ry = 1;            _Rz =1;
	  _alphax = 1;       _alphay = 1;        _alphaz = 100; 
	  _beltax = 100;   _beltay = 100;        _beltaz = 20000000;
	  _gamax =  50000000; _gamay = 100000000;  _gamaz = 200;
	  _Rthetax = 1; _Rthetay = 1;
	  _alphathetax =1; _alphathetay = 1;
	  _beltathetax = 1000; _beltathetay = 1000;
        } 
	else
	{DPRINTF("Errorrrrrrrr for IK\n");}
	
	
	
	
	
	
	
	// time cost consumption
	_tcpu.setZero(_nsum);
	_tcpu_iterative.setZero(_nsum);
	_tcpu_prepara.setZero(_nsum);
	_tcpu_prepara2.setZero(_nsum);
	_tcpu_qp.setZero(_nsum);

	
	_pvu_2 = _pvu.transpose()*_pvu;
	_ppu_2 = _ppu.transpose()*_ppu;

	
	_loop = 2;
	


	
	
///////////////////////////////////////////////////////////////=====================================//////////////////////////////////////
//////////// next code block just run once	
	  A_unit.setIdentity(_nh,_nh);
	  C_unit.setIdentity(_nstep,_nstep);		
	

	  // optimization objective function 	
	  _WX.setZero(_nh,_nh);
	  _WY.setZero(_nh,_nh);
	  _WZ.setZero(_nh,_nh);
	  _WthetaX.setZero(_nh,_nh);
	  _WthetaY.setZero(_nh,_nh);
	  _PHIX.setZero(_nstep,_nstep);
	  _PHIY.setZero(_nstep,_nstep);
	  _PHIZ.setZero(_nstep,_nstep);
	  _Q_goal.setZero(_Nt,_Nt);
	  _q_goal.setZero(_Nt,1);
	  _Q_goal1.setZero(_Nt,_Nt);
	  _q_goal1.setZero(_Nt,1);	
	  
	  _WX = _Rx*0.5 * A_unit + _alphax*0.5 * _pvu_2 + _beltax*0.5 * _ppu_2;	  
	  _WY = _Ry/2 * A_unit + _alphay/2 * _pvu_2 + _beltay/2 * _ppu_2;
	  _WZ = _Rz/2 * A_unit + _alphaz/2 * _pvu_2 + _beltaz/2 * _ppu_2;  
	  _WthetaX = _Rthetax/2 * A_unit + _alphathetax/2 * _pvu_2 + _beltathetax/2 * _ppu_2;
	  _WthetaY = _Rthetay/2 * A_unit + _alphathetay/2 * _pvu_2 + _beltathetay/2 * _ppu_2;
	  _PHIX  = _gamax/2 * C_unit;
	  _PHIY  = _gamay/2 * C_unit;
	  _PHIZ  = _gamaz/2 * C_unit;
		  
	  _Q_goal.block< Dynamic, Dynamic>(0, 0,_nh, _nh) = _WX;
	  _Q_goal.block< Dynamic, Dynamic>(_nh, _nh,_nh, _nh) = _WY;
	  _Q_goal.block< Dynamic, Dynamic>(2*_nh, 2*_nh,_nh, _nh) = _WZ;
	  _Q_goal.block< Dynamic, Dynamic>(3*_nh, 3*_nh,_nh, _nh) = _WthetaX;
	  _Q_goal.block< Dynamic, Dynamic>(4*_nh, 4*_nh,_nh, _nh) = _WthetaY;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh, 5*_nh,_nstep,_nstep) = _PHIX;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh+_nstep, 5*_nh+_nstep,_nstep, _nstep) = _PHIY;
	  _Q_goal.block< Dynamic, Dynamic>(5*_nh+2*_nstep, 5*_nh+2*_nstep,_nstep, _nstep) = _PHIZ;	  
	
	  _Q_goal1 = 2 * _Q_goal;	
	
	
	  // constraints
	  _Sjx.setZero(_nh,_Nt);
	  _Sjy.setZero(_nh,_Nt);
	  _Sjz.setZero(_nh,_Nt);
	  _Sjthetax.setZero(_nh,_Nt);
	  _Sjthetay.setZero(_nh,_Nt);
	  _Sjx.block< Dynamic, Dynamic>(0, 0,_nh, _nh) = A_unit;
	  _Sjy.block< Dynamic, Dynamic>(0, _nh,_nh, _nh) = A_unit;
	  _Sjz.block< Dynamic, Dynamic>(0, 2*_nh,_nh, _nh) = A_unit;	
	  _Sjthetax.block< Dynamic, Dynamic>(0, 3*_nh,_nh, _nh) = A_unit;
	  _Sjthetay.block< Dynamic, Dynamic>(0, 4*_nh,_nh, _nh) = A_unit;
	  
	  // ZMP boundary preparation
	  _H_q_upx.setZero(_nh,_Nt);
	  _F_zmp_upx.setZero(_nh,1);
	  _H_q_lowx.setZero(_nh,_Nt);
	  _F_zmp_lowx.setZero(_nh,1);
	  _H_q_upy.setZero(_nh,_Nt);
	  _F_zmp_upy.setZero(_nh,1);
	  _H_q_lowy.setZero(_nh,_Nt);
	  _F_zmp_lowy.setZero(_nh,1);

	  _phi_i_x_up.setZero(_Nt,_Nt);
	  _p_i_x_t_up.setZero(_Nt,_nh);
	  _del_i_x_up.setZero(1,_nh);
	  _phi_i_x_low.setZero(_Nt,_Nt);
	  _p_i_x_t_low.setZero(_Nt,_nh);
	  _del_i_x_low.setZero(1,_nh);
	  _phi_i_y_up.setZero(_Nt,_Nt);
	  _p_i_y_t_up.setZero(_Nt,_nh);
	  _del_i_y_up.setZero(1,_nh);
	  _phi_i_y_low.setZero(_Nt,_Nt);
	  _p_i_y_t_low.setZero(_Nt,_nh);
	  _del_i_y_low.setZero(1,_nh);	  
	  
	  // angle boundary preparation
	  _q_upx.setZero(_nh,_Nt);
	  _qq_upx.setZero(_nh,1);
	  _q_lowx.setZero(_nh,_Nt);
	  _qq_lowx.setZero(_nh,1);
	  _q_upy.setZero(_nh,_Nt);
	  _qq_upy.setZero(_nh,1);
	  _q_lowy.setZero(_nh,_Nt);
	  _qq_lowy.setZero(_nh,1);
	  
	  _qq1_upx.setZero(_nh,1);
	  _qq1_lowx.setZero(_nh,1);
	  _qq1_upy.setZero(_nh,1);
	  _qq1_lowy.setZero(_nh,1);	  

	  // torque bondary preparation
	  _t_upx.setZero(_nh,_Nt);
	  _tt_upx.setZero(_nh,1);
	  _t_lowx.setZero(_nh,_Nt);
	  _tt_lowx.setZero(_nh,1);
	  _t_upy.setZero(_nh,_Nt);
	  _tt_upy.setZero(_nh,1);
	  _t_lowy.setZero(_nh,_Nt);
	  _tt_lowy.setZero(_nh,1);
	  
	  _tt1_upx.setZero(_nh,1);
	  _tt1_lowx.setZero(_nh,1);
	  _tt1_upy.setZero(_nh,1);
	  _tt1_lowy.setZero(_nh,1);
	  
	  // CoM height boundary preparation
	  _H_h_upz.setZero(_nh,_Nt);
	  _F_h_upz.setZero(_nh,1);
	  _H_h_lowz.setZero(_nh,_Nt);
	  _F_h_lowz.setZero(_nh,1);
	  _delta_footz_up.setZero(_nh,1);
	  _delta_footz_low.setZero(_nh,1);

	  // CoM height acceleration boundary preparation
	  _H_hacc_lowz.setZero(_nh,_Nt);
	  _F_hacc_lowz.setZero(_nh,1);
	  _delta_footzacc_up.setZero(_nh,1);	  

	  
	  //swing foot velocity constraints
	  _Footvx_max.setZero(1,_Nt);
	  _Footvx_min.setZero(1,_Nt);
	  _Footvy_max.setZero(1,_Nt);
	  _Footvy_min.setZero(1,_Nt);
	  _footubxv.setZero(1,1);
	  _footlbxv.setZero(1,1);
	  _footubyv.setZero(1,1);
	  _footlbyv.setZero(1,1);
	  
	  // foot vertical location-equality constraints
	  _H_q_footz.setZero(1, _Nt);
	  _F_footz.setZero(1, 1);
	  
	  // CoMZ height-equality constraints
	  _h_h.setZero(_nh, _Nt);
	  _F_footz.setZero(_nh, 1);	  

	  // body inclination-equality constraints
	  _a_hx.setZero(_nh, _Nt);
	  _a_hxx.setZero(_nh, 1);
	  _a_hy.setZero(_nh, _Nt);
	  _a_hyy.setZero(_nh, 1);
	  
	  
	  // foot location constraints
	  _Sfoot.setZero(1,2);
	  _Sfoot(0,0) = -1;
	  _Sfoot(0,1) = 1;
	  
	  _S1.setZero(1,_nh);
	  _S1(0,0) = 1;	  
	  
	// offline calulated the ZMP constraints coefficient==================================
	  //////////initiallize vector
// 	vector <Eigen::MatrixXd> x_offline1(_nh)  ;
	
	Eigen::MatrixXd mmmm;
	mmmm.setZero(_Nt, _Nt);
	vector <Eigen::MatrixXd> x_offline1(_nh)  ;
	for (int j=0;j<_nh; j++)
	{
	  x_offline1[j]= mmmm;
	}

	ZMPx_constraints_offfline = x_offline1;
	ZMPy_constraints_offfline = x_offline1;	
	
	_phi_i_x_up_est = x_offline1;
	_phi_i_y_up_est = x_offline1;
	

	
	Eigen::MatrixXd mmmmn;
	mmmmn.setZero(_nh, _Nt);
	
	vector <Eigen::MatrixXd> x_offline2(_nh)  ;
	for (int j=0;j<_nh; j++)
	{
	  x_offline2[j]= mmmmn;
	}	

		

	ZMPx_constraints_half = x_offline2;	
	ZMPy_constraints_half = x_offline2;
		
	for(int jxx=1; jxx<=_nh; jxx++)
	{
	  _Si.setZero(1,_nh);
	  _Si(0,jxx-1) = 1;
	  // ZMP constraints	      
		 
         ZMPx_constraints_offfline[jxx-1] = (_Si * _ppu * _Sjx).transpose() * _Si * _pau * _Sjz - (_Si * _pau * _Sjx).transpose() * _Si * _ppu * _Sjz;
	 ZMPx_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
	  
	  
         ZMPy_constraints_offfline[jxx-1] = (_Si * _ppu * _Sjy).transpose() * _Si * _pau * _Sjz - (_Si * _pau * _Sjy).transpose() * _Si * _ppu * _Sjz;
	 ZMPy_constraints_half[jxx-1] = - (_Si).transpose() * _Si * _pau * _Sjz;
	      
	}
        
   
  
  
  
  /////////// initialize each variable
        _bjxx_1 = 0; 
	_bjxx = 0; 
	_t_f.setZero(_nh);
	_bjx1 = 0;
	_bjx2 = 0;
        _mx = 0;
        _tnx.setZero(5);
  
	
	
        _Footx_global_relative =0;
  
        _Footy_global_relative =0;
	  _v_i.setZero(_nh,1);
	  _VV_i.setZero(_nh, _nstep);	
	_n_vis =0; 
	xxx = 0; xxx1=0;
	

	  // foot location constraints: be careful that the step number is change: so should be intialized in each whole loop
	  _H_q_footx_up.setZero(_nstep,_Nt);
	  _F_foot_upx.setZero(_nstep,1);
	  _H_q_footx_low.setZero(_nstep,_Nt);
	  _F_foot_lowx.setZero(_nstep,1);
	  _H_q_footy_up.setZero(_nstep,_Nt);
	  _F_foot_upy.setZero(_nstep,1);
	  _H_q_footy_low.setZero(_nstep,_Nt);
	  _F_foot_lowy.setZero(_nstep,1);


	  // boundary initialzation
	  _Sfx.setZero(_nstep,_Nt);
	  _Sfy.setZero(_nstep,_Nt);
	  _Sfz.setZero(_nstep,_Nt);	  
   
	  
         _Si.setZero(1,_nh);	  

         _ZMPx_constraints_half2.setZero(_Nt, _nh);	 
	 _ZMPy_constraints_half2.setZero(_Nt, _nh);
	 _phi_i_x_up1.setZero(_Nt,_Nt);
	 _phi_i_y_up1.setZero(_Nt,_Nt);
	 
	 CoMMM_ZMP_foot.setZero(28,_nsum);
	 
	 ///QP initiallize
	 if (_method_flag ==0)
	 {
	  int nVars = _Nt;
	  int nEqCon = 1+3*_nh;
	  int nIneqCon = 5*_nh + 4*_nstep+4;  
	  resizeQP(nVars, nEqCon, nIneqCon);		   
	}
	 else if (_method_flag ==1)
	 {
	  int nVars = _Nt;
	  int nEqCon = 1+_nh;
	  int nIneqCon = 13*_nh + 4*_nstep +4;
	  resizeQP(nVars, nEqCon, nIneqCon);	   
	}
	 else
	 {
	  int nVars = _Nt;
	  int nEqCon = 1;
	  int nIneqCon = 15*_nh + 4*_nstep +4;
	  resizeQP(nVars, nEqCon, nIneqCon);		   
	}
// 	cout<<_tx<<endl;
	

	_t_end_walking = _tx(_footstepsnumber-1)- 9*_tstep/4;
// 	cout << _t_end_walking<<endl;
	_n_end_walking = round(_t_end_walking/_dt);
// 	cout <<_n_end_walking<<endl;
	
	_comy_matrix_inv.setZero(4,1);
	
}




/////////////////////// local coordinate CoM solution---modified---------------------------------
void MPCClass::CoM_foot_trajection_generation_local(int i, Eigen::Matrix<double,18,1> estimated_state, Eigen::Vector3d _Rfoot_location_feedback, Eigen::Vector3d _Lfoot_location_feedback, bool _stopwalking)
{
 
       if (i<_n_end_walking)   ///////normal walking
       {
	 

	  /// modified the footy_min
	  if (i==(round(2*_ts(1)/_dt))+1) ///update the footy_limit
	  {
	      _footy_min.setConstant(1, RobotParaClass::FOOT_WIDTH()+0.01);
	  }
  //         cout<<"_footy_min:"<<_footy_min<<endl;


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
  // 	================ iterative calulcation: predicitve control_tracking with time-varying height+angular momentum	  
	    clock_t t_start,t_start1, t_start2, t_start3, t_start4,t_finish,t_finish1;	  
	    
	    /// run once
	    ///////////////////////////////////////////////////////
	    ////////////////////////////////////// time_clock0////////////////////////
	    t_start = clock();
		    
  // 	  Indexfind((i-1)*_dt,_tx,xyz1);
	    Indexfind((i-1)*_dt,xyz1);	  
	    
	    _bjxx_1 = _j_period+1;
	    _j_period = 0;
	    
	    
	    Indexfind(i*_dt,xyz1);
	    _bjxx = _j_period+1;  //coincidence with matlab 
	    _j_period = 0;	  
	    
	    // com_center_ref = ZMP_center_ref = v_i*f + V_i*L_ref
	    //solve the following steps: 1.5s may couve 2 or three  steps, so  one/two following steps
	    
	    _t_f.setLinSpaced(_nh,(i+1)*_dt, (i+_nh)*_dt);
	    
	    Indexfind(_t_f(0),xyz1);
	    _bjx1 = _j_period+1;
	    _j_period = 0;
	    
	    Indexfind(_t_f(_nh-1),xyz1);
	    _bjx2 = _j_period+1;
	    _j_period = 0;	  
	    
	    
	    ////================================================================
	  /// judge if stop walking enable
	  if(_stopwalking)  
	  {
	    
	    for (int i_t = _bjx1; i_t < _footstepsnumber; i_t++) {
	      _steplength(i_t) = 0;	
	      _footx_ref(i_t) = _footx_ref(i_t-1) + _steplength(i_t-1); 
	      
	      _stepwidth1(0) = _stepwidth(0)/2;
	      _footy_ref1(i_t) = _footy_ref1(i_t-1) + (int)pow(-1,i_t-1)*_stepwidth1(i_t-1);	
	    }	  
	  }
	  else if (_bjx2==_footstepsnumber-4)	  
	  {
/*	    for (int i_t = _bjx1; i_t < _footstepsnumber; i_t++) {	    
	      _stepwidth1(0) = _stepwidth(0)/2;
	      _footy_ref1(i_t) = _footy_ref1(i_t-1) + (int)pow(-1,i_t-1)*_stepwidth1(i_t-1);	
	    } */ 
	  }
	    
	    
	    
	    
	    _mx = _bjx2 - _bjx1 +1;
	    
	    for (int j=1;j<_mx; j++)
	    {
	      Indexfind(_tx(_bjx1+j-1),xyz2);
	      _tnx(j-1) = _j_period; 
	      _j_period = 0;
	      
	    }
	    
  // 	  _v_i.setZero(_nh,1);
  // 	  _VV_i.setZero(_nh, _nstep);
	    _v_i.setZero();
	    _VV_i.setZero();
	    // be careful that the position is from 0;;;;;         
	    if (fabs(_tnx(0) - _nT) <=0.00001)
	    {
	      _n_vis =2;
	      for (int jjj = 1; jjj <= _mx; jjj++)
	      {
		if (jjj == 1)
		{
		  xxx = _tnx(0);
		  _VV_i.block< Dynamic,1>(0, 0, xxx, 1).setOnes();
		}
		else
		{	
		  xxx1 = _nh-_tnx(0);
		  _VV_i.block< Dynamic,1>(xxx+1, 1, xxx1, 1).setOnes();		
		}
	      }	    
	    }
	    else
	    {	
	      xxx = _tnx(0);
	      _v_i.block< Dynamic,1>(0, 0, xxx, 1).setOnes();

	      if (abs(_mx - 2) <=0.00001)
	      {
		_n_vis = 1;
		_VV_i.block< Dynamic,1>(xxx, 0, _nh -xxx, 1).setOnes();
	      }
	      else
	      {
		_n_vis = 2;
		xxx2 = _nh -_tnx(1);
		_VV_i.block< Dynamic,1>(xxx, 0, _nT, 1).setOnes();
		_VV_i.block< Dynamic,1>(_tnx(1), 1, xxx2, 1).setOnes();	      	      
	      }
	    }
	    

    
	    _flag(i-1,0)= _n_vis;
	    
	    _flag_global(i-1,0) = _bjxx;
	    _flag_global(i,0) = _bjx1;	  

  ////////////////////////////////////////////////////////////////////////////////////	  
  //============================================================//	  
	  
  
	    
  //////////////////// relative state: 	  
	    ///// pass the actual stage into the control loop	  	  
	      
	    if (i>1)
	    {

	      _Footx_global_relative = _xk(0,i-1) + _fxx_global(0);
	      _Footy_global_relative = _yk(0,i-1) + _fyy_global(0);	    
	    }
	    // current foot location

  /*	  _fx.setZero(1);
	    _fy.setZero(1);*/	  
	    _fx(0) =0;
	    _fy(0) = 0;	  
	    
	    _fxx_global(0) = _footx_real(_bjxx-1);
	    _fyy_global(0) = _footy_real(_bjxx-1);
	    

	    if (_n_vis ==1)
	    {
	      _Lx_ref(0) = _footx_ref(_bjx2-1) - _fxx_global(0);
	      _Ly_ref(0) = _footy_ref(_bjx2-1) - _fyy_global(0);
	      _Lz_ref(0) = _footz_ref(_bjx2-1);
	      _Lx_ref(1) = 0;
	      _Ly_ref(1) = 0;
	      _Lz_ref(1) = 0;	    
	    }
	    else
	    {
	      _Lx_ref(0) = _footx_ref(_bjx2-2) - _fxx_global(0);
	      _Ly_ref(0) = _footy_ref(_bjx2-2) - _fyy_global(0);
	      _Lz_ref(0) = _footz_ref(_bjx2-2);
	      _Lx_ref(1) = _footx_ref(_bjx2-1) - _fxx_global(0);
	      _Ly_ref(1) = _footy_ref(_bjx2-1) - _fyy_global(0);
	      _Lz_ref(1) = _footz_ref(_bjx2-1);	    
	    }

	    
/*	  if(_stopwalking)  
	  {
	    
	  _comy_center_ref = _v_i*_fy + _VV_i*_Ly_ref/2;
	  
	  }
	  else if (_bjx2>= _footstepsnumber-5)
	  {	  
	    _comy_center_ref = _v_i*_fy + _VV_i*_Ly_ref/2;

	  }
	  else
	  {
	    _comy_center_ref = _v_i*_fy + _VV_i*_Ly_ref;	  
	  }
		    
	*/    
	    
	    
	  // com_center_ref
	    _comx_center_ref = _v_i*_fx + _VV_i*_Lx_ref;
  	  _comy_center_ref = _v_i*_fy + _VV_i*_Ly_ref;
	    _comz_center_ref = _Zsc.block< Dynamic,1>(i, 0,_nh, 1) + _Hcom1;

	    
	    /// hot start
	    if (i==1)
	    {
	      _V_ini(5*_nh,0) = _footx_ref(1);
	      _V_ini(5*_nh+1,0) = _footy_ref(1);	    
	    }
	    else
	    {
	      _V_ini.topRows(5*_nh) = _V_optimal.block< Dynamic, 1>(0, i-2, 5*_nh, 1);
	      if (_n_vis > _flag(i-2))
	      { 
		_V_ini(_Nt -1-1,0) = _V_optimal(5*_nh+6-1, i-2);
		_V_ini(_Nt -3-1,0) = _V_optimal(5*_nh+6-2-1, i-2);
		_V_ini(_Nt -5-1,0) = _V_optimal(5*_nh+6-4-1, i-2); 
	      }
	      else
	      {
		if (_n_vis < _flag(i-2))
		{ 
		  _V_ini(_Nt-1,0) = _V_optimal(5*_nh+6-1, i-2);
		  _V_ini(_Nt -1-1,0) = _V_optimal(5*_nh+6-2-1, i-2);
		  _V_ini(_Nt -2-1,0) = _V_optimal(5*_nh+6-4-1, i-2); 
		}	   
		else
		{
		  if (_n_vis ==1)
		  { 
		    _V_ini(_Nt-1,0) = _V_optimal(5*_nh+6-1, i-2);
		    _V_ini(_Nt -1-1,0) = _V_optimal(5*_nh+6-2-1, i-2);
		    _V_ini(_Nt -2-1,0) = _V_optimal(5*_nh+6-4-1, i-2); 
		  }
		  else
		  {
		    _V_ini.bottomRows(6) = _V_optimal.block< Dynamic, 1>(5*_nh+6-6, i-2, 6, 1);
		  }
		}
	      }
	    }
	      
	    

  // 	      // relative state switch	      
	    if (i>1)
	    {
	      if (_flag_global(i-2,0) < _flag_global(i-1,0) )
	      {
		/// reference relative state switch
	      
		_xk(0,i-1) = _Footx_global_relative - _fxx_global(0); 
		_yk(0,i-1) = _Footy_global_relative - _fyy_global(0);
		
	      }

	    }

	    

	    ///////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////////////////////////////////////============================================================//////////////////////////////////////////

	    // foot location constraints: be careful that the step number is change: so should be intialized in each whole loop
	    _H_q_footx_up.setZero();
	    _F_foot_upx.setZero();
	    _H_q_footx_low.setZero();
	    _F_foot_lowx.setZero();
	    _H_q_footy_up.setZero();
	    _F_foot_upy.setZero();
	    _H_q_footy_low.setZero();
	    _F_foot_lowy.setZero();
	    

	    // boundary initialzation
	    _Sfx.setZero();
	    _Sfy.setZero();
	    _Sfz.setZero();
	    
	    if (_n_vis ==1)
	    {
	      _Sfx(0,5*_nh) = 1;
	      _Sfy(0,5*_nh+_nstep) = 1;
	      _Sfz(0,5*_nh+2*_nstep) = 1;
	    }  
	    else
	    {
	      _Sfx(0,5*_nh) = 1;
	      _Sfx(1,5*_nh+1) = 1;
	      _Sfy(0,5*_nh+_nstep) = 1;
	      _Sfy(1,5*_nh+_nstep+1) = 1;
	      _Sfz(0,5*_nh+2*_nstep) = 1;	 
	      _Sfz(1,5*_nh+2*_nstep+1) = 1;	
	      
	    }
		
		  
	    

	    //SQP MOdels	 
	      _q_goal.block< Dynamic, 1>(0, 0,_nh, 1) = _alphax * _pvu.transpose() * _pvs * _xk.col(i-1) + _beltax * _ppu.transpose() * _pps * _xk.col(i-1) - _beltax * _ppu.transpose() * _comx_center_ref;
	      _q_goal.block< Dynamic, 1>(_nh, 0,_nh, 1) = _alphay * _pvu.transpose() * _pvs * _yk.col(i-1) + _beltay * _ppu.transpose() * _pps * _yk.col(i-1) - _beltay * _ppu.transpose() * _comy_center_ref;
	      _q_goal.block< Dynamic, 1>(2*_nh, 0,_nh, 1) = _alphaz * _pvu.transpose() * _pvs * _zk.col(i-1) + _beltaz * _ppu.transpose() * _pps * _zk.col(i-1) - _beltaz * _ppu.transpose() * _comz_center_ref;
	      _q_goal.block< Dynamic, 1>(3*_nh, 0,_nh, 1) = _alphathetax * _pvu.transpose() * _pvs * _thetaxk.col(i-1) + _beltathetax * _ppu.transpose() * _pps * _thetaxk.col(i-1) - _beltathetax * _ppu.transpose() * _thetax_center_ref;
	      _q_goal.block< Dynamic, 1>(4*_nh, 0,_nh, 1) = _alphathetay * _pvu.transpose() * _pvs * _thetayk.col(i-1) + _beltathetay * _ppu.transpose() * _pps * _thetayk.col(i-1) - _beltathetay * _ppu.transpose() * _thetay_center_ref;
	      _q_goal.block< Dynamic, 1>(5*_nh, 0,_nstep, 1) = -_gamax * _Lx_ref;
	      _q_goal.block< Dynamic, 1>(5*_nh+_nstep, 0,_nstep, 1) = -_gamay * _Ly_ref;
	      _q_goal.block< Dynamic, 1>(5*_nh+2*_nstep, 0,_nstep, 1) = -_gamaz * _Lz_ref;

	      
		
	    ///// the following block only run once in each loop   
	      for(int jxx=1; jxx<=_nh; jxx++)
	      {
  // 	      _Si.setZero(1,_nh);
		_Si.setZero();
		_Si(0,jxx-1) = 1;
		// ZMP constraints
		// x-ZMP upper boundary

  // 	      Eigen::MatrixXd _p_i_x_t_up1 = (((_Si * _pps * _xk.col(i-1)).transpose() *_Si*_pau*_Sjz + (_Si*_pas*_zk.col(i-1)).transpose()*_Si*_ppu*_Sjx + _g*_Si*_ppu*_Sjx - ((_Si * _pps * _zk.col(i-1)).transpose() *_Si*_pau*_Sjx + _Si * _pas * _xk.col(i-1)* _Si* _ppu* _Sjz) + _Zsc.row(i+jxx-1)*_Si*_pau*_Sjx - ((_Si * _pas * _zk.col(i-1)).transpose() *_Si*_VV_i*_Sfx + (_Si * _v_i * _fx).transpose() *_Si*_pau*_Sjz) - _g*_Si*_VV_i*_Sfx - _zmpx_ub.row(i+jxx-1)*_Si*_pau*_Sjz).transpose());                                         
		_p_i_x_t_up.col(jxx-1) = _mass * (((_Si * _pps * _xk.col(i-1)).transpose() *_Si*_pau*_Sjz + (_Si*_pas*_zk.col(i-1)).transpose()*_Si*_ppu*_Sjx + _g*_Si*_ppu*_Sjx - ((_Si * _pps * _zk.col(i-1)).transpose() *_Si*_pau*_Sjx + _Si * _pas * _xk.col(i-1)* _Si* _ppu* _Sjz) + _Zsc.row(i+jxx-1)*_Si*_pau*_Sjx - ((_Si * _pas * _zk.col(i-1)).transpose() *_Si*_VV_i*_Sfx + (_Si * _v_i * _fx).transpose() *_Si*_pau*_Sjz) - _g*_Si*_VV_i*_Sfx - _zmpx_ub.row(i+jxx-1)*_Si*_pau*_Sjz).transpose()) - (_j_ini * _Si*_pau * _Sjthetay).transpose();
		
  // 	      Eigen::MatrixXd _del_i_x_up1 = ((_Si * _pps * _xk.col(i-1)).transpose() *_Si*_pas*_zk.col(i-1) + _g*_Si * _pps * _xk.col(i-1) - (_Si * _pas * _xk.col(i-1)).transpose() *_Si*_pps*_zk.col(i-1) + (_Si * _pas * _xk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - (_Si * _v_i * _fx).transpose() *_Si*_pas*_zk.col(i-1) - _g *_Si * _v_i * _fx - _zmpx_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1) - _g *_zmpx_ub.row(i+jxx-1));
		_del_i_x_up.col(jxx-1) = _mass * ((_Si * _pps * _xk.col(i-1)).transpose() *_Si*_pas*_zk.col(i-1) + _g*_Si * _pps * _xk.col(i-1) - (_Si * _pas * _xk.col(i-1)).transpose() *_Si*_pps*_zk.col(i-1) + (_Si * _pas * _xk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - (_Si * _v_i * _fx).transpose() *_Si*_pas*_zk.col(i-1) - _g *_Si * _v_i * _fx - _zmpx_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1) - _g *_zmpx_ub.row(i+jxx-1)) - _j_ini * _Si*_pas * _thetayk.col(i-1);


		// x-ZMP low boundary
  // 	      _phi_i_x_low = _phi_i_x_up;
		_p_i_x_t_low.col(jxx-1) = (_p_i_x_t_up.col(jxx-1).transpose() + _mass * _zmpx_ub.row(i+jxx-1)*_Si*_pau*_Sjz - _mass * _zmpx_lb.row(i+jxx-1)*_Si*_pau*_Sjz).transpose();	      
		_del_i_x_low.col(jxx-1) = _del_i_x_up.col(jxx-1) +_mass*_zmpx_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1)+  _mass * _g*_zmpx_ub.row(i+jxx-1) - _mass*_zmpx_lb.row(i+jxx-1)*_Si*_pas*_zk.col(i-1)-_mass * _g*_zmpx_lb.row(i+jxx-1);
		
		// y-ZMP upper boundary
    
  // 	      Eigen::MatrixXd _p_i_y_t_up1 = (((_Si * _pps * _yk.col(i-1)).transpose() *_Si*_pau*_Sjz + (_Si*_pas*_zk.col(i-1)).transpose()*_Si*_ppu*_Sjy + _g*_Si*_ppu*_Sjy - ((_Si * _pps * _zk.col(i-1)).transpose() *_Si*_pau*_Sjy + _Si * _pas * _yk.col(i-1)* _Si* _ppu* _Sjz) + _Zsc.row(i+jxx-1)*_Si*_pau*_Sjy - ((_Si * _pas * _zk.col(i-1)).transpose() *_Si*_VV_i*_Sfy + (_Si * _v_i * _fy).transpose() *_Si*_pau*_Sjz) - _g*_Si*_VV_i*_Sfy - _zmpy_ub.row(i+jxx-1)*_Si*_pau*_Sjz).transpose());                                         
		_p_i_y_t_up.col(jxx-1) = _mass * (((_Si * _pps * _yk.col(i-1)).transpose() *_Si*_pau*_Sjz + (_Si*_pas*_zk.col(i-1)).transpose()*_Si*_ppu*_Sjy + _g*_Si*_ppu*_Sjy - ((_Si * _pps * _zk.col(i-1)).transpose() *_Si*_pau*_Sjy + _Si * _pas * _yk.col(i-1)* _Si* _ppu* _Sjz) + _Zsc.row(i+jxx-1)*_Si*_pau*_Sjy - ((_Si * _pas * _zk.col(i-1)).transpose() *_Si*_VV_i*_Sfy + (_Si * _v_i * _fy).transpose() *_Si*_pau*_Sjz) - _g*_Si*_VV_i*_Sfy - _zmpy_ub.row(i+jxx-1)*_Si*_pau*_Sjz).transpose()) + (_j_ini * _Si*_pau * _Sjthetax).transpose();
		
  //	      Eigen::MatrixXd _del_i_y_up1 = ((_Si * _pps * _yk.col(i-1)).transpose() *_Si*_pas*_zk.col(i-1) + _g*_Si * _pps * _yk.col(i-1) - (_Si * _pas * _yk.col(i-1)).transpose() *_Si*_pps*_zk.col(i-1) + (_Si * _pas * _yk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - (_Si * _v_i * _fy).transpose() *_Si*_pas*_zk.col(i-1) - _g *_Si * _v_i * _fy - _zmpy_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1) - _g *_zmpy_ub.row(i+jxx-1));
		_del_i_y_up.col(jxx-1) = _mass * ((_Si * _pps * _yk.col(i-1)).transpose() *_Si*_pas*_zk.col(i-1) + _g*_Si * _pps * _yk.col(i-1) - (_Si * _pas * _yk.col(i-1)).transpose() *_Si*_pps*_zk.col(i-1) + (_Si * _pas * _yk.col(i-1)).transpose() *_Zsc.row(i+jxx-1) - (_Si * _v_i * _fy).transpose() *_Si*_pas*_zk.col(i-1) - _g *_Si * _v_i * _fy - _zmpy_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1) - _g *_zmpy_ub.row(i+jxx-1)); + _j_ini * _Si*_pas * _thetaxk.col(i-1);	      
	      
		// y-ZMP low boundary
	      _phi_i_y_low = _phi_i_y_up;  
		_p_i_y_t_low.col(jxx-1) = (_p_i_y_t_up.col(jxx-1).transpose() + _mass * _zmpy_ub.row(i+jxx-1)*_Si*_pau*_Sjz - _mass * _zmpy_lb.row(i+jxx-1)*_Si*_pau*_Sjz).transpose();	      
		_del_i_y_low.col(jxx-1) = _del_i_y_up.col(jxx-1) +_mass*_zmpy_ub.row(i+jxx-1)*_Si*_pas*_zk.col(i-1)+  _mass * _g*_zmpy_ub.row(i+jxx-1) - _mass*_zmpy_lb.row(i+jxx-1)*_Si*_pas*_zk.col(i-1)-_mass * _g*_zmpy_lb.row(i+jxx-1);	      	      	     

		
		
		//angle range constraints
		_q_upx.row(jxx-1) = _Si* _ppu* _Sjthetax;
		_q_lowx.row(jxx-1) = -_q_upx.row(jxx-1);	         

		_qq1_upx.row(jxx-1) = _Si* _pps* _thetaxk.col(i-1)-_thetax_max;
		
		_q_upy.row(jxx-1) = _Si* _ppu* _Sjthetay;
		_q_lowy.row(jxx-1) = -_q_upy.row(jxx-1);	
		
		_qq1_upy.row(jxx-1) = _Si* _pps* _thetayk.col(i-1)-_thetay_max;
    

		//torque range constraints
		_t_upx.row(jxx-1) = _Si* _pau* _Sjthetax;
		_t_lowx.row(jxx-1) = -_t_upx.row(jxx-1);
		
		_tt1_upx.row(jxx-1) = _Si* _pas* _thetaxk.col(i-1)-_torquex_max;
		
		_t_upy.row(jxx-1) = _Si* _pau* _Sjthetay;
		_t_lowy.row(jxx-1) = -_t_upy.row(jxx-1);	 
		
		_tt1_upy.row(jxx-1) = _Si* _pas* _thetayk.col(i-1)-_torquey_max;
		
		
		// body height constraints
		_H_h_upz.row(jxx-1) = _Si* _ppu* _Sjz;
		_H_h_lowz.row(jxx-1) = -_Si* _ppu* _Sjz;   	
		_delta_footz_up.row(jxx-1) = _Si*_pps*_zk.col(i-1) - _comz_center_ref.row(jxx-1) - _z_max.row(i+jxx-1);
		
		// body height acceleration constraints	      
		_H_hacc_lowz.row(jxx-1) = -_Si* _pau* _Sjz;   
		_delta_footzacc_up.row(jxx-1) = _Si*_pas*_zk.col(i-1) + _ggg;
	      }
			      
    
    
    
    
  ///       only one-time caluation	  
	    _ZMPx_constraints_half2 = (_VV_i* _Sfx).transpose();
	    _ZMPy_constraints_half2 = (_VV_i* _Sfy).transpose();
	    

	    for(int jxx=1; jxx<=_nh; jxx++)
	    {
	      // ZMP constraints
	      // x-ZMP upper boundary
	      _phi_i_x_up1 = ZMPx_constraints_offfline[jxx-1] + _ZMPx_constraints_half2 * ZMPx_constraints_half[jxx-1];	      
	      _phi_i_x_up_est[jxx-1] = _mass * (_phi_i_x_up1 + _phi_i_x_up1.transpose())/2;
    
	      // y-ZMP upper boundary
	      _phi_i_y_up1 = ZMPy_constraints_offfline[jxx-1] + _ZMPy_constraints_half2 * ZMPy_constraints_half[jxx-1];
	      _phi_i_y_up_est[jxx-1] = _mass * (_phi_i_y_up1 + _phi_i_y_up1.transpose())/2;   
	    }

	    
	    // constraints: only once 
	    _Footvx_max = _Sfx.row(0);
	    _Footvx_min = -_Sfx.row(0);
	    _Footvy_max = _Sfy.row(0);
	    _Footvy_min = -_Sfy.row(0);	  
		    
	    ///////////// equality equation	    
	    //equality constraints
	    _H_q_footz = _Sfz.row(0);
	    
	    _h_h = _ppu * _Sjz;
	    
	    _a_hx = _ppu * _Sjthetax;
	    _a_hy = _ppu * _Sjthetay;
	    
	    
		    
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
	    // SEQUENCE QUADARTIC PROGRAMMING: lOOP_until the maximal loops reaches	
	    // SEQUENCE QUADARTIC PROGRAMMING: lOOP_until the maximal loops reaches	
	    t_start1 = clock();
    
	  // hot start
  /*	  
	    _V_ini.setZero(_Nt,1);*/
	    
	    for (int xxxx=1; xxxx <= _loop; xxxx++)
	    {	
    
	      _q_goal1 = _Q_goal1 * _V_ini + _q_goal;	  
		      
  ///////////// inequality equation   
	      t_start2 = clock();
	      
	      /// time consuming process	    
	      
	      for(int jxx=1; jxx<=_nh; jxx++)
	      {
		// ZMP constraints
		_phi_i_x_up = _phi_i_x_up_est[jxx-1];
		_H_q_upx.row(jxx-1) = (2*_phi_i_x_up*_V_ini + _p_i_x_t_up.col(jxx-1)).transpose(); 
		_F_zmp_upx.row(jxx-1) = -((_V_ini.transpose() * _phi_i_x_up + _p_i_x_t_up.col(jxx-1).transpose()) * _V_ini + _del_i_x_up.col(jxx-1)); 

		// x-ZMP low boundary
		_phi_i_x_low = _phi_i_x_up;	      	      
		_H_q_lowx.row(jxx-1) = (- _p_i_x_t_low.col(jxx-1) + _p_i_x_t_up.col(jxx-1)).transpose() - _H_q_upx.row(jxx-1);  
		_F_zmp_lowx.row(jxx-1) = (_p_i_x_t_low.col(jxx-1)-_p_i_x_t_up.col(jxx-1)).transpose() * _V_ini + _del_i_x_low.col(jxx-1) -  _del_i_x_up.col(jxx-1)-_F_zmp_upx.row(jxx-1); 
		
		
		// y-ZMP upper boundary	      
		_phi_i_y_up = _phi_i_y_up_est[jxx-1];	      
		_H_q_upy.row(jxx-1) = (2*_phi_i_y_up*_V_ini + _p_i_y_t_up.col(jxx-1)).transpose(); 
		_F_zmp_upy.row(jxx-1) = -((_V_ini.transpose() * _phi_i_y_up + _p_i_y_t_up.col(jxx-1).transpose()) * _V_ini + _del_i_y_up.col(jxx-1)); 
		
		// y-ZMP low boundary
		_phi_i_y_low = _phi_i_y_up;  	      	      
		_H_q_lowy.row(jxx-1) = (- _p_i_y_t_low.col(jxx-1) + _p_i_y_t_up.col(jxx-1)).transpose() - _H_q_upy.row(jxx-1); 
		_F_zmp_lowy.row(jxx-1) = (_p_i_y_t_low.col(jxx-1)-_p_i_y_t_up.col(jxx-1)).transpose() * _V_ini + _del_i_y_low.col(jxx-1) - _del_i_y_up.col(jxx-1)-_F_zmp_upy.row(jxx-1);      	      
	    

		
		//angle range constraints
		_qq_upx.row(jxx-1) = -(_q_upx.row(jxx-1)* _V_ini + _qq1_upx.row(jxx-1));
		_qq_lowx.row(jxx-1) = _thetax_max - _thetax_min - _qq_upx.row(jxx-1);	      

				
		_qq_upy.row(jxx-1) = -(_q_upy.row(jxx-1)* _V_ini + _qq1_upy.row(jxx-1));
		_qq_lowy.row(jxx-1) = _thetay_max - _thetay_min - _qq_upy.row(jxx-1);	    

		//torque range constraints	      
		_tt_upx.row(jxx-1) = -(_t_upx.row(jxx-1)* _V_ini +  _tt1_upx.row(jxx-1));
		_tt_lowx.row(jxx-1) = _torquex_max - _torquex_min - _tt_upx.row(jxx-1);	
		
		_tt_upy.row(jxx-1) = -(_t_upy.row(jxx-1)* _V_ini +  _tt1_upy.row(jxx-1));
		_tt_lowy.row(jxx-1) = _torquey_max - _torquey_min - _tt_upy.row(jxx-1);		      
		
		// body height constraints	      
		_F_h_upz.row(jxx-1) = -(_H_h_upz.row(jxx-1)*_V_ini + _delta_footz_up.row(jxx-1));
		_F_h_lowz.row(jxx-1) = _z_max.row(i+jxx-1) - _z_min.row(i+jxx-1)-_F_h_upz.row(jxx-1);	      
		
		// body height acceleration constraints	      
		_F_hacc_lowz.row(jxx-1) = (-_H_hacc_lowz.row(jxx-1)*_V_ini + _delta_footzacc_up.row(jxx-1));	      	      
	      }

	      
	      
	      
	      
	      
	      t_start3 = clock();
	      
	      // foot location constraints
	      if (_n_vis == 1)  //one next steo
	      {
		_H_q_footx_up.row(0) = _Sfx.row(0);
		_F_foot_upx.row(0) = -(_H_q_footx_up.row(0) * _V_ini -_fx - _footx_max); 
		_H_q_footx_low.row(0) = -_Sfx.row(0);
		_F_foot_lowx.row(0) = (_H_q_footx_up.row(0) * _V_ini -_fx - _footx_min);
		
		
		
		// footy location constraints
		if (_bjxx % 2 == 0) //odd
		{
		  _H_q_footy_up.row(0) = _Sfy.row(0);
		  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini -_fy + _footy_min); 
		  _H_q_footy_low.row(0) = -_Sfy.row(0);
		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini -_fy + _footy_max);		
		}
		else
		{
		  _H_q_footy_up.row(0) = _Sfy.row(0);
		  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini -_fy - _footy_max); 
		  _H_q_footy_low.row(0) = -_Sfy.row(0);
		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini -_fy - _footy_min);			
		}	 
		
		
	      }
	      else   //two next steps
	      {
		_H_q_footx_up.row(0) = _Sfx.row(0);
		_F_foot_upx.row(0) = -(_H_q_footx_up.row(0) * _V_ini -_fx - _footx_max); 
		_H_q_footx_low.row(0) = -_Sfx.row(0);
		_F_foot_lowx.row(0) = (_H_q_footx_up.row(0) * _V_ini -_fx - _footx_min);
		
		// footy location constraints
		if (_bjxx % 2 == 0) //odd
		{
		  _H_q_footy_up.row(0) = _Sfy.row(0);
		  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini -_fy + _footy_min); 
		  _H_q_footy_low.row(0) = -_Sfy.row(0);
		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini -_fy + _footy_max);		
		}
		else
		{
		  _H_q_footy_up.row(0) = _Sfy.row(0);
		  _F_foot_upy.row(0) = -(_H_q_footy_up.row(0) * _V_ini -_fy - _footy_max); 
		  _H_q_footy_low.row(0) = -_Sfy.row(0);
		  _F_foot_lowy.row(0) = (_H_q_footy_up.row(0) * _V_ini -_fy - _footy_min);			
		}
		
		// the next two steps 
		_H_q_footx_up.row(1) = _Sfoot* _Sfx;
		_F_foot_upx.row(1) = -(_H_q_footx_up.row(1) * _V_ini - _footx_max); 	      
		_H_q_footx_low.row(1) = -_Sfoot* _Sfx;
		_F_foot_lowx.row(1) = (_H_q_footx_up.row(1) * _V_ini - _footx_min);
		
		// footy location constraints
		if (_bjxx % 2 == 0) //odd
		{
		  _H_q_footy_up.row(1) = _Sfoot* _Sfy;
		  _F_foot_upy.row(1) = -(_H_q_footy_up.row(1) * _V_ini - _footy_max); 
		  _H_q_footy_low.row(1) = -_H_q_footy_up.row(1);
		  _F_foot_lowy.row(1) = (_H_q_footy_up.row(1) * _V_ini - _footy_min);				
		}
		else
		{
		  _H_q_footy_up.row(1) = _Sfoot* _Sfy;
		  _F_foot_upy.row(1) = -(_H_q_footy_up.row(1) * _V_ini + _footy_min); 
		  _H_q_footy_low.row(1) = -_H_q_footy_up.row(1);
		  _F_foot_lowy.row(1) = (_H_q_footy_up.row(1) * _V_ini + _footy_max);		
		}	      	     	      
	      }

	      
	      //swing foot veloctiy boundary
	      if (i ==1)
	      {
		_footubxv = -(_Sfx.row(0) * _V_ini - _footx_max - _footx_real_next.row(i+_nT-2));
		_footlbxv = (_Sfx.row(0) * _V_ini - _footx_min - _footx_real_next.row(i+_nT-2));
		_footubyv = -(_Sfy.row(0) * _V_ini - _footy_max - _footy_real_next.row(i+_nT-2));
		_footlbyv = (_Sfy.row(0) * _V_ini - _footy_min - _footy_real_next.row(i+_nT-2));	      
	      }
	      else
	      {
		if (fabs(i*_dt - _tx(_bjxx-1))<=0.01)
		{	
  /*		_Footvx_max.setZero(1,_Nt);  _Footvx_min.setZero(1,_Nt);  _Footvy_max.setZero(1,_Nt);  _Footvy_min.setZero(1,_Nt);
		  
		  _footubxv.setZero(1,1); _footlbxv.setZero(1,1); _footubyv.setZero(1,1); _footlbyv.setZero(1,1);	*/	
		  _Footvx_max.setZero();  _Footvx_min.setZero();  _Footvy_max.setZero();  _Footvy_min.setZero();
		  
		  _footubxv.setZero(); _footlbxv.setZero(); _footubyv.setZero(); _footlbyv.setZero();			
		}
		else
		{
		  _footubxv = -(_Sfx.row(0) * _V_ini - _footx_vmax*_dt - _footx_real_next.row(i+_nT-2));
		  _footlbxv = (_Sfx.row(0) * _V_ini - _footx_vmin*_dt - _footx_real_next.row(i+_nT-2));
		  _footubyv = -(_Sfy.row(0) * _V_ini - _footy_vmax*_dt - _footy_real_next.row(i+_nT-2));
		  _footlbyv = (_Sfy.row(0) * _V_ini - _footy_vmin*_dt - _footy_real_next.row(i+_nT-2));		
		}
	      }


	    ///////////// equality equation	    
	    //equality constraints
	    _F_footz = _Sfz.row(0)*_V_ini - _Lz_ref.row(0);	    

	    _hhhx = _h_h*_V_ini + _pps * _zk.col(i-1) - _Hcom1;	

	    _a_hxx = _a_hx * _V_ini + _pps * _thetaxk.col(i-1);
	    _a_hyy = _a_hy * _V_ini + _pps * _thetayk.col(i-1);
		    

	      
  //////////////////////////////===========////////////////////////////////////////////////////	    
  // 	    // quadratic program GetSolution
	      t_start4 = clock();
	      
	      if (_method_flag ==0)
	      {
		solve_reactive_step(); ///merely the reactive steps
	      }
	      else
	      {
		if (_method_flag ==1)
		{
		  solve_reactive_step_body_inclination(); /// the reactive steps + body inclination
		}
		else
		{
		  solve_reactive_step_body_inclination_CoMz(); /// the reactive steps + body inclination + height variance
		}
	      }

	      
	      t_finish1 = clock();
	      _tcpu_prepara(0,i-1) = (double)(t_finish1 - t_start2)/CLOCKS_PER_SEC ;
    
	      _tcpu_prepara2(0,i-1) = (double)(t_finish1 - t_start3)/CLOCKS_PER_SEC ;
	      _tcpu_qp(0,i-1) = (double)(t_finish1 - t_start4)/CLOCKS_PER_SEC ;	    
	    }

	    
	    
	    t_finish = clock();	  
  /////////////////////////////////////////////////////////
  /////////===================================================%%%%%	  
	  // results postprocessed:	  
	    if (_n_vis == 1)
	    {
	      _V_optimal.block< Dynamic, 1>(0, i-1, 5*_nh, 1)  = _V_ini.topRows(5*_nh);
	      _V_optimal.block< Dynamic, 1>(5*_nh, i-1, 2, 1)  = _V_inix.setConstant(2,1,_V_ini(5*_nh,0));
	      _V_optimal.block< Dynamic, 1>(5*_nh+2, i-1, 2, 1)  = _V_inix.setConstant(2,1,_V_ini(5*_nh+1,0));
	      _V_optimal.block< Dynamic, 1>(5*_nh+4, i-1, 2, 1)  = _V_inix.setConstant(2,1,_V_ini(5*_nh+2,0));
	      
	    }
	    else
	    {
	      _V_optimal.col(i-1) = _V_ini;
	    }
	    
	    //next step location
	    _x_vacc_k.col(i-1) = _V_ini.row(0);
	    _footx_real.row(_bjxx) = _V_ini.row(5*_nh);
	    _footx_real_next.row(i+_nT -1) = _V_ini.row(5*_nh);
	    _xk.col(i) = _a * _xk.col(i-1)  + _b* _x_vacc_k.col(i-1);
	    _comx(0,i)=_xk(0,i); _comvx(0,i) = _xk(1,i); _comax(0,i)=_xk(2,i); 	  
	    
	    _y_vacc_k.col(i-1) = _V_ini.row(0+_nh);
	    _footy_real.row(_bjxx) = _V_ini.row(5*_nh + _nstep);
	    _footy_real_next.row(i+_nT -1) = _V_ini.row(5*_nh + _nstep);
	    _yk.col(i) = _a * _yk.col(i-1)  + _b* _y_vacc_k.col(i-1);
	    _comy(0,i)=_yk(0,i); _comvy(0,i) = _yk(1,i); _comay(0,i)=_yk(2,i); 
	    
	    
	    
	    _comx(0,i) +=  _fxx_global(0);
	    _comy(0,i) +=  _fyy_global(0);	  
	    _footx_real(_bjxx) = _footx_real(_bjxx) + _fxx_global(0);
	    _footy_real(_bjxx) = _footy_real(_bjxx) + _fyy_global(0);	  	  
	    _footx_real_next1(i+_nT -1) = _footx_real_next(i+_nT -1) + _fxx_global(0);	  
	    _footy_real_next1(i+_nT -1) = _footy_real_next(i+_nT -1) + _fyy_global(0);

	    
	    
	    _z_vacc_k.col(i-1) = _V_ini.row(0+2*_nh);
	    _footz_real.row(_bjxx) = _V_ini.row(5*_nh + 2*_nstep);
	    _footz_real_next.row(i+_nT -1) = _V_ini.row(5*_nh + 2*_nstep);	  
	    _zk.col(i) = _a * _zk.col(i-1)  + _b* _z_vacc_k.col(i-1);
	    _comz(0,i)=_zk(0,i); _comvz(0,i) = _zk(1,i); _comaz(0,i)=_zk(2,i); 

	    _thetax_vacc_k.col(i-1) = _V_ini.row(0+3*_nh);
	    _thetaxk.col(i) = _a * _thetaxk.col(i-1)  + _b* _thetax_vacc_k.col(i-1);
	    _thetax(0,i) = _thetaxk(0,i); _thetavx(0,i) = _thetaxk(1,i); _thetaax(0,i)=_thetaxk(2,i); 	


	    _thetay_vacc_k.col(i-1) = _V_ini.row(0+4*_nh);
	    _thetayk.col(i) = _a * _thetayk.col(i-1)  + _b* _thetay_vacc_k.col(i-1);
	    _thetay(0,i) = _thetayk(0,i); _thetavy(0,i) = _thetayk(1,i); _thetaay(0,i)=_thetayk(2,i); 	
	    
	    
	    // reference relative state
	    
  // /*  	  estimated_state(0,0) =  _comx(0,i-1) - _fxx_global(0); 
  //  	  estimated_state(1,0) =  _comvx(0,i-1);
  // 	  estimated_state(2,0) =  _comax(0,i-1); */  
	    

	    // reference relative state
  // 	  estimated_state(3,0) =  _comy(0,i-1) - _fyy_global(0);  
  // 	  estimated_state(4,0) =  _comvy(0,i-1);  
  // 	  estimated_state(5,0) =  _comay(0,i-1);
		  
	    
	    /// /// relative state to the actual foot lcoation: very good
	    if (_bjxx % 2 == 0)  // odd : left support
	    {
	      estimated_state(0,0) =  estimated_state(0,0) - _Lfoot_location_feedback(0);
	      estimated_state(3,0) =  estimated_state(3,0) - _Lfoot_location_feedback(1);	 	    
	    }
	    else
	    {
	      estimated_state(0,0) =  estimated_state(0,0) - _Rfoot_location_feedback(0);
	      estimated_state(3,0) =  estimated_state(3,0) - _Rfoot_location_feedback(1);	  
	    }
	    
  // 	      cout << "i:"<<i<<endl;

  //////============================================================================================================================	      
  ////////////////////////////// state modified:====================================================================================

	    if (_method_flag <2)
	    {
	      estimated_state(6,0) =  _comz(0,i-1);  
	      estimated_state(7,0) =  _comvz(0,i-1);  
	      estimated_state(8,0) =  _comaz(0,i-1);
	    }	  
	  
	    
	    if (_method_flag <=0)
	    {
	      estimated_state(9,0) =  _thetax(0,i-1);  
	      estimated_state(10,0) =  _thetavx(0,i-1);  
	      estimated_state(11,0) =  _thetaax(0,i-1);	  
	      estimated_state(12,0) =  _thetay(0,i-1);  
	      estimated_state(13,0) =  _thetavy(0,i-1);  
	      estimated_state(14,0) =  _thetaay(0,i-1);	    
	    }
		
		
		
		
  /*	  ///////////////================state feedback=========================================////
	    /// model0================COMX+COMY feedback

	    _xk(0,i) = (estimated_state(0,0)+2*_xk(0,i))/3;             
	    _xk(1,i) = (estimated_state(1,0)+2*_xk(1,i))/3; 
	    _xk(2,i) = (estimated_state(2,0)+2*_xk(2,i))/3;
	    

	    if (i<round(2*_ts(1)/_dt))
	    {
		
	    _yk(0,i) = (estimated_state(3,0)+24*_yk(0,i))/25; _yk(1,i) = (estimated_state(4,0)+24*_yk(1,i))/25;
	    _yk(2,i) = (estimated_state(5,0)+24*_yk(2,i))/25;		    
	      
	    }
	    else
	    {
		
	    _yk(0,i) = (estimated_state(3,0)+2*_yk(0,i))/3; _yk(1,i) = (estimated_state(4,0)+2*_yk(1,i))/3;
	    _yk(2,i) = (estimated_state(5,0)+2*_yk(2,i))/3;		    
	    }

	    ///model1===================COMX+COMY + thetax+ thetay feedback
	    _thetaxk(0,i) = (estimated_state(9,0)+1*_thetaxk(0,i))/2;
	    _thetaxk(1,i) = (estimated_state(10,0)+1*_thetaxk(1,i))/2; 
	    _thetaxk(2,i) = (estimated_state(11,0)+1*_thetaxk(2,i))/2;	  
	    _thetayk(0,i) = (estimated_state(12,0)+1*_thetayk(0,i))/2; 
	    _thetayk(1,i) = (estimated_state(13,0)+1*_thetayk(1,i))/2; 
	    _thetayk(2,i) = (estimated_state(14,0)+1*_thetayk(2,i))/2;	
	    
	    //model2====================COMX+COMY + thetax+ thetay + CoMz feedback
	    _zk(0,i) = (estimated_state(6,0)+2*_zk(0,i))/3; 
	    _zk(1,i) = (estimated_state(7,0)+2*_zk(1,i))/3; 
	    _zk(2,i) = (estimated_state(8,0)+2*_zk(2,i))/3;*/	
	    
	    
	    
	    
	    
    

  ////////////////===============================================================================================	  
	  /// next two sample time:	actually the preictive value is not reliable  
	    _comx(0,i+1) = _comx(0,i) + _dt * _comvx(0,i); 	  	  
	    _comy(0,i+1) = _comy(0,i) + _dt * _comvy(0,i); 	 
	    _comz(0,i+1) = _comz(0,i) + _dt * _comvz(0,i); 	 
	    _thetax(0,i+1) = _thetax(0,i)+ _dt * _thetavx(0,i); 	
	    _thetay(0,i+1) = _thetay(0,i)+ _dt * _thetavy(0,i); 	
	    
	    
	    _torquex_real.col(i) = _j_ini * _thetaax.col(i);
	    _torquey_real.col(i) = _j_ini * _thetaay.col(i);
	    
	    _zmpx_real(0,i) = _comx(0,i) - (_comz(0,i) - _Zsc(i,0))/(_comaz(0,i)+_g)*_comax(0,i) - _j_ini * _thetavy(0,i)/(_mass * (_g + _comaz(0,i)));
	    _zmpy_real(0,i) = _comy(0,i) - (_comz(0,i) - _Zsc(i,0))/(_comaz(0,i)+_g)*_comay(0,i) + _j_ini * _thetavx(0,i)/(_mass * (_g + _comaz(0,i)));
	    
	    t_finish = clock();
	    _tcpu(0,i-1) = (double)(t_finish - t_start)/CLOCKS_PER_SEC ;
	    _tcpu_iterative(0,i-1) = (double)(t_finish - t_start1)/CLOCKS_PER_SEC ;
	    
	    _footxyz_real.row(0) = _footx_real.transpose();
	    _footxyz_real.row(1) = _footy_real.transpose();	  
	    _footxyz_real.row(2) = _footz_real.transpose();
	  
	  
	  if (i>=1)
	  {	  
	    _Rfootx(0) = _Rfootx(1);
	    _Lfootx(0) = _Lfootx(1);
	    _Rfooty(0) = _Rfooty(1);
	    _Lfooty(0) = _Lfooty(1);
	    _comx(0) = _comx(1);	
	    _comy(0) = _comy(1);
	    _comz(0) = _comz(1);	  
	  }	 
	  
	 
	 
	 
	 
	 
      }
       else
       {
	 if(i <= _n_end_walking+round(_tstep/_dt/2)){
	   Eigen::Matrix<double, 4, 1> _comy_temp;
	   _comy_temp.setZero();
	   _comy_temp(0) = _comy(0,_n_end_walking-2);
	   _comy_temp(1) = _comy(0,_n_end_walking-1);
	   _comy_temp(2) = 0;
	   _comy_temp(3) = 0;
	      
	   
	   Eigen::Matrix<double, 1, 4> _t_temp;
	   

	   _t_temp(0) = pow(i-_n_end_walking+1, 3);
	   _t_temp(1) = pow(i-_n_end_walking+1, 2);
	   _t_temp(2) = pow(i-_n_end_walking+1, 1);
	   _t_temp(3) = pow(i-_n_end_walking+1, 0);
	   
	   Eigen::Matrix<double, 1, 4> _t_temp1;

	   _t_temp1(0) = pow(i-_n_end_walking+2, 3);
	   _t_temp1(1) = pow(i-_n_end_walking+2, 2);
	   _t_temp1(2) = pow(i-_n_end_walking+2, 1);
	   _t_temp1(3) = pow(i-_n_end_walking+2, 0);	   
	   
	   
	   Eigen::Matrix<double, 4, 4> _comy_matrix;
	   
	   int ix_temp1 = -1;
	   int ix_temp2 = 0;	   
	   int ix_temp3 = round(_tstep/_dt)/2+1;
	   
	  
	  Eigen::Matrix4d AAA_inv;
	  
	  double abx1, abx2, abx3, abx4;
	  abx1 = ((ix_temp1 - ix_temp2)*pow(ix_temp1 - ix_temp3, 2));
	  abx2 = ((ix_temp1 - ix_temp2)*pow(ix_temp2 - ix_temp3, 2));
	  abx3 =(pow(ix_temp1 - ix_temp3, 2)*pow(ix_temp2 - ix_temp3, 2));
	  abx4 = ((ix_temp1 - ix_temp3)*(ix_temp2 - ix_temp3));
	  

	  AAA_inv(0,0) = 1/ abx1;
	  AAA_inv(0,1) =  -1/ abx2;
	  AAA_inv(0,2) = (ix_temp1 + ix_temp2 - 2*ix_temp3)/ abx3;
	  AAA_inv(0,3) = 1/ abx4;
	  
	  AAA_inv(1,0) = -(ix_temp2 + 2*ix_temp3)/ abx1;
	  AAA_inv(1,1) = (ix_temp1 + 2*ix_temp3)/ abx2;
	  AAA_inv(1,2) = -(pow(ix_temp1, 2) + ix_temp1*ix_temp2 + pow(ix_temp2, 2) - 3*pow(ix_temp3, 2))/ abx3;
	  AAA_inv(1,3) = -(ix_temp1 + ix_temp2 + ix_temp3)/ abx4;
	  
	  AAA_inv(2,0) = (ix_temp3*(2*ix_temp2 + ix_temp3))/ abx1;
	  AAA_inv(2,1) = -(ix_temp3*(2*ix_temp1 + ix_temp3))/ abx2;
	  AAA_inv(2,2) = (ix_temp3*(2*pow(ix_temp1, 2) + 2*ix_temp1*ix_temp2 - 3*ix_temp3*ix_temp1 + 2*pow(ix_temp2, 2) - 3*ix_temp3*ix_temp2))/ abx3;
	  AAA_inv(2,3) = (ix_temp1*ix_temp2 + ix_temp1*ix_temp3 + ix_temp2*ix_temp3)/ abx4;
	  
	  AAA_inv(3,0) = -(ix_temp2*pow(ix_temp3, 2))/ abx1;
	  AAA_inv(3,1) = (ix_temp1*pow(ix_temp3, 2))/ abx2;
	  AAA_inv(3,2) = (ix_temp1*ix_temp2*(ix_temp1*ix_temp2 - 2*ix_temp1*ix_temp3 - 2*ix_temp2*ix_temp3 + 3*pow(ix_temp3, 2)))/ abx3;
	  AAA_inv(3,3) = -(ix_temp1*ix_temp2*ix_temp3)/ abx4;
	  	  
	  	   
	   
	   
	   
	   
	   _comy_matrix_inv = AAA_inv * _comy_temp;
	   

	   
	   _comy.col(i) = _t_temp* _comy_matrix_inv;
	   _comy.col(i+1) = _t_temp1* _comy_matrix_inv;	   
	   
	   _comx(0,i) = _comx(0,_n_end_walking-1);
	   _comx(0,i+1) = _comx(0,_n_end_walking-1);
	   _comz(0,i) = _comz(0,_n_end_walking-1);
	   _comz(0,i+1) = _comz(0,_n_end_walking-1);	   
	   _thetax(0,i) = _thetax(0,_n_end_walking-1);
	   _thetax(0,i+1) = _thetax(0,_n_end_walking-1);	   
	   _thetay(0,i) = _thetay(0,_n_end_walking-1);
	   _thetay(0,i+1) = _thetay(0,_n_end_walking-1);


	   
	  _footx_real_next.row(i+_nT -1)=_footx_real_next.row(_n_end_walking-1+_nT -1)	;
	  _footy_real_next.row(i+_nT -1)=_footy_real_next.row(_n_end_walking-1-_nT -1)	;	   
	   
	   for (int jxxx = _bjxx+1; jxxx<_footstepsnumber; jxxx++){
	   	    _footx_real(jxxx) = _footx_real(_bjxx) ;
	   	    _footy_real(jxxx) = _footy_real(_bjxx-2);		     
	     
	  }
	    _footxyz_real.row(0) = _footx_real.transpose();
	    _footxyz_real.row(1) = _footy_real.transpose();	  
	    _footxyz_real.row(2) = _footz_real.transpose();
  
	   
	}
	 else
	 {
	   
	   _comy(0,i) = _comy(0,_n_end_walking+round(_tstep/_dt/2));
	   _comy(0,i+1) = _comy(0,_n_end_walking+round(_tstep/_dt/2));
	   
	   _comx(0,i) = _comx(0,_n_end_walking-1);
	   _comx(0,i+1) = _comx(0,_n_end_walking-1);
	   _comz(0,i) = _comz(0,_n_end_walking-1);
	   _comz(0,i+1) = _comz(0,_n_end_walking-1);	   
	   _thetax(0,i) = _thetax(0,_n_end_walking-1);
	   _thetax(0,i+1) = _thetax(0,_n_end_walking-1);	   
	   _thetay(0,i) = _thetay(0,_n_end_walking-1);
	   _thetay(0,i+1) = _thetay(0,_n_end_walking-1);
	   
	  _footx_real_next.row(i+_nT -1)=_footx_real_next.row(_n_end_walking-1+_nT -1)	;
	  _footy_real_next.row(i+_nT -1)=_footy_real_next.row(_n_end_walking-1-_nT -1)	;	   
	   
	   for (int jxxx = _bjxx+1; jxxx<_footstepsnumber; jxxx++){
	   	    _footx_real(jxxx) = _footx_real(_bjxx) ;
	   	    _footy_real(jxxx) = _footy_real(_bjxx-2);		     
	     
	  }
	    _footxyz_real.row(0) = _footx_real.transpose();
	    _footxyz_real.row(1) = _footy_real.transpose();	  
	    _footxyz_real.row(2) = _footz_real.transpose();
   
	}
	 
       }
        
  


}



//////////////////////////// modified
int MPCClass::Indexfind(double goalvari, int xyz)
{
        _j_period = 0;
	if (xyz<0.05)
	{
	  while (goalvari >= _tx(_j_period))
	  {
	    _j_period++;
	  }
	  
	  _j_period = _j_period-1;	  
	}
	else
	{
	  while ( fabs(goalvari - _t_f(_j_period)) >0.0001 )
	  {
	    _j_period++;
	  }
	  	  
	}	  
	
// 	return j;
  
}



////////////////////////////modified
int MPCClass::is_sigular(int num)
{
	if ( num & 1)
	{	  
	  return 1; //sigular 
	}
	else
	{
	  return 0;	// odd  	 
	}	   
}

///// only walking once when initialize

Eigen::MatrixXd  MPCClass::Matrix_ps(Eigen::MatrixXd a, int nh,Eigen::MatrixXd cxps)
{
//   Eigen::MatrixXd matrixps(nh,3);
  Eigen::MatrixXd matrixps;
  matrixps.setZero(nh,3);  
  
  
  
  Eigen::MatrixXd A;
//   A.setIdentity(a.rows(),a.cols());
  
  for (int i = 0; i < nh; i++) {
    A.setIdentity(a.rows(),a.cols());
    for (int j = 1; j < i+2; j++)
    {
      A = A*a;
    }  
    
     matrixps.middleRows(i, 1)= cxps * A;      
  }
    
  return matrixps;
}


Eigen::MatrixXd MPCClass::Matrix_pu(Eigen::MatrixXd a, Eigen::MatrixXd b, int nh, Eigen::MatrixXd cxpu)
{
  Eigen::MatrixXd matrixpu;
  matrixpu.setZero(nh,nh);
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd Tempxx;
  
  
  for (int i = 1; i < nh+1; i++) {
    for (int j = 1; j < i+1; j++)
    { 
      A.setIdentity(a.rows(),a.cols());      
      if (j==i)
      {
	Tempxx = cxpu * A * b;
	matrixpu(i-1,j-1) = Tempxx(0,0);
      }
      else
      {	
	for (int k = 1; k < i-j+1; k++)
	{
	  A = A*a;
	}
	Tempxx = cxpu * A * b;
	matrixpu(i-1,j-1) = Tempxx(0,0);
      }          
    }       
  }
    
  return matrixpu;  
}



///DATA SAVING:modified=========================================================

void MPCClass::File_wl()
{
        
// 	CoMMM_ZMP_foot.setZero();
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(0, 0,1,_comx.cols()) = _comx;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(1, 0,1,_comx.cols()) = _comy;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(2, 0,1,_comx.cols()) = _comz;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(3, 0,1,_comx.cols()) = _zmpx_real;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(4, 0,1,_comx.cols()) = _zmpy_real;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(5, 0,1,_comx.cols()) = _thetax;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(6, 0,1,_comx.cols()) = _thetay;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(7, 0,1,_comx.cols()) = _torquex_real;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(8, 0,1,_comx.cols()) = _torquey_real;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(9, 0,1,_comx.cols()) = _footx_real_next1.transpose();	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(10, 0,1,_comx.cols()) = _footy_real_next1.transpose();	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(11, 0,1,_comx.cols()) = _footz_real_next.transpose();
	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(12, 0,1,_comx.cols()) = _Lfootx;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(13, 0,1,_comx.cols()) = _Lfooty;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(14, 0,1,_comx.cols()) = _Lfootz;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(15, 0,1,_comx.cols()) = _Rfootx;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(16, 0,1,_comx.cols()) = _Rfooty;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(17, 0,1,_comx.cols()) = _Rfootz;

	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(18, 0,1,_comx.cols()) = _comvx;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(19, 0,1,_comx.cols()) = _comax;
	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(20, 0,1,_comx.cols()) = _comvy;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(21, 0,1,_comx.cols()) = _comay;
	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(22, 0,1,_comx.cols()) = _comvz;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(23, 0,1,_comx.cols()) = _comaz;	

	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(24, 0,1,_comx.cols()) = _thetavx;
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(25, 0,1,_comx.cols()) = _thetaax;
	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(26, 0,1,_comx.cols()) = _thetavy;	
	CoMMM_ZMP_foot.block< Dynamic, Dynamic>(27, 0,1,_comx.cols()) = _thetaay;

	
	
	
	
  
	std::string fileName = "C++_NMPC2018_3robut3_runtime.txt" ;
	std::ofstream outfile( fileName.c_str() ) ; // file name and the operation type. 
       
        for(int i=0; i<_tcpu.rows(); i++){
           for(int j=0; j<_tcpu.cols(); j++){
                 outfile << (double) _tcpu(i,j) << " " ; 
           }
           outfile << std::endl;       // a   newline
        }
        outfile.close();
	
        for(int i=0; i<_tcpu.rows(); i++){
           for(int j=0; j<_tcpu.cols(); j++){
                 outfile << (double) _tcpu(i,j) << " " ; 
           }
           outfile << std::endl;       // a   newline
        }
        outfile.close();	


	std::string fileName1 = "C++_NMPC2018_3robut3_optimal_trajectory.txt" ;
	std::ofstream outfile1( fileName1.c_str() ) ; // file name and the operation type.        
	
        for(int i=0; i<CoMMM_ZMP_foot.rows(); i++){
           for(int j=0; j<CoMMM_ZMP_foot.cols(); j++){
                 outfile1 << (double) CoMMM_ZMP_foot(i,j) << " " ; 
           }
           outfile1 << std::endl;       // a   newline
        }
        outfile1.close();	
	
	
	
}


///// three model MPC solution :modified================================================================

void MPCClass::solve_reactive_step()
{
/*  int nVars = _Nt;
  int nEqCon = 1+3*_nh;
  int nIneqCon = 5*_nh + 4*_nstep+4;  
  resizeQP(nVars, nEqCon, nIneqCon);	*/    

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;
	    
  

  _CI.block< Dynamic, Dynamic>(0,0, _Nt,_nh) = _H_q_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,_nh, _Nt,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,2*_nh, _Nt,_nh) = _H_q_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,3*_nh, _Nt,_nh) = _H_q_lowy.transpose() * (-1); 
  _CI.block< Dynamic, Dynamic>(0,4*_nh, _Nt,_nh) = _H_hacc_lowz.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh, _Nt,_nstep) = _H_q_footx_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+_nstep, _Nt,_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+2*_nstep, _Nt,_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+3*_nstep, _Nt,_nstep) = _H_q_footy_low.transpose() * (-1);  

  
  _CI.block< Dynamic, Dynamic>(0,5*_nh+4*_nstep, _Nt,1) = _Footvx_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+4*_nstep+1, _Nt,1) = _Footvx_min.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+4*_nstep+2, _Nt,1) = _Footvy_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+4*_nstep+3, _Nt,1) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block< Dynamic, Dynamic>(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block< Dynamic, Dynamic>(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block< Dynamic, Dynamic>(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block< Dynamic, Dynamic>(3*_nh, 0,_nh,1) = _F_zmp_lowy;
  
  _ci0.block< Dynamic, Dynamic>(4*_nh, 0,_nh,1) = _F_hacc_lowz;
//   
  _ci0.block< Dynamic, Dynamic>(5*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block< Dynamic, Dynamic>(5*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block< Dynamic, Dynamic>(5*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block< Dynamic, Dynamic>(5*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;
//   
  _ci0.block< Dynamic, Dynamic>(5*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block< Dynamic, Dynamic>(5*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block< Dynamic, Dynamic>(5*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block< Dynamic, Dynamic>(5*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  _CE.block< Dynamic, Dynamic>(0,0, _Nt,1) = _H_q_footz.transpose();
  _CE.block< Dynamic, Dynamic>(0,1, _Nt,_nh) = _h_h.transpose();
  _CE.block< Dynamic, Dynamic>(0,_nh+1, _Nt,_nh) = _a_hx.transpose();
  _CE.block< Dynamic, Dynamic>(0,2*_nh+1, _Nt,_nh) = _a_hy.transpose();
  
  _ce0.block< Dynamic, Dynamic>(0,0, 1,1) = _F_footz;
  _ce0.block< Dynamic, Dynamic>(1,0, _nh,1) = _hhhx;  
  _ce0.block< Dynamic, Dynamic>(1+_nh,0, _nh,1) = _a_hxx;  
  _ce0.block< Dynamic, Dynamic>(1+2*_nh,0, _nh,1) = _a_hyy;    
  
  
  
  Solve();  

}


void MPCClass::solve_reactive_step_body_inclination()
{
/*  int nVars = _Nt;
  int nEqCon = 1+_nh;
  int nIneqCon = 13*_nh + 4*_nstep +4;
  resizeQP(nVars, nEqCon, nIneqCon);*/	    

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;

  
  _CI.block< Dynamic, Dynamic>(0,0, _Nt,_nh) = _H_q_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,_nh, _Nt,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,2*_nh, _Nt,_nh) = _H_q_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,3*_nh, _Nt,_nh) = _H_q_lowy.transpose() * (-1);
  
//   _CI.block< Dynamic, Dynamic>(0,4*_nh, _Nt,_nh) = _H_h_upz.transpose() * (-1);
//   _CI.block< Dynamic, Dynamic>(0,5*_nh, _Nt,_nh) = _H_h_lowz.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,4*_nh, _Nt,_nh) = _H_hacc_lowz.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,5*_nh, _Nt,_nstep) = _H_q_footx_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+_nstep, _Nt,_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+2*_nstep, _Nt,_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh+3*_nstep, _Nt,_nstep) = _H_q_footy_low.transpose() * (-1);	    

  _CI.block< Dynamic, Dynamic>(0,5*_nh+4*_nstep, _Nt,_nh) = _q_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,6*_nh+4*_nstep, _Nt,_nh) = _q_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,7*_nh+4*_nstep, _Nt,_nh) = _q_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,8*_nh+4*_nstep,_Nt,_nh) = _q_lowy.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,9*_nh+4*_nstep, _Nt,_nh) = _t_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,10*_nh+4*_nstep, _Nt,_nh) = _t_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,11*_nh+4*_nstep, _Nt,_nh) = _t_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,12*_nh+4*_nstep, _Nt,_nh) = _t_lowy.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,13*_nh+4*_nstep, _Nt,1) = _Footvx_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,13*_nh+4*_nstep+1, _Nt,1) = _Footvx_min.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,13*_nh+4*_nstep+2, _Nt,1) = _Footvy_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,13*_nh+4*_nstep+3, _Nt,1) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block< Dynamic, Dynamic>(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block< Dynamic, Dynamic>(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block< Dynamic, Dynamic>(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block< Dynamic, Dynamic>(3*_nh, 0,_nh,1) = _F_zmp_lowy;
  
//   _ci0.block< Dynamic, Dynamic>(4*_nh, 0,_nh,1) = _F_h_upz;
//   _ci0.block< Dynamic, Dynamic>(5*_nh, 0,_nh,1) = _F_h_lowz;
  
  _ci0.block< Dynamic, Dynamic>(4*_nh, 0,_nh,1) = _F_hacc_lowz;
  
  _ci0.block< Dynamic, Dynamic>(5*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block< Dynamic, Dynamic>(5*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block< Dynamic, Dynamic>(5*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block< Dynamic, Dynamic>(5*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;	    

  _ci0.block< Dynamic, Dynamic>(5*_nh+4*_nstep, 0,_nh,1) = _qq_upx;
  _ci0.block< Dynamic, Dynamic>(6*_nh+4*_nstep, 0,_nh,1) = _qq_lowx;
  _ci0.block< Dynamic, Dynamic>(7*_nh+4*_nstep, 0,_nh,1) = _qq_upy;
  _ci0.block< Dynamic, Dynamic>(8*_nh+4*_nstep, 0,_nh,1) = _qq_lowy;
  
  _ci0.block< Dynamic, Dynamic>(9*_nh+4*_nstep, 0,_nh,1) = _tt_upx;
  _ci0.block< Dynamic, Dynamic>(10*_nh+4*_nstep, 0,_nh,1) = _tt_lowx;
  _ci0.block< Dynamic, Dynamic>(11*_nh+4*_nstep, 0,_nh,1) = _tt_upy;
  _ci0.block< Dynamic, Dynamic>(12*_nh+4*_nstep, 0,_nh,1) = _tt_lowy;
  
  _ci0.block< Dynamic, Dynamic>(13*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block< Dynamic, Dynamic>(13*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block< Dynamic, Dynamic>(13*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block< Dynamic, Dynamic>(13*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  _CE.block< Dynamic, Dynamic>(0,0, _Nt,1) = _H_q_footz.transpose();
  _CE.block< Dynamic, Dynamic>(0,1, _Nt,_nh) = _h_h.transpose();
  
  _ce0.block< Dynamic, Dynamic>(0,0, 1,1) = _F_footz;
  _ce0.block< Dynamic, Dynamic>(1,0, _nh,1) = _hhhx;  

  Solve();  

}


void MPCClass::solve_reactive_step_body_inclination_CoMz()
{
/*  int nVars = _Nt;
  int nEqCon = 1;
  int nIneqCon = 15*_nh + 4*_nstep +4;
  resizeQP(nVars, nEqCon, nIneqCon);	*/    

  _G = _Q_goal1;
  _g0 = _q_goal1;
  _X = _V_ini;

  
  _CI.block< Dynamic, Dynamic>(0,0, _Nt,_nh) = _H_q_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,_nh, _Nt,_nh) = _H_q_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,2*_nh, _Nt,_nh) = _H_q_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,3*_nh, _Nt,_nh) = _H_q_lowy.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,4*_nh, _Nt,_nh) = _H_h_upz.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,5*_nh, _Nt,_nh) = _H_h_lowz.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,6*_nh, _Nt,_nh) = _H_hacc_lowz.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,7*_nh, _Nt,_nstep) = _H_q_footx_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,7*_nh+_nstep, _Nt,_nstep) = _H_q_footx_low.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,7*_nh+2*_nstep, _Nt,_nstep) = _H_q_footy_up.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,7*_nh+3*_nstep, _Nt,_nstep) = _H_q_footy_low.transpose() * (-1);	    

  _CI.block< Dynamic, Dynamic>(0,7*_nh+4*_nstep, _Nt,_nh) = _q_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,8*_nh+4*_nstep, _Nt,_nh) = _q_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,9*_nh+4*_nstep, _Nt,_nh) = _q_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,10*_nh+4*_nstep,_Nt,_nh) = _q_lowy.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,11*_nh+4*_nstep, _Nt,_nh) = _t_upx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,12*_nh+4*_nstep, _Nt,_nh) = _t_lowx.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,13*_nh+4*_nstep, _Nt,_nh) = _t_upy.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,14*_nh+4*_nstep, _Nt,_nh) = _t_lowy.transpose() * (-1);
  
  _CI.block< Dynamic, Dynamic>(0,15*_nh+4*_nstep, _Nt,1) = _Footvx_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,15*_nh+4*_nstep+1, _Nt,1) = _Footvx_min.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,15*_nh+4*_nstep+2, _Nt,1) = _Footvy_max.transpose() * (-1);
  _CI.block< Dynamic, Dynamic>(0,15*_nh+4*_nstep+3, _Nt,1) = _Footvy_min.transpose() * (-1);

  

  
  _ci0.block< Dynamic, Dynamic>(0, 0,_nh,1) = _F_zmp_upx;
  _ci0.block< Dynamic, Dynamic>(_nh, 0,_nh,1) = _F_zmp_lowx;
  _ci0.block< Dynamic, Dynamic>(2*_nh, 0,_nh,1) = _F_zmp_upy;
  _ci0.block< Dynamic, Dynamic>(3*_nh, 0,_nh,1) = _F_zmp_lowy;
  
  _ci0.block< Dynamic, Dynamic>(4*_nh, 0,_nh,1) = _F_h_upz;
  _ci0.block< Dynamic, Dynamic>(5*_nh, 0,_nh,1) = _F_h_lowz;
  
  _ci0.block< Dynamic, Dynamic>(6*_nh, 0,_nh,1) = _F_hacc_lowz;
  
  _ci0.block< Dynamic, Dynamic>(7*_nh, 0,_nstep,1) = _F_foot_upx;
  _ci0.block< Dynamic, Dynamic>(7*_nh+_nstep, 0,_nstep,1) = _F_foot_lowx;
  _ci0.block< Dynamic, Dynamic>(7*_nh+2*_nstep, 0,_nstep,1) = _F_foot_upy;
  _ci0.block< Dynamic, Dynamic>(7*_nh+3*_nstep, 0,_nstep,1) = _F_foot_lowy;	    

  _ci0.block< Dynamic, Dynamic>(7*_nh+4*_nstep, 0,_nh,1) = _qq_upx;
  _ci0.block< Dynamic, Dynamic>(8*_nh+4*_nstep, 0,_nh,1) = _qq_lowx;
  _ci0.block< Dynamic, Dynamic>(9*_nh+4*_nstep, 0,_nh,1) = _qq_upy;
  _ci0.block< Dynamic, Dynamic>(10*_nh+4*_nstep, 0,_nh,1) = _qq_lowy;
  
  _ci0.block< Dynamic, Dynamic>(11*_nh+4*_nstep, 0,_nh,1) = _tt_upx;
  _ci0.block< Dynamic, Dynamic>(12*_nh+4*_nstep, 0,_nh,1) = _tt_lowx;
  _ci0.block< Dynamic, Dynamic>(13*_nh+4*_nstep, 0,_nh,1) = _tt_upy;
  _ci0.block< Dynamic, Dynamic>(14*_nh+4*_nstep, 0,_nh,1) = _tt_lowy;
  
  _ci0.block< Dynamic, Dynamic>(15*_nh+4*_nstep, 0,1,1) = _footubxv;
  _ci0.block< Dynamic, Dynamic>(15*_nh+4*_nstep+1, 0,1,1) = _footlbxv;
  _ci0.block< Dynamic, Dynamic>(15*_nh+4*_nstep+2, 0,1,1) = _footubyv;
  _ci0.block< Dynamic, Dynamic>(15*_nh+4*_nstep+3, 0,1,1) = _footlbyv;

  
  
  
  
  
  _CE = _H_q_footz.transpose();
  _ce0 = _F_footz;
  
  
  
  Solve();  

}



void MPCClass::Solve()
{
// min 0.5 * x G x + g0 x
// _s.t.
// 		CE^T x + ce0 = 0
// 		CI^T x + ci0 >= 0
		solveQP();
		if (_X.rows() == _Nt)
		{
		  _V_ini += _X;
		}

}









///////////////////////////////////////////////////////////////////////////////////=================swing foot trajectory==============================================////////////////////////
/////////////============================================================================================================================////////////////////////////////////////////////


//// foot trajectory solve--------polynomial ================================================
void MPCClass::Foot_trajectory_solve(int j_index,bool _stopwalking)
{
  // maximal swing foot height: 
//   double  Footz_ref = _lift_height;
  
  //////////////external push mode==========
//   double  Footz_ref = 0.1;  
  ///// stairs :============================
//    double  Footz_ref = 0.07;  //0.02m
//     double  Footz_ref = 0.09;  //0.05m  
//    double  Footz_ref = 0.1;      //0.06m 
 
//    double  Footz_ref = 0.12;      //0.07m 
//    double  Footz_ref = 0.13;      //0.08m,0.09m 
//    double  Footz_ref = 0.14;      //0.1m,

//// judge if stop  
        if(_stopwalking)  
	{
	  
	  for (int i_t = _bjx1+1; i_t < _footstepsnumber; i_t++) {	  
	    _lift_height_ref(i_t) = 0;  
	  }	  

	}  
  
  
  
  _footxyz_real(1,0) = -_stepwidth(0);
  
//   foot trajectory generation:
  if (_bjx1 >= 2)
  {
//     cout << "_bjx1 >= 2"<<endl;
    if (_bjx1 % 2 == 0)           //odd:left support
    {
//     no change on the left support location
      _Lfootx(j_index) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Lfooty(j_index) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Lfootz(j_index) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);
      
      _Lfootx(j_index+1) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Lfooty(j_index+1) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Lfootz(j_index+1) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);    
  
      /// right swing
      if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double support
      {
// 	cout << "dsp"<<endl;
	_Rfootx(j_index) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);	
	
	_Rfootx(j_index+1) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfooty(j_index+1) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Rfootz(j_index+1) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);	
	
	
      }
      else
      {
	
// 	cout << "ssp"<<endl;
	//initial state and final state and the middle state
	double t_des = (j_index +1 - round(_tx(_bjx1-1)/_dt) +1)*_dt;
	Eigen::Vector3d t_plan;
	t_plan(0) = t_des - _dt;
	t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
	t_plan(2) = _ts(_bjx1-1) + 0.0001;
	
	if (abs(t_des - _ts(_bjx1-1)) <= (_dt + 0.0005))
	{
	  _Rfootx(j_index) = _footxyz_real(0,_bjxx); 
	  _Rfooty(j_index) = _footxyz_real(1,_bjxx);
	  _Rfootz(j_index) = _footxyz_real(2,_bjxx); 	
	  
	  _Rfootx(j_index+1) = _footxyz_real(0,_bjxx); 
	  _Rfooty(j_index+1) = _footxyz_real(1,_bjxx);
	  _Rfootz(j_index+1) = _footxyz_real(2,_bjxx); 	  
	  
	}
	else
	{
	  Eigen::Matrix<double,7,7> AAA;
	  AAA.setZero();
	  Eigen::Matrix<double, 1, 7> aaaa;
	  aaaa.setZero();
	  
	  aaaa(0) = 6*pow(t_plan(0), 5);   aaaa(1) =  5*pow(t_plan(0), 4);  aaaa(2) =  4*pow(t_plan(0), 3);   aaaa(3) =  3*pow(t_plan(0), 2);
	  aaaa(4) = 2*pow(t_plan(0), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
	  AAA.row(0) = aaaa;
	  
	  aaaa(0) = 30*pow(t_plan(0), 4);  aaaa(1) =  20*pow(t_plan(0), 3);  aaaa(2) =  12*pow(t_plan(0), 2);   aaaa(3) =  6*pow(t_plan(0), 1);
	  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
	  AAA.row(1) = aaaa;
	  
	  aaaa(0) = pow(t_plan(0), 6);     aaaa(1) =  pow(t_plan(0), 5);     aaaa(2) =  pow(t_plan(0), 4);   aaaa(3) =  pow(t_plan(0), 3);
	  aaaa(4) = pow(t_plan(0), 2);     aaaa(5) = pow(t_plan(0), 1);      aaaa(6) =  1;	  
	  AAA.row(2) = aaaa;
	  
	  aaaa(0) = pow(t_plan(1), 6);     aaaa(1) =  pow(t_plan(1), 5);     aaaa(2) =  pow(t_plan(1), 4);   aaaa(3) =  pow(t_plan(1), 3);
	  aaaa(4) = pow(t_plan(1), 2);     aaaa(5) = pow(t_plan(1), 1);      aaaa(6) =  1;
	  AAA.row(3) = aaaa;
	  
	  aaaa(0) = pow(t_plan(2), 6);     aaaa(1) =  pow(t_plan(2), 5);     aaaa(2) =  pow(t_plan(2), 4);   aaaa(3) =  pow(t_plan(2), 3);
	  aaaa(4) = pow(t_plan(2), 2);     aaaa(5) = pow(t_plan(2), 1);      aaaa(6) =  1;
	  AAA.row(4) = aaaa;	  
	  
	  aaaa(0) = 6*pow(t_plan(2), 5);   aaaa(1) =  5*pow(t_plan(2), 4);  aaaa(2) =  4*pow(t_plan(2), 3);   aaaa(3) =  3*pow(t_plan(2), 2);
	  aaaa(4) = 2*pow(t_plan(2), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
	  AAA.row(5) = aaaa;
	  
	  aaaa(0) = 30*pow(t_plan(2), 4);  aaaa(1) =  20*pow(t_plan(2), 3);  aaaa(2) =  12*pow(t_plan(2), 2);   aaaa(3) =  6*pow(t_plan(2), 1);
	  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
	  AAA.row(6) = aaaa;
	  
	  Eigen::Matrix<double,7,7> AAA_inv;
// 	  AAA_inv.setZero(); 
	  AAA_inv = AAA.inverse();
	 	  
	  Eigen::Matrix<double, 1, 7> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_des, 6);   t_a_plan(1) = pow(t_des, 5);   t_a_plan(2) = pow(t_des, 4);  t_a_plan(3) = pow(t_des, 3);
	  t_a_plan(4) = pow(t_des, 2);   t_a_plan(5) = pow(t_des, 1);   t_a_plan(6) = 1;
	  

	  Eigen::Matrix<double, 1, 7> t_a_planv;
	  t_a_planv.setZero();
	  t_a_planv(0) = 6*pow(t_des, 5);   t_a_planv(1) = 5*pow(t_des, 4);   t_a_planv(2) = 4*pow(t_des, 3);  t_a_planv(3) = 3*pow(t_des, 2);
	  t_a_planv(4) = 2*pow(t_des, 1);   t_a_planv(5) = 1;                 t_a_planv(6) = 0;
	  
	  
	  Eigen::Matrix<double, 1, 7> t_a_plana;
	  t_a_plana.setZero();
	  t_a_plana(0) = 30*pow(t_des, 4);   t_a_plana(1) = 20*pow(t_des, 3);   t_a_plana(2) = 12*pow(t_des, 2);  t_a_plana(3) = 6*pow(t_des, 1);
	  t_a_plana(4) = 2;                  t_a_plana(5) = 0;                  t_a_plana(6) = 0;
	  
// 	  cout<<"AAA="<<endl<<AAA<<endl;
// 	  cout <<"AAA_inverse="<<endl<<AAA.inverse()<<endl;	  
// 	  cout <<"t_des="<<endl<<t_des<<endl;
// 	  cout <<"t_plan="<<endl<<t_plan<<endl;
	  
	  ////////////////////////////////////////////////////////////////////////////
	  Eigen::Matrix<double, 7, 1> Rfootx_plan;
	  Rfootx_plan.setZero();	
	  Rfootx_plan(0) = _Rfootvx(j_index-1);     Rfootx_plan(1) = _Rfootax(j_index-1); Rfootx_plan(2) = _Rfootx(j_index-1); Rfootx_plan(3) = _Lfootx(j_index);
	  Rfootx_plan(4) = _footxyz_real(0,_bjxx);  Rfootx_plan(5) = 0;                   Rfootx_plan(6) = 0;
	  
	  
	  Eigen::Matrix<double, 7, 1> Rfootx_co;
	  Rfootx_co.setZero();
	  Rfootx_co = AAA_inv * Rfootx_plan;
	  
	  _Rfootx(j_index) = t_a_plan * Rfootx_co;
	  _Rfootvx(j_index) = t_a_planv * Rfootx_co;
	  _Rfootax(j_index) = t_a_plana * Rfootx_co;
	  
	  /////////////////////////////////////////////////////////////////////////////
	  if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
	  {
	    _ry_left_right = (_footxyz_real(1,_bjxx) + _footxyz_real(1,_bjxx-2))/2;
	  }
	  
	  Eigen::Matrix<double, 7, 1> Rfooty_plan;
	  Rfooty_plan.setZero();	
	  Rfooty_plan(0) = _Rfootvy(j_index-1);     Rfooty_plan(1) = _Rfootay(j_index-1); Rfooty_plan(2) = _Rfooty(j_index-1); Rfooty_plan(3) = _ry_left_right;
	  Rfooty_plan(4) = _footxyz_real(1,_bjxx);  Rfooty_plan(5) = 0;                   Rfooty_plan(6) = 0;	    
	  
	  Eigen::Matrix<double, 7, 1> Rfooty_co;
	  Rfooty_co.setZero();
	  Rfooty_co = AAA_inv * Rfooty_plan;
	  
	  _Rfooty(j_index) = t_a_plan * Rfooty_co;
	  _Rfootvy(j_index) = t_a_planv * Rfooty_co;
	  _Rfootay(j_index) = t_a_plana * Rfooty_co;	
	  
	  
	  //////////////////////////////////////////////////////////
	  Eigen::Matrix<double, 7, 1> Rfootz_plan;
	  Rfootz_plan.setZero();	
	  Rfootz_plan(0) = _Rfootvz(j_index-1);     Rfootz_plan(1) = _Rfootaz(j_index-1); Rfootz_plan(2) = _Rfootz(j_index-1); Rfootz_plan(3) = _Lfootz(j_index)+_lift_height_ref(_bjx1-1);
	  Rfootz_plan(4) = _footxyz_real(2,_bjxx);  Rfootz_plan(5) = 0;                   Rfootz_plan(6) = 0;	
	  
	  Eigen::Matrix<double, 7, 1> Rfootz_co;
	  Rfootz_co.setZero();
	  Rfootz_co = AAA_inv * Rfootz_plan;
	  
	  _Rfootz(j_index) = t_a_plan * Rfootz_co;
	  _Rfootvz(j_index) = t_a_planv * Rfootz_co;
	  _Rfootaz(j_index) = t_a_plana * Rfootz_co;	
	 
	  
	  _Rfootx(j_index+1) = _Rfootx(j_index)+_dt * _Rfootvx(j_index);
	  _Rfooty(j_index+1) = _Rfooty(j_index)+_dt * _Rfootvy(j_index);
	  _Rfootz(j_index+1) = _Rfootz(j_index)+_dt * _Rfootvz(j_index);
	  
	  
	  
	  //////spline for Rfoot_Z
/*	  std::vector<double> TXRFOOTZ(3), XRFOOTZ(3);
	  TXRFOOTZ[0] = t_plan(0); TXRFOOTZ[1] = t_plan(1); TXRFOOTZ[2] = t_plan(2);	  
	  XRFOOTZ[0] = Rfootx_plan(2);     XRFOOTZ[1] = Rfootx_plan(3);     XRFOOTZ[2] = Rfootx_plan(4);
	  
	  spline s1;
	  s1.set_points(TXRFOOTZ, XRFOOTZ);
	  _Rfootx(j_index)= s1(t_des);	  

	  
	  XRFOOTZ[0] = Rfooty_plan(2);     XRFOOTZ[1] = Rfooty_plan(3);     XRFOOTZ[2] = Rfooty_plan(4);	  
	  s1.set_points(TXRFOOTZ, XRFOOTZ);
	  _Rfooty(j_index) = s1(t_des);	  
	  	  
	  XRFOOTZ[0] = Rfootz_plan(2);     XRFOOTZ[1] = Rfootz_plan(3);     XRFOOTZ[2] = Rfootz_plan(4);	  
	  s1.set_points(TXRFOOTZ, XRFOOTZ);
	  _Rfootz(j_index) = s1(t_des);	*/	  
	  
	  
	}
      }
 
      
    }
    
    
    else                       //right support
    {
//       no change on right support
      _Rfootx(j_index) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfooty(j_index) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfootz(j_index) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);
      
      _Rfootx(j_index+1) = _Rfootx(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfooty(j_index+1) = _Rfooty(round(_tx(_bjx1-1)/_dt) -1-1);
      _Rfootz(j_index+1) = _Rfootz(round(_tx(_bjx1-1)/_dt) -1-1);      
  
      /// left swing
      if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1))  // j_index and _bjx1 coincident with matlab: double suppot
      {
// 	cout << "dsp"<<endl;
	_Lfootx(j_index) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);

	_Lfootx(j_index+1) = _Lfootx(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfooty(j_index+1) = _Lfooty(round(_tx(_bjx1-1)/_dt) -1-1);
	_Lfootz(j_index+1) = _Lfootz(round(_tx(_bjx1-1)/_dt) -1-1);
	
      }
      else
      {
// 	cout << "ssp"<<endl;
	//initial state and final state and the middle state
	double t_des = (j_index +1 - round(_tx(_bjx1-1)/_dt) +1)*_dt;
	Eigen::Vector3d t_plan;
	t_plan(0) = t_des - _dt;
	t_plan(1) = (_td(_bjx1-1) + _ts(_bjx1-1))/2 + 0.0001;
	t_plan(2) = _ts(_bjx1-1) + 0.0001;
	
	if (abs(t_des - _ts(_bjx1-1)) <= (_dt + 0.0005))
	{
	  
	  _Lfootx(j_index) = _footxyz_real(0,_bjxx); 
	  _Lfooty(j_index) = _footxyz_real(1,_bjxx);
	  _Lfootz(j_index) = _footxyz_real(2,_bjxx); 

	  _Lfootx(j_index+1) = _footxyz_real(0,_bjxx); 
	  _Lfooty(j_index+1) = _footxyz_real(1,_bjxx);
	  _Lfootz(j_index+1) = _footxyz_real(2,_bjxx); 
	  
	}
	else
	{
	  Eigen::Matrix<double, 7, 7> AAA;
	  AAA.setZero();
	  Eigen::Matrix<double, 1, 7> aaaa;
	  aaaa.setZero();
	  
	  aaaa(0) = 6*pow(t_plan(0), 5);   aaaa(1) =  5*pow(t_plan(0), 4);  aaaa(2) =  4*pow(t_plan(0), 3);   aaaa(3) =  3*pow(t_plan(0), 2);
	  aaaa(4) = 2*pow(t_plan(0), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
	  AAA.row(0) = aaaa;
	  
	  aaaa(0) = 30*pow(t_plan(0), 4);  aaaa(1) =  20*pow(t_plan(0), 3);  aaaa(2) =  12*pow(t_plan(0), 2);   aaaa(3) =  6*pow(t_plan(0), 1);
	  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
	  AAA.row(1) = aaaa;
	  
	  aaaa(0) = pow(t_plan(0), 6);     aaaa(1) =  pow(t_plan(0), 5);     aaaa(2) =  pow(t_plan(0), 4);   aaaa(3) =  pow(t_plan(0), 3);
	  aaaa(4) = pow(t_plan(0), 2);     aaaa(5) = pow(t_plan(0), 1);      aaaa(6) =  1;	  
	  AAA.row(2) = aaaa;
	  
	  aaaa(0) = pow(t_plan(1), 6);     aaaa(1) =  pow(t_plan(1), 5);     aaaa(2) =  pow(t_plan(1), 4);   aaaa(3) =  pow(t_plan(1), 3);
	  aaaa(4) = pow(t_plan(1), 2);     aaaa(5) = pow(t_plan(1), 1);      aaaa(6) =  1;
	  AAA.row(3) = aaaa;
	  
	  aaaa(0) = pow(t_plan(2), 6);     aaaa(1) =  pow(t_plan(2), 5);     aaaa(2) =  pow(t_plan(2), 4);   aaaa(3) =  pow(t_plan(2), 3);
	  aaaa(4) = pow(t_plan(2), 2);     aaaa(5) = pow(t_plan(2), 1);      aaaa(6) =  1;
	  AAA.row(4) = aaaa;	  
	  
	  aaaa(0) = 6*pow(t_plan(2), 5);   aaaa(1) =  5*pow(t_plan(2), 4);  aaaa(2) =  4*pow(t_plan(2), 3);   aaaa(3) =  3*pow(t_plan(2), 2);
	  aaaa(4) = 2*pow(t_plan(2), 1);   aaaa(5) = 1;                     aaaa(6) =  0;
	  AAA.row(5) = aaaa;
	  
	  aaaa(0) = 30*pow(t_plan(2), 4);  aaaa(1) =  20*pow(t_plan(2), 3);  aaaa(2) =  12*pow(t_plan(2), 2);   aaaa(3) =  6*pow(t_plan(2), 1);
	  aaaa(4) = 2;                     aaaa(5) = 0;                      aaaa(6) =  0;
	  AAA.row(6) = aaaa;
	  
          
	  Eigen::Matrix<double,7,7> AAA_inv;
// 	  AAA_inv.setZero(); 
	  AAA_inv = AAA.inverse();
	  
	  Eigen::Matrix<double, 1, 7> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_des, 6);   t_a_plan(1) = pow(t_des, 5);   t_a_plan(2) = pow(t_des, 4);  t_a_plan(3) = pow(t_des, 3);
	  t_a_plan(4) = pow(t_des, 2);   t_a_plan(5) = pow(t_des, 1);   t_a_plan(6) = 1;
	  

	  Eigen::Matrix<double, 1, 7> t_a_planv;
	  t_a_planv.setZero();
	  t_a_planv(0) = 6*pow(t_des, 5);   t_a_planv(1) = 5*pow(t_des, 4);   t_a_planv(2) = 4*pow(t_des, 3);  t_a_planv(3) = 3*pow(t_des, 2);
	  t_a_planv(4) = 2*pow(t_des, 1);   t_a_planv(5) = 1;                 t_a_planv(6) = 0;
	  
	  
	  Eigen::Matrix<double, 1, 7> t_a_plana;
	  t_a_plana.setZero();
	  t_a_plana(0) = 30*pow(t_des, 4);   t_a_plana(1) = 20*pow(t_des, 3);   t_a_plana(2) = 12*pow(t_des, 2);  t_a_plana(3) = 6*pow(t_des, 1);
	  t_a_plana(4) = 2;                  t_a_plana(5) = 0;                  t_a_plana(6) = 0;
	  
	  
	  ////////////////////////////////////////////////////////////////////////////
	  Eigen::Matrix<double, 7, 1> Lfootx_plan;
	  Lfootx_plan.setZero();	
	  Lfootx_plan(0) = _Lfootvx(j_index-1);     Lfootx_plan(1) = _Lfootax(j_index-1); Lfootx_plan(2) = _Lfootx(j_index-1); Lfootx_plan(3) = _Rfootx(j_index);
	  Lfootx_plan(4) = _footxyz_real(0,_bjxx);  Lfootx_plan(5) = 0;                   Lfootx_plan(6) = 0;	  
	  
	  
	  Eigen::Matrix<double, 7, 1> Lfootx_co;
	  Lfootx_co.setZero();
	  Lfootx_co = AAA_inv * Lfootx_plan;
	  
	  _Lfootx(j_index) = t_a_plan * Lfootx_co;
	  _Lfootvx(j_index) = t_a_planv * Lfootx_co;
	  _Lfootax(j_index) = t_a_plana * Lfootx_co;
	  
	  /////////////////////////////////////////////////////////////////////////////
	  if ((j_index +1 - round(_tx(_bjx1-1)/_dt))*_dt < _td(_bjx1-1)+_dt)
	  {
	    _ry_left_right = (_footxyz_real(1,_bjxx) + _footxyz_real(1,_bjxx-2))/2;
	  }
	  
	  Eigen::Matrix<double, 7, 1> Lfooty_plan;
	  Lfooty_plan.setZero();	
	  Lfooty_plan(0) = _Lfootvy(j_index-1);     Lfooty_plan(1) = _Lfootay(j_index-1); Lfooty_plan(2) = _Lfooty(j_index-1); Lfooty_plan(3) = _ry_left_right;
	  Lfooty_plan(4) = _footxyz_real(1,_bjxx);  Lfooty_plan(5) = 0;                   Lfooty_plan(6) = 0;		  
	  
	  
	  Eigen::Matrix<double, 7, 1> Lfooty_co;
	  Lfooty_co.setZero();
	  Lfooty_co = AAA_inv * Lfooty_plan;
	  
	  _Lfooty(j_index) = t_a_plan * Lfooty_co;
	  _Lfootvy(j_index) = t_a_planv * Lfooty_co;
	  _Lfootay(j_index) = t_a_plana * Lfooty_co;	
	  
	  
	  //////////////////////////////////////////////////////////
	  Eigen::Matrix<double, 7, 1> Lfootz_plan;
	  Lfootz_plan.setZero();		
	  Lfootz_plan(0) = _Lfootvz(j_index-1);     Lfootz_plan(1) = _Lfootaz(j_index-1); Lfootz_plan(2) = _Lfootz(j_index-1); Lfootz_plan(3) = _Rfootz(j_index)+_lift_height_ref(_bjx1-1);
	  Lfootz_plan(4) = _footxyz_real(2,_bjxx);  Lfootz_plan(5) = 0;                   Lfootz_plan(6) = 0;		  
	  
	  
	  Eigen::Matrix<double, 7, 1> Lfootz_co;
	  Lfootz_co.setZero();
	  Lfootz_co = AAA_inv * Lfootz_plan;
	  
	  _Lfootz(j_index) = t_a_plan * Lfootz_co;
	  _Lfootvz(j_index) = t_a_planv * Lfootz_co;
	  _Lfootaz(j_index) = t_a_plana * Lfootz_co;
	  
	  
	  _Lfootx(j_index+1) = _Lfootx(j_index)+_dt * _Lfootvx(j_index);
	  _Lfooty(j_index+1) = _Lfooty(j_index)+_dt * _Lfootvy(j_index);
	  _Lfootz(j_index+1) = _Lfootz(j_index)+_dt * _Lfootvz(j_index);
	  
	  
/*	  //////spline for Rfoot_Z
	  std::vector<double> TXRFOOTZ(3), XRFOOTZ(3);
	  TXRFOOTZ[0] = t_plan(0); TXRFOOTZ[1] = t_plan(1); TXRFOOTZ[2] = t_plan(2);	  
	  XRFOOTZ[0] = Lfootz_plan(2);     XRFOOTZ[1] = Lfootz_plan(3);     XRFOOTZ[2] = Lfootz_plan(4);
	  
	  spline s;
	  s.set_points(TXRFOOTZ, XRFOOTZ);
	  _Lfootz(j_index)= s(t_des);	  

	  
	  XRFOOTZ[0] = Lfooty_plan(2);     XRFOOTZ[1] = Lfooty_plan(3);     XRFOOTZ[2] = Lfooty_plan(4);	  
	  s.set_points(TXRFOOTZ, XRFOOTZ);
	  _Lfooty(j_index) = s(t_des);	  
	  	  
	  XRFOOTZ[0] = Lfootx_plan(2);     XRFOOTZ[1] = Lfootx_plan(3);     XRFOOTZ[2] = Lfootx_plan(4);	  
	  s.set_points(TXRFOOTZ, XRFOOTZ);
	  _Lfootx(j_index) = s(t_des);	*/  
	  
	  
	}
      }

    }
      
  }
  else
  {
    _Rfooty(j_index) = -_stepwidth(0);
    _Lfooty(j_index) = _stepwidth(0);
  }
    
// 	cout << "Rfooty_generated:"<<_Rfooty(j_index)<<endl;
  
}








///////////////////////// ODE - sampling time maximal ========================================
int MPCClass::Get_maximal_number_reference()
{
  return (_nsum -_nh-1);
}



int MPCClass::Get_maximal_number(double dtx)
{
  
  return (_nsum -_nh-1)*floor(_dt/dtx);
}






////====================================================================================================================
/////////////////////////// using the lower-level control-loop  sampling time as the reference: every 5ms;  at the same time: just using the next one position + next one velocity

Vector3d MPCClass::XGetSolution_CoM_position(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
  //reference com position
        _CoM_position_optimal.row(0) = _comx;
	_CoM_position_optimal.row(1) = _comy;
	_CoM_position_optimal.row(2) = _comz;
	_comz(0) = RobotParaClass::Z_C();
	_comz(1) = RobotParaClass::Z_C();
	_comz(2) = RobotParaClass::Z_C();
	_comz(3) = RobotParaClass::Z_C();
	_comz(4) = RobotParaClass::Z_C();
	
	
	
	Vector3d com_inte;	
	
	if (walktime>=2)
	{
	  int t_int; 
	  t_int = floor(walktime / (_dt / dt_sample) );

	  ///// chage to be relative time
	  double t_cur;
	  t_cur = walktime * dt_sample ;
	  

	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);

// 	  Eigen::MatrixXd AAA1;	
// 
// 	  AAA1.setZero(4,4);	
// 	  AAA1(0,0) = pow(t_plan(0), 3); AAA1(0,1) = pow(t_plan(0), 2); AAA1(0,2) = pow(t_plan(0), 1); AAA1(0,3) = pow(t_plan(0), 0); 
// 	  AAA1(1,0) = pow(t_plan(1), 3); AAA1(1,1) = pow(t_plan(1), 2); AAA1(1,2) = pow(t_plan(1), 1); AAA1(1,3) = pow(t_plan(0), 0); 
// 	  AAA1(2,0) = pow(t_plan(2), 3); AAA1(2,1) = pow(t_plan(2), 2); AAA1(2,2) = pow(t_plan(2), 1); AAA1(2,3) = pow(t_plan(0), 0); 
// 	  AAA1(3,0) = 3*pow(t_plan(2), 2); AAA1(3,1) = 2*pow(t_plan(2), 1); AAA1(3,2) = pow(t_plan(2), 0); AAA1(3,3) = 0;  


	  Eigen::Matrix4d AAA_inv;
	  
	  double abx1, abx2, abx3, abx4;
	  abx1 = ((t_plan(0) - t_plan(1))*pow(t_plan(0) - t_plan(2), 2));
	  abx2 = ((t_plan(0) - t_plan(1))*pow(t_plan(1) - t_plan(2), 2));
	  abx3 =(pow(t_plan(0) - t_plan(2), 2)*pow(t_plan(1) - t_plan(2), 2));
	  abx4 = ((t_plan(0) - t_plan(2))*(t_plan(1) - t_plan(2)));
	  

	  AAA_inv(0,0) = 1/ abx1;
	  AAA_inv(0,1) =  -1/ abx2;
	  AAA_inv(0,2) = (t_plan(0) + t_plan(1) - 2*t_plan(2))/ abx3;
	  AAA_inv(0,3) = 1/ abx4;
	  
	  AAA_inv(1,0) = -(t_plan(1) + 2*t_plan(2))/ abx1;
	  AAA_inv(1,1) = (t_plan(0) + 2*t_plan(2))/ abx2;
	  AAA_inv(1,2) = -(pow(t_plan(0), 2) + t_plan(0)*t_plan(1) + pow(t_plan(1), 2) - 3*pow(t_plan(2), 2))/ abx3;
	  AAA_inv(1,3) = -(t_plan(0) + t_plan(1) + t_plan(2))/ abx4;
	  
	  AAA_inv(2,0) = (t_plan(2)*(2*t_plan(1) + t_plan(2)))/ abx1;
	  AAA_inv(2,1) = -(t_plan(2)*(2*t_plan(0) + t_plan(2)))/ abx2;
	  AAA_inv(2,2) = (t_plan(2)*(2*pow(t_plan(0), 2) + 2*t_plan(0)*t_plan(1) - 3*t_plan(2)*t_plan(0) + 2*pow(t_plan(1), 2) - 3*t_plan(2)*t_plan(1)))/ abx3;
	  AAA_inv(2,3) = (t_plan(0)*t_plan(1) + t_plan(0)*t_plan(2) + t_plan(1)*t_plan(2))/ abx4;
	  
	  AAA_inv(3,0) = -(t_plan(1)*pow(t_plan(2), 2))/ abx1;
	  AAA_inv(3,1) = (t_plan(0)*pow(t_plan(2), 2))/ abx2;
	  AAA_inv(3,2) = (t_plan(0)*t_plan(1)*(t_plan(0)*t_plan(1) - 2*t_plan(0)*t_plan(2) - 2*t_plan(1)*t_plan(2) + 3*pow(t_plan(2), 2)))/ abx3;
	  AAA_inv(3,3) = -(t_plan(0)*t_plan(1)*t_plan(2))/ abx4;
	  	  
	  	  
	  
	  

	
	  
	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);   t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);   t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1);  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  // COM&&foot trajectory interpolation
	  
	  Eigen::Vector3d  x10;
	  Eigen::Vector3d  x11;
	  Eigen::Vector3d  x12;
// 	  Eigen::Vector3d  x13;


/*	  x10(0) = COM_IN(0,walktime-2); x10(1) = COM_IN(1,walktime-2); x10(2) = COM_IN(2,walktime-2);
	  x11(0) = COM_IN(0,walktime-1); x11(1) = COM_IN(1,walktime-1); x11(2) = COM_IN(2,walktime-1);	 */ 
	  x10 = body_in1; 
	  x11 = body_in2;  
	  x12 = _CoM_position_optimal.col(t_int);
// 	  x13 = _CoM_position_optimal.col(t_int+1);
	  
	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = x10(0); temp(1) = x11(0); temp(2) = x12(0); temp(3) = _comvx(t_int);	  
	  com_inte(0) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(1); temp(1) = x11(1); temp(2) = x12(1); temp(3) = _comvy(t_int);	  
	  com_inte(1) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(2); temp(1) = x11(2); temp(2) = x12(2); temp(3) = _comvz(t_int);	  
	  com_inte(2) = t_a_plan *(AAA_inv)*temp;

  
	  
	  
	}
	else
	{
	  com_inte(0) = body_in3(0);	  
	  com_inte(1) = body_in3(1);	  	  
	  com_inte(2) = body_in3(2);	
	  
	}

 	return com_inte;
	
}


Vector3d MPCClass::XGetSolution_body_inclination(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
  //reference com position
        _torso_angle_optimal.row(0) = _thetax;
	_torso_angle_optimal.row(1) = _thetay;
	_torso_angle_optimal.row(2) = _thetaz;
	
	
	Vector3d com_inte;	
	
	if (walktime>=2)
	{
	  int t_int; 
	  t_int = floor(walktime / (_dt / dt_sample) );

	  
	  double t_cur;
	  t_cur = walktime * dt_sample;
	  

	  
	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);

// 	  Eigen::MatrixXd AAA;	
// 
// 	  AAA.setZero(4,4);	
// 	  AAA(0,0) = pow(t_plan(0), 3); AAA(0,1) = pow(t_plan(0), 2); AAA(0,2) = pow(t_plan(0), 1); AAA(0,3) = pow(t_plan(0), 0); 
// 	  AAA(1,0) = pow(t_plan(1), 3); AAA(1,1) = pow(t_plan(1), 2); AAA(1,2) = pow(t_plan(1), 1); AAA(1,3) = pow(t_plan(0), 0); 
// 	  AAA(2,0) = pow(t_plan(2), 3); AAA(2,1) = pow(t_plan(2), 2); AAA(2,2) = pow(t_plan(2), 1); AAA(2,3) = pow(t_plan(0), 0); 
// 	  AAA(3,0) = 3*pow(t_plan(2), 2); AAA(3,1) = 2*pow(t_plan(2), 1); AAA(3,2) = pow(t_plan(2), 0); AAA(3,3) = 0;  


	  Eigen::Matrix4d AAA_inv;
	  
	  double abx1, abx2, abx3, abx4;
	  abx1 = ((t_plan(0) - t_plan(1))*pow(t_plan(0) - t_plan(2), 2));
	  abx2 = ((t_plan(0) - t_plan(1))*pow(t_plan(1) - t_plan(2), 2));
	  abx3 =(pow(t_plan(0) - t_plan(2), 2)*pow(t_plan(1) - t_plan(2), 2));
	  abx4 = ((t_plan(0) - t_plan(2))*(t_plan(1) - t_plan(2)));
	  

	  AAA_inv(0,0) = 1/ abx1;
	  AAA_inv(0,1) =  -1/ abx2;
	  AAA_inv(0,2) = (t_plan(0) + t_plan(1) - 2*t_plan(2))/ abx3;
	  AAA_inv(0,3) = 1/ abx4;
	  
	  AAA_inv(1,0) = -(t_plan(1) + 2*t_plan(2))/ abx1;
	  AAA_inv(1,1) = (t_plan(0) + 2*t_plan(2))/ abx2;
	  AAA_inv(1,2) = -(pow(t_plan(0), 2) + t_plan(0)*t_plan(1) + pow(t_plan(1), 2) - 3*pow(t_plan(2), 2))/ abx3;
	  AAA_inv(1,3) = -(t_plan(0) + t_plan(1) + t_plan(2))/ abx4;
	  
	  AAA_inv(2,0) = (t_plan(2)*(2*t_plan(1) + t_plan(2)))/ abx1;
	  AAA_inv(2,1) = -(t_plan(2)*(2*t_plan(0) + t_plan(2)))/ abx2;
	  AAA_inv(2,2) = (t_plan(2)*(2*pow(t_plan(0), 2) + 2*t_plan(0)*t_plan(1) - 3*t_plan(2)*t_plan(0) + 2*pow(t_plan(1), 2) - 3*t_plan(2)*t_plan(1)))/ abx3;
	  AAA_inv(2,3) = (t_plan(0)*t_plan(1) + t_plan(0)*t_plan(2) + t_plan(1)*t_plan(2))/ abx4;
	  
	  AAA_inv(3,0) = -(t_plan(1)*pow(t_plan(2), 2))/ abx1;
	  AAA_inv(3,1) = (t_plan(0)*pow(t_plan(2), 2))/ abx2;
	  AAA_inv(3,2) = (t_plan(0)*t_plan(1)*(t_plan(0)*t_plan(1) - 2*t_plan(0)*t_plan(2) - 2*t_plan(1)*t_plan(2) + 3*pow(t_plan(2), 2)))/ abx3;
	  AAA_inv(3,3) = -(t_plan(0)*t_plan(1)*t_plan(2))/ abx4;
	  	  
	  
	  
	  
	
	  
	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);   t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);   t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1);  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  // COM&&foot trajectory interpolation
	  
	  Eigen::Vector3d  x10;
	  Eigen::Vector3d  x11;
	  Eigen::Vector3d  x12;
	  Eigen::Vector3d  x13;


/*	  x10(0) = COM_IN(0,walktime-2); x10(1) = COM_IN(1,walktime-2); x10(2) = COM_IN(2,walktime-2);
	  x11(0) = COM_IN(0,walktime-1); x11(1) = COM_IN(1,walktime-1); x11(2) = COM_IN(2,walktime-1);	 */ 
	  x10 = body_in1; 
	  x11 = body_in2;  
	  x12 = _torso_angle_optimal.col(t_int);
	  x13 = _torso_angle_optimal.col(t_int+1);
	  
	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = x10(0); temp(1) = x11(0); temp(2) = x12(0); temp(3) = _thetavx(t_int);	  
	  com_inte(0) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(1); temp(1) = x11(1); temp(2) = x12(1); temp(3) = _thetavy(t_int);	  
	  com_inte(1) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(2); temp(1) = x11(2); temp(2) = x12(2); temp(3) = _thetavz(t_int);	  
	  com_inte(2) = t_a_plan *(AAA_inv)*temp;

	  
	  
	}
	else
	{
	  com_inte(0) = body_in3(0);	  
	  com_inte(1) = body_in3(1);	  	  
	  com_inte(2) = body_in3(2);	
	  
	}

 	return com_inte;
	
  
  
  
}


Vector3d MPCClass::XGetSolution_Foot_positionR(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
	
        _R_foot_optition_optimal.row(0) = _Rfootx;
	_R_foot_optition_optimal.row(1) = _Rfooty;
	_R_foot_optition_optimal.row(2) = _Rfootz;
	
	Vector3d com_inte;	
	
	if (walktime>=2)
	{
	  int t_int; 
	  t_int = floor(walktime / (_dt / dt_sample) );

	  
	  double t_cur;
	  t_cur = walktime * dt_sample;
	  

	  
	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);

// 	  Eigen::MatrixXd AAA;
// 
// 	  
// 	  AAA.setZero(4,4);	
// 	  AAA(0,0) = pow(t_plan(0), 3); AAA(0,1) = pow(t_plan(0), 2); AAA(0,2) = pow(t_plan(0), 1); AAA(0,3) = pow(t_plan(0), 0); 
// 	  AAA(1,0) = pow(t_plan(1), 3); AAA(1,1) = pow(t_plan(1), 2); AAA(1,2) = pow(t_plan(1), 1); AAA(1,3) = pow(t_plan(0), 0); 
// 	  AAA(2,0) = pow(t_plan(2), 3); AAA(2,1) = pow(t_plan(2), 2); AAA(2,2) = pow(t_plan(2), 1); AAA(2,3) = pow(t_plan(0), 0); 
// 	  AAA(3,0) = 3*pow(t_plan(2), 2); AAA(3,1) = 2*pow(t_plan(2), 1); AAA(3,2) = pow(t_plan(2), 0); AAA(3,3) = 0;  


	  
	  Eigen::Matrix4d AAA_inv;
	  
	  double abx1, abx2, abx3, abx4;
	  abx1 = ((t_plan(0) - t_plan(1))*pow(t_plan(0) - t_plan(2), 2));
	  abx2 = ((t_plan(0) - t_plan(1))*pow(t_plan(1) - t_plan(2), 2));
	  abx3 =(pow(t_plan(0) - t_plan(2), 2)*pow(t_plan(1) - t_plan(2), 2));
	  abx4 = ((t_plan(0) - t_plan(2))*(t_plan(1) - t_plan(2)));
	  

	  AAA_inv(0,0) = 1/ abx1;
	  AAA_inv(0,1) =  -1/ abx2;
	  AAA_inv(0,2) = (t_plan(0) + t_plan(1) - 2*t_plan(2))/ abx3;
	  AAA_inv(0,3) = 1/ abx4;
	  
	  AAA_inv(1,0) = -(t_plan(1) + 2*t_plan(2))/ abx1;
	  AAA_inv(1,1) = (t_plan(0) + 2*t_plan(2))/ abx2;
	  AAA_inv(1,2) = -(pow(t_plan(0), 2) + t_plan(0)*t_plan(1) + pow(t_plan(1), 2) - 3*pow(t_plan(2), 2))/ abx3;
	  AAA_inv(1,3) = -(t_plan(0) + t_plan(1) + t_plan(2))/ abx4;
	  
	  AAA_inv(2,0) = (t_plan(2)*(2*t_plan(1) + t_plan(2)))/ abx1;
	  AAA_inv(2,1) = -(t_plan(2)*(2*t_plan(0) + t_plan(2)))/ abx2;
	  AAA_inv(2,2) = (t_plan(2)*(2*pow(t_plan(0), 2) + 2*t_plan(0)*t_plan(1) - 3*t_plan(2)*t_plan(0) + 2*pow(t_plan(1), 2) - 3*t_plan(2)*t_plan(1)))/ abx3;
	  AAA_inv(2,3) = (t_plan(0)*t_plan(1) + t_plan(0)*t_plan(2) + t_plan(1)*t_plan(2))/ abx4;
	  
	  AAA_inv(3,0) = -(t_plan(1)*pow(t_plan(2), 2))/ abx1;
	  AAA_inv(3,1) = (t_plan(0)*pow(t_plan(2), 2))/ abx2;
	  AAA_inv(3,2) = (t_plan(0)*t_plan(1)*(t_plan(0)*t_plan(1) - 2*t_plan(0)*t_plan(2) - 2*t_plan(1)*t_plan(2) + 3*pow(t_plan(2), 2)))/ abx3;
	  AAA_inv(3,3) = -(t_plan(0)*t_plan(1)*t_plan(2))/ abx4;
	  	  
	  
	  	  
	  

	
	  
	  
	  Eigen::Matrix<double, 1, 4> t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);  
	  t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);  
	  t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1);  
	  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  // COM&&foot trajectory interpolation
	  
	  Eigen::Vector3d  x10;
	  Eigen::Vector3d  x11;
	  Eigen::Vector3d  x12;
	  Eigen::Vector3d  x13;


/*	  x10(0) = COM_IN(0,walktime-2); x10(1) = COM_IN(1,walktime-2); x10(2) = COM_IN(2,walktime-2);
	  x11(0) = COM_IN(0,walktime-1); x11(1) = COM_IN(1,walktime-1); x11(2) = COM_IN(2,walktime-1);	 */ 
	  x10 = body_in1; 
	  x11 = body_in2;  
	  x12 =  _R_foot_optition_optimal.col(t_int);
	  x13 =  _R_foot_optition_optimal.col(t_int+1);
	  
	  
	  Eigen::Matrix<double, 4, 1>  temp;
	  temp.setZero();
	  temp(0) = x10(0); temp(1) = x11(0); temp(2) = x12(0); temp(3) = _Rfootvx(t_int);	  
	  com_inte(0) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(1); temp(1) = x11(1); temp(2) = x12(1); temp(3) = _Rfootvy(t_int);

// 	  cout << "Rfooty:"<<temp<<endl;
	  com_inte(1) = t_a_plan * (AAA_inv)*temp;
	  temp(0) = x10(2); temp(1) = x11(2); temp(2) = x12(2); temp(3) = _Rfootvz(t_int);	  
	  com_inte(2) = t_a_plan *(AAA_inv)*temp;
	  
////=====================================================////////////////////////////////////////////////
// 	  //////spline for Rfoot_
	  
	  
	  
	  
	  
	  
	  
	  
	  
	}
	else
	{
	  com_inte(0) = body_in3(0);	  
	  com_inte(1) = body_in3(1);	  	  
	  com_inte(2) = body_in3(2);

// 	  com_inte = Rfoot_IN.col(walktime);	  
	}
	
// 	cout << "Rfooty_generated:"<<com_inte(1)<<endl;

 	return com_inte;
	
	
	
	
}

Vector3d MPCClass::XGetSolution_Foot_positionL(int walktime, double dt_sample, Eigen::Vector3d body_in1, Eigen::Vector3d body_in2, Eigen::Vector3d body_in3)
{
        _L_foot_optition_optimal.row(0) = _Lfootx;
	_L_foot_optition_optimal.row(1) = _Lfooty;
	_L_foot_optition_optimal.row(2) = _Lfootz;
	
	Vector3d com_inte;	
	
	if (walktime>=2)
	{
	  int t_int; 
	  t_int = floor(walktime / (_dt / dt_sample) );

	  
	  double t_cur;
	  t_cur = walktime * dt_sample;
	  

	  
	  Eigen::Matrix<double, 4, 1> t_plan;
	  t_plan.setZero();
	  t_plan(0) = t_cur - 2*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(1) = t_cur - 1*dt_sample-( t_cur - 2*dt_sample);
	  t_plan(2) = (t_int + 1) *_dt-( t_cur - 2*dt_sample);
	  t_plan(3) = (t_int + 2) *_dt-( t_cur - 2*dt_sample);


/*	  Eigen::MatrixXd  AAAaaa; /// should be Marix4d
	  AAAaaa.setZero(4,4);
	  
	  AAAaaa(0,0) = pow(t_plan(0), 3); AAAaaa(0,1) = pow(t_plan(0), 2); AAAaaa(0,2) = pow(t_plan(0), 1); AAAaaa(0,3) = pow(t_plan(0), 0); 
	  AAAaaa(1,0) = pow(t_plan(1), 3); AAAaaa(1,1) = pow(t_plan(1), 2); AAAaaa(1,2) = pow(t_plan(1), 1); AAAaaa(1,3) = pow(t_plan(1), 0); 
	  AAAaaa(2,0) = pow(t_plan(2), 3); AAAaaa(2,1) = pow(t_plan(2), 2); AAAaaa(2,2) = pow(t_plan(2), 1); AAAaaa(2,3) = pow(t_plan(2), 0); 
         // AAAaaa(3,0) = pow(t_plan(3), 3); AAAaaa(3,1) = pow(t_plan(3), 2); AAAaaa(3,2) = pow(t_plan(3), 1); AAAaaa(3,3) = pow(t_plan(3), 0); /// using the next to positioin would cause over-fitting	  
	  AAAaaa(3,0) = 3*pow(t_plan(2), 2); AAAaaa(3,1) = 2*pow(t_plan(2), 1); AAAaaa(3,2) = pow(t_plan(2), 0); AAAaaa(3,3) = 0.0;  	  
 */	  
//       MatrixXd.inverse( != Matrix4d (under Xd =4. so write the inverse of Matrix4d explicitly): for the time being)
	  
	  
	  Eigen::Matrix4d AAA_inv;
	  
	  double abx1, abx2, abx3, abx4;
	  abx1 = ((t_plan(0) - t_plan(1))*pow(t_plan(0) - t_plan(2), 2));
	  abx2 = ((t_plan(0) - t_plan(1))*pow(t_plan(1) - t_plan(2), 2));
	  abx3 =(pow(t_plan(0) - t_plan(2), 2)*pow(t_plan(1) - t_plan(2), 2));
	  abx4 = ((t_plan(0) - t_plan(2))*(t_plan(1) - t_plan(2)));
	  

	  AAA_inv(0,0) = 1/ abx1;
	  AAA_inv(0,1) =  -1/ abx2;
	  AAA_inv(0,2) = (t_plan(0) + t_plan(1) - 2*t_plan(2))/ abx3;
	  AAA_inv(0,3) = 1/ abx4;
	  
	  AAA_inv(1,0) = -(t_plan(1) + 2*t_plan(2))/ abx1;
	  AAA_inv(1,1) = (t_plan(0) + 2*t_plan(2))/ abx2;
	  AAA_inv(1,2) = -(pow(t_plan(0), 2) + t_plan(0)*t_plan(1) + pow(t_plan(1), 2) - 3*pow(t_plan(2), 2))/ abx3;
	  AAA_inv(1,3) = -(t_plan(0) + t_plan(1) + t_plan(2))/ abx4;
	  
	  AAA_inv(2,0) = (t_plan(2)*(2*t_plan(1) + t_plan(2)))/ abx1;
	  AAA_inv(2,1) = -(t_plan(2)*(2*t_plan(0) + t_plan(2)))/ abx2;
	  AAA_inv(2,2) = (t_plan(2)*(2*pow(t_plan(0), 2) + 2*t_plan(0)*t_plan(1) - 3*t_plan(2)*t_plan(0) + 2*pow(t_plan(1), 2) - 3*t_plan(2)*t_plan(1)))/ abx3;
	  AAA_inv(2,3) = (t_plan(0)*t_plan(1) + t_plan(0)*t_plan(2) + t_plan(1)*t_plan(2))/ abx4;
	  
	  AAA_inv(3,0) = -(t_plan(1)*pow(t_plan(2), 2))/ abx1;
	  AAA_inv(3,1) = (t_plan(0)*pow(t_plan(2), 2))/ abx2;
	  AAA_inv(3,2) = (t_plan(0)*t_plan(1)*(t_plan(0)*t_plan(1) - 2*t_plan(0)*t_plan(2) - 2*t_plan(1)*t_plan(2) + 3*pow(t_plan(2), 2)))/ abx3;
	  AAA_inv(3,3) = -(t_plan(0)*t_plan(1)*t_plan(2))/ abx4;
	  	  
	  
	  
	  
	  
	
// // 	  Eigen::Matrix4d  AAAaaa1_inv=AAA_inv;
	  
// 	  cout<< t_plan<<endl;
// 	  cout<< AAAaaa<<endl;
// 	  cout<< AAA_inv - AAAaaa1_inv<<endl;
// 	  cout<< AAA_inv - AAAaaa1.inverse()<<endl;
	  
	  Eigen::RowVector4d t_a_plan;
	  t_a_plan.setZero();
	  t_a_plan(0) = pow(t_cur-( t_cur - 2*dt_sample), 3);  
	  t_a_plan(1) = pow(t_cur-( t_cur - 2*dt_sample), 2);   
	  t_a_plan(2) = pow(t_cur-( t_cur - 2*dt_sample), 1); 
	  t_a_plan(3) = pow(t_cur-( t_cur - 2*dt_sample), 0); 


	  
	  // COM&&foot trajectory interpolation
	  
	  Eigen::Vector3d  x10;
	  Eigen::Vector3d  x11;
	  Eigen::Vector3d  x12;
// 	  Eigen::Vector3d  x13;


/*	  x10(0) = COM_IN(0,walktime-2); x10(1) = COM_IN(1,walktime-2); x10(2) = COM_IN(2,walktime-2);
	  x11(0) = COM_IN(0,walktime-1); x11(1) = COM_IN(1,walktime-1); x11(2) = COM_IN(2,walktime-1);	 */ 
	  x10 = body_in1; 
	  x11 = body_in2;  
	  x12 =  _L_foot_optition_optimal.col(t_int);
// 	  x13 =  _L_foot_optition_optimal.col(t_int+1); /// using next two position would caused over-fitting
	  
	  
	  Eigen::Vector4d  temp;
	  temp.setZero();
	  temp(0) = x10(0); temp(1) = x11(0); temp(2) = x12(0);
	  temp(3) = _Lfootvx(t_int);
// 	  	  temp(3) = x13(0);
	  Eigen::Vector4d tmp111 = AAA_inv*temp;
	  com_inte(0) = t_a_plan * tmp111;
	  temp(0) = x10(1); temp(1) = x11(1); temp(2) = x12(1);
	  temp(3) = _Lfootvy(t_int);	 
/*          temp(3) = x13(1);	*/  
	  tmp111 = AAA_inv*temp;  
	  com_inte(1) = t_a_plan * tmp111;
	  temp(0) = x10(2); temp(1) = x11(2); temp(2) = x12(2); 
	  temp(3) = _Lfootvz(t_int);	
// 	  temp(3) = x13(2);
	  tmp111 = AAA_inv*temp;	  
	  com_inte(2) = t_a_plan * tmp111;

////================================================not use spline.h=====////////////////////////////////////////////////


	  
	}
	else
	{
	  com_inte(0) = body_in3(0);	  
	  com_inte(1) = body_in3(1);	  	  
	  com_inte(2) = body_in3(2);
	}

 	return com_inte;
	
	cout << "Lfooty_generated:"<<com_inte(1)<<endl;
}






