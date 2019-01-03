/*****************************************************************************
ODERobotBaseClass.h

Description:	Header file of ODERobotBaseClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2016/03/18
@Update:	2016/03/28
*****************************************************************************/
#ifndef ODE_ROBOT_CLASS_H
#define ODE_ROBOT_CLASS_H

#include "ComanODE/ODEsettings.h"
#include "RobotPara/RobotParaClass.h"

#include <vector>
#include <deque>

#include "utils/IOClass.h"

#ifdef PUBLISH_TO_ROS
#include "utils/ChartClass.h"
#endif

#include "RTControl/MpcRTControlClass.h"


enum joint_type
{
	BALL,
	HINGE,
	SLIDER,
	HINGE2,
	FIXED,
};

enum body_type
{
	BOX,
	SPHERE,
	CAPSULE,
	CYLINDER,
};

typedef struct {
	bool IsExist = false;
	std::string name;
	dBodyID id;
	dGeomID gid;
	int 	type; 			//!< collsion body type
	int 	ds_type; 		//!< drawstuff body type
	dReal   px, py, pz; 	//!< position (center of gravity) x,y,z
	dReal   lx, ly, lz; 	//!< length x axis, y axis, z axis
	dReal 	cgx, cgy, cgz; 	//!< center of gravity position in body frame (x,y,z)
	dReal   I11, I22, I33;
	dReal   I12, I13, I23;
	dReal   r;         		//!< radius
	dReal   m = JOINT_MASS;         		//!< total mass of the rigid body
	dReal	color_r, color_g, color_b, color_a; //!< red/green/blue/alpha, each in the range of [0,1]
	int direction; //!< The cylinder's long axis is oriented along the body's x, y or z axis according to the value of direction (1=x, 2=y, 3=z)
} RobotLink;

typedef struct {
	bool IsExist = false;
	std::string name;
	dJointID id;
	int   type;                   //!< 0:ball 1:hinge 2:slider 3:hinge2
	int	  link[2];            	  //!< Two links should be attached. 0:parent, 1:child
	RobotLink* links[2];
	dReal px, py, pz;             //!< anchor position x,y,z
	dReal axis_x, axis_y, axis_z; //!< rotation axis
	dReal lo_stop, hi_stop;       //!< low stop, high stop
	dReal fmax;                   //!< max force
	dReal vmax;                   //!< max velocity
	dReal fudge_factor;           //!< fudge factor
	dReal bounce;                 //!< boucyness of the stop
} RobotJoint;


class ODERobotBaseClass
{
public:
	ODERobotBaseClass(const double &z, const double &y, const double &x, const double &yaw);
	ODERobotBaseClass(const RobotParaClass& robotpara, const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	virtual ~ODERobotBaseClass() {};

	virtual void controlMotors(std::string method);
	virtual void controlMotor(std::string method, const int& joint_name);

	virtual void makeRobot(const dWorldID &world, const dSpaceID &space);
	static void RotZ(double& px, double& py, const double& yaw);

	const RobotParaClass& RobotPara() {return mComanWBS->RobotPara();};;

	bool IsFinishHoming;

	double dt;
	double g;
	double start_x, start_y, start_z;
	double start_yaw;
	double JointNUM;
	std::string name;

	std::vector<RobotLink>	rlink;
	std::vector<RobotJoint>	rjoint;
	std::vector<double>		joint_angle_ref, joint_torque_ref, joint_vel_ref, joint_torque_comp_ref;
	std::vector<double> 	l_foot_FT, r_foot_FT, l_hand_FT, r_hand_FT; // force torque snesor

	virtual void UpdateSensorFeedback();
	virtual inline void KeyBoardControl(const int& cmd) {RTControl.KeyBoardControl(cmd);};
	virtual inline void savedata() {RTControl.savedata();};

	inline const RobotStateClass& getRobotState() {return mComanWBS->getRobotState();}
	inline void UpdateFallState(const FallState &state) {mComanWBS->UpdateFallState(state);};

	Eigen::Vector3d gcom, gdcom, gcom_old;

	virtual void drawRobot() = 0;

	void makeLink(RobotLink& robotlink, const dWorldID &world, const dSpaceID &space);
	void drawLink(const RobotLink& robotlink);

	static std::vector<double> getLocalJointForceTorque(const RobotJoint& robotjoint);

	inline bool IsInitRTctrl() {return RTControl.IsInit;};
	inline void AddObstacle(const Eigen::Vector3d& obs_center_pos, const Eigen::Vector3d& obs_center_ori) {RTControl.AddObstacle(obs_center_pos, obs_center_ori);};

	MpcRTControlClass RTControl;

private:
	void Init(const RobotParaClass& robotpara, const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	void readLinkParam(const RobotParaClass& robotpara);
	void readJointParam(const RobotParaClass& robotpara);

	std::vector<double> getJointPosition(const RobotParaClass& robotpara, const int& joint_name);


	std::vector<double> torque_diff_old;
	dJointFeedback *fb;

	std::vector<FilterClass> torque_fb_filter;

	double totalmass;

#ifdef PUBLISH_TO_ROS
	boost::shared_ptr<ChartClass> data_publisher;
#endif

	void ReadForeTorqueSensor();
	void CalcGlobalCOM();
	void PublishToRos();
	void UpdateWBS();

protected:
	WholeBodySensingClass *mComanWBS;

	std::deque<Eigen::Vector3d> COMqueue;
	std::vector<dJointFeedback> jnt_fb;
	// virtual void updateLinkParam(const RobotParaClass& robotpara) = 0;
	// virtual void updateJointParam(const RobotParaClass& robotpara) = 0;
	void updateRobotYaw();

	void PositionMotorControl(const int& joint_name, dReal joint_angle);
	void TorqueMotorControl(const int& joint_name, const dReal& torque);
	void ImpedanceMotorControl(const int& joint_name, dReal joint_angle, dReal joint_vel = 0.0, const dReal& tau_ff = 0.0);
	void VelocityMotorControl(const int& joint_name, const dReal& joint_angular_velocity);
	void AdmittanceMotorControl(const int& joint_name, dReal joint_angle);
	void drawCOPLine();

};


#endif
