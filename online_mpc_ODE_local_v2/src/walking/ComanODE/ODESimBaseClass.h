/*****************************************************************************
ODESimBaseClass.h

Description:	Header file of ODESimBaseClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2016/03/21
@Update:	2016/03/21
*****************************************************************************/
#ifndef ODE_SIM_BASE_CLASS_H
#define ODE_SIM_BASE_CLASS_H

#include "ComanODE/ODERobotBaseClass.h"

class ExternalForceProfile
{
public:
	ExternalForceProfile()
		: _IsApplying(false)
		, _start_time(0.0)
		, _duration(0.0)
		, _amplitude(0.0, 0.0, 0.0)
		, _mode("force")
	{};
	virtual ~ExternalForceProfile() {};

	void init(const dBodyID& link_id, const std::string mode, const double& startT, const Eigen::Vector3d& amplitude, const double& duration)
	{
		if (!_IsApplying) {
			_IsApplying = true;
			_link_id 	= link_id;
			_mode		= mode;
			_start_time = startT;
			_amplitude	= amplitude;
			_duration 	= duration;
		}
	};

	void ApplyImpulse(const double& currentT)
	{
		Eigen::Vector3d force = GenerateImpulse(currentT, _IsApplying, _start_time, _amplitude, _duration);
		if (_mode == "force") {
			dBodyAddForce(_link_id, force(0), force(1), force(2));
		}
		else if (_mode == "torque") {
			dBodyAddTorque(_link_id, force(0), force(1), force(2));
		}
	}

private:
	dBodyID _link_id;
	bool 	_IsApplying;
	double 	_start_time;
	double 	_duration;
	std::string _mode;
	Eigen::Vector3d _amplitude;

	Eigen::Vector3d GenerateImpulse(const double& currentT, bool& flag, const double& start_time, const Eigen::Vector3d& goal, const double& duration)
	{
		Eigen::Vector3d result(0.0, 0.0, 0.0);

		if (flag) {
			if ((currentT <= (start_time + duration)) && duration != 0) {
				double passed_time = currentT - start_time;
				result(0) = goal(0) * sin(1 * M_PI / duration * passed_time);
				result(1) = goal(1) * sin(1 * M_PI / duration * passed_time);
				result(2) = goal(2) * sin(1 * M_PI / duration * passed_time);
			}
			else {
				flag = false;
			}
		}
		return result;
	}

};



class ODESimBaseClass
{
public:
	ODESimBaseClass();
	virtual ~ODESimBaseClass() {};

	int main(dsFunctions &fn);

	std::vector<boost::shared_ptr<ODERobotBaseClass> > oderobot;

	virtual void start();
	virtual void command(int cmd);
	virtual void nearCallback(void *data, dGeomID o1, dGeomID o2);
	virtual void simLoop(int pause) = 0;
	virtual void UpdateODE();

	inline dSpaceID getSpaceID() {return space;};


protected:
	dWorldID world;
	dSpaceID space;
	dGeomID  ground;
	dJointGroupID contactgroup;

	double loop;
	double sim_t;
	double dt;
	double g;

	virtual void makeRobot();
	void changeCameraPosition(dBodyID id, double x, double y, double z, double h, double p, double r);

	std::vector<RobotLink> balls;
	/**
	 * @brief      throw a ball to the robot link
	 *
	 * @param[in]  link_id  { goal link body id }
	 * @param[in]  rel_pos  { ball initial position relative to the link }
	 * @param      obj      { need to specify obj.m obj.r obj.color and obj.lx for initial linear velocity}
	 */
	void BallAttack(const dBodyID &link_id, const Eigen::Vector3d &rel_pos, RobotLink& obj);

	std::vector<ExternalForceProfile> extforce;

private:
	// ros::NodeHandle _nh;

};


#endif
