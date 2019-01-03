/*****************************************************************************
ComanSimClass.h

Description:	Header file of ComanSimClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2016/03/21
@Update:	2016/03/21
*****************************************************************************/
#ifndef COMAN_SIM_CLASS_H
#define COMAN_SIM_CLASS_H

#include "ComanODE/ODESimBaseClass.h"
#include "ComanODE/ComanODERobotClass.h"
#include "ComanODE/BigmanODERobotClass.h"
#include "ComanODE/CogimonODERobotClass.h"



class ComanSimClass: public ODESimBaseClass
{
public:
	ComanSimClass();
	~ComanSimClass();

	void command(int cmd);
	void nearCallback(void *data, dGeomID o1, dGeomID o2);
	void simLoop(int pause);
	void UpdateODE();

private:

	double tn, fn, Amp, offset, wn, cpg, anti_cpg;

	int IsStart;

	void init();

	void makeRobot();
	void restartODE();
	void initExternalForce(const int& link_name, const Eigen::Vector3d& amp);
	void initExternalTorque(const int& link_name, const Eigen::Vector3d& amp);
	void applyExternalForce();
	void applyBallAttack();

	void drawToughTerrain();
	void makeToughTerrain();

	bool IsDropDownRobot;
	std::vector<dJointID> fixedhead;

	void makeObj();
	std::vector<dJointID> fixedjntID;
	std::vector<RobotJoint> objfixedjnts;
	std::vector<RobotLink> objs;
	std::vector<dJointFeedback> objs_jnt_fb;

	void makeObj2();
};


#endif
