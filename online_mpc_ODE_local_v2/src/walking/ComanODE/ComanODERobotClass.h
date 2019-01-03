/*****************************************************************************
ComanODERobotClass.h

Description:	Header file of ComanODERobotClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2016/03/21
@Update:	2016/03/21
*****************************************************************************/
#ifndef COMAN_ODE_ROBOT_CLASS_H
#define COMAN_ODE_ROBOT_CLASS_H

#include "ComanODE/ODERobotBaseClass.h"


class ComanODERobotClass: public ODERobotBaseClass
{
public:
	// ComanODERobotClass(const RobotParaClass& robotpara, const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	ComanODERobotClass(const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	~ComanODERobotClass();

	void drawRobot();

private:
	void updateLinkParam(const RobotParaClass& robotpara);
	void updateJointParam(const RobotParaClass& robotpara);
};


#endif
