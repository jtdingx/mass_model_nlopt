/*****************************************************************************
CogimonODERobotClass.h

Description:	Header file of CogimonODERobotClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	Thu 30 Nov 2017 08:29:55 PM CET
@Update:	Thu 30 Nov 2017 08:30:00 PM CET
*****************************************************************************/
#pragma once

#include "ComanODE/ODERobotBaseClass.h"


class CogimonODERobotClass: public ODERobotBaseClass
{
public:
	CogimonODERobotClass(const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	~CogimonODERobotClass();

	void drawRobot();

private:
	void updateLinkParam(const RobotParaClass& robotpara);
	void updateJointParam(const RobotParaClass& robotpara);
};
