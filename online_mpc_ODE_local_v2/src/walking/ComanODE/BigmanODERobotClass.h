/*****************************************************************************
BigmanODERobotClass.h

Description:	Header file of BigmanODERobotClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	Wed 31 Aug 2016 06:40:23 PM CEST
@Update:	Wed 31 Aug 2016 06:40:19 PM CEST
*****************************************************************************/
#pragma once

#include "ComanODE/ODERobotBaseClass.h"
// #include "ComanODE/ODERobotBaseEmptyClass.h"


class BigmanODERobotClass: public ODERobotBaseClass
{
public:
	BigmanODERobotClass(const double &z = 0.6, const double &y = 0, const double &x = 0, const double &yaw = 0);
	~BigmanODERobotClass();

	void drawRobot();

private:
	void updateLinkParam(const RobotParaClass& robotpara);
	void updateJointParam(const RobotParaClass& robotpara);
};
