/***************************************************************************** 
ODEsettings.h

Description:	Header file of ODEsettings

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2016/03/18
@Update:	2016/03/18
*****************************************************************************/ 
#ifndef ODE_SETTINGS_H
#define ODE_SETTINGS_H

#include <stdio.h>
#include <sys/types.h>
#include "ode/ode.h"
#include "drawstuff/drawstuff.h"

#ifdef dDOUBLE
#define dsDrawBox      dsDrawBoxD
#define dsDrawSphere   dsDrawSphereD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawCapsule  dsDrawCapsuleD
#define dsDrawLine    dsDrawLineD
#endif


#define JOINT_SIZE	      0.0001
#define JOINT_MASS	      0.0001
#define JOINT_RADIUS	  0.05
#define FUDGE_FACTOR      0.1

#define HEAD_RADIUS       0.1

#define HAND_LENGTH     0.04
#define HAND_WIDTH      0.01
#define HAND_HEIGHT     0.06

// const double DELTAtime = deltaT;
const double GlobalERP = 0.2;
const double GlobalCFM = 0.0001;







#endif
