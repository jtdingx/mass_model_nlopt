/*****************************************************************************
Comanoderobot[0]Class.cpp

Description:  cpp file of Comanoderobot[0]Class

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: 2016/03/18
@Update:  2016/03/20
*****************************************************************************/

#include "ComanODE/ODESimBaseClass.h"

ODESimBaseClass::ODESimBaseClass()
  : loop(0.0)
  , sim_t(0.0)
{
  // oderobot.push_back(boost::shared_ptr<ComanODERobotClass>(new ComanODERobotClass(0.55, 0, 0, 0)));
  // dt = oderobot[0]->dt;
  // g = oderobot[0]->g;

// #ifdef PUBLISH_TO_ROS
//   int argc = 0;
//   char **argv = NULL;
//   ros::init(argc, argv, robotpara.name);
// #endif

};

void ODESimBaseClass::UpdateODE()
{
  if (!oderobot.empty()) {
    for (auto it : oderobot)
    {
      it->controlMotors("position");
    }
  }
  // dWorldQuickStep(world, dt);
  dWorldStep(world, dt);
  dJointGroupEmpty(contactgroup);

  loop++;

  if (!oderobot.empty()) {
    for (auto it : oderobot) {
      it->drawRobot();
    }
  }
}

void ODESimBaseClass::nearCallback(void *data, dGeomID o1, dGeomID o2)
{
  const int contact_no = 5;

  dBodyID b1 = dGeomGetBody(o1);
  dBodyID b2 = dGeomGetBody(o2);
  if (b1 && b2 && dAreConnected(b1, b2)) return;

  dContact contact[contact_no];
  int numc = dCollide(o1, o2, contact_no, &contact[0].geom, sizeof(dContact));

  for (int i = 0; i < numc; i++)
  {
    contact[i].surface.mode     = dContactSoftERP | dContactSoftCFM | dContactApprox1;
    contact[i].surface.soft_erp   = GlobalERP;
    contact[i].surface.soft_cfm   = GlobalCFM;
    contact[i].surface.mu     = dInfinity;
    contact[i].surface.mu2      = 0;
    contact[i].surface.bounce   = 0.0;      // bouncing the objects
    contact[i].surface.bounce_vel = 0.0;      // bouncing velocity

    dJointID c = dJointCreateContact(world, contactgroup, &contact[i]);
    dJointAttach(c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
  }
}


void ODESimBaseClass::command(int cmd)
{
  switch (cmd)
  {
  case 'q':
    dsStop(); break;
  }
}

void ODESimBaseClass::changeCameraPosition(dBodyID id, double x, double y, double z, double h, double p, double r)
{
  const dReal *pos;
  pos = dBodyGetPosition(id);

  float xyz[3] = {float(pos[0] + x), float(pos[1] + y), float(pos[2] + z)};
  float hpr[3] = {float(h), float(p), float(r)};

  dsSetViewpoint(xyz, hpr);
  dsSetSphereQuality(3);
}

void ODESimBaseClass::BallAttack(const dBodyID &link_id, const Eigen::Vector3d &rel_pos, RobotLink& obj)
{
  obj.id = dBodyCreate(world);

  const dReal *pos = dBodyGetPosition(link_id);
  dBodySetPosition(obj.id, pos[0] + rel_pos(0), pos[1] + rel_pos(1), pos[2] + rel_pos(2));
  dMass m;
  dMassSetSphereTotal(&m, obj.m, obj.r);
  dMassAdjust(&m, obj.m);
  dBodySetMass(obj.id, &m);
  obj.gid = dCreateSphere(space, obj.r);
  dGeomSetBody(obj.gid, obj.id);

  // Eigen::Vector3d vel = obj.lx * (-rel_pos / rel_pos.norm());
  Eigen::Vector3d vel(obj.lx, obj.ly, 0.0);

  dBodySetLinearVel(obj.id, vel(0), vel(1), vel(2));

  // dsSetColor(obj.color_r, obj.color_g, obj.color_b);
  // dsDrawSphere(dBodyGetPosition(obj.id), dBodyGetRotation(obj.id), obj.r);

}

void ODESimBaseClass::start()
{
  float xyz[3] = {  1.6f, -1.0f, 0.5f};
  float hpr[3] = {150.0f,  0.0f, 0.0f};

  //// along x axis
  //double xyz[3] = {  1.4f, 0.0f, 0.5f};
  //double hpr[3] = {180.0f,   0.0f, 0.0f};

  dsSetViewpoint(xyz, hpr);
  dsSetSphereQuality(3);
}

void ODESimBaseClass::makeRobot()
{
  if (!oderobot.empty()) {
    for (auto it : oderobot)
    {
      it->makeRobot(world, space);
    }
  }
}

int ODESimBaseClass::main(dsFunctions &fn)
{
  int argc = 0;
  char **argv = NULL;
  // ros::init(argc, argv, robotpara.name);
  // ros::NodeHandle n;
  // _nh = &n;

  dInitODE();
  world  = dWorldCreate();
  space  = dHashSpaceCreate(0);
  ground = dCreatePlane(space, 0, 0, 1, 0);
  contactgroup = dJointGroupCreate(0);


  dWorldSetERP(world, GlobalERP);
  dWorldSetCFM(world, GlobalCFM);

  makeRobot();

  dWorldSetGravity(world, 0, 0, -g);

  dsSimulationLoop(argc, argv, 600, 800, &fn);
  // dsSimulationLoop(argc, argv, 1280, 720, &fn);

  dJointGroupDestroy(contactgroup);
  dSpaceDestroy(space);
  dWorldDestroy(world);
  dCloseODE();
  return 0;
}



