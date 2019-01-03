/*****************************************************************************
ComanSimClass.cpp

Description:  cpp file of ComanSimClass

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: 2016/03/18
@Update:  2016/03/20
*****************************************************************************/

#include "ComanODE/ComanSimClass.h"

ComanSimClass::ComanSimClass()
  : ODESimBaseClass()
  , IsDropDownRobot(false)
  , IsStart(1)
{
  objs_jnt_fb.resize(4);

};

ComanSimClass::~ComanSimClass()
{
  std::cout << "=================== ODE Simulation Ended, 888 ================" << std::endl;
};
double x_push = 0.0;  // x: +-130 is the fall limit when applied at torso
double y_push = 0.0;  // y: +-200 is the fall limit when applied at torso
double z_push = 0.0;
//////////////// Main Loop Function /////////////////////////
// #define WhichRobot 16 // coman x no stabilizer with stepping 0.2 m foot length

#if WhichRobot < 20
static double Dt_imp = 0.112; // coman
#else
static double Dt_imp = 0.1558; // walkman
#endif

// coman x
#if WhichRobot == 16
static double Df_imp = 61.6163;
#endif

void ComanSimClass::simLoop(int pause)
{
  if (IsStart == 1)
  {
    applyExternalForce();
    sim_t = loop * dt;

    if (loop == 0) {
      changeCameraPosition(oderobot[0]->rlink[RIGHT_FOOT].id, 0.0, -1.0, 0.5, 90, 0, 0); // side view
    }

    if (loop == int(round(0.3 / dt))) {
      command('t');
    }
    
 ////////////////////push for CoMAN 
//     if (loop == int(round(3.5/ dt))) {
// /*      x_push = 300;  //maximum   */ 
// 
//       x_push = 250;   
//       y_push = 0;
//       command('1');
//       cout<<"x_push:"<<x_push<<"N"<<endl;
//       cout<<"y_push:"<<y_push<<"N"<<endl;
//     }
//  
//     if (loop == int(round(6.6/ dt))) {
// //       x_push = -250; //maximum
//       x_push = -150; 
//       y_push = 0;
//       command('1');
//       cout<<"x_push:"<<x_push<<"N"<<endl;
//       cout<<"y_push:"<<y_push<<"N"<<endl;
//     }
//     
//    if (loop == int(round(9.8/ dt))) {
//       x_push = 0;
//       y_push = -100;
//       command('1');
//       cout<<"x_push:"<<x_push<<"N"<<endl;
//       cout<<"y_push:"<<y_push<<"N"<<endl;
//     }
//     
//     if (loop == int(round(12.4/ dt))) {
//       x_push = 0;
//       y_push = 100;
//       command('1');
//       cout<<"x_push:"<<x_push<<"N"<<endl;
//       cout<<"y_push:"<<y_push<<"N"<<endl;
//     }   
//       


    if (loop == int(round((8.0) / dt))) {
//       command('q');
    }
    
    
    

    for (auto it : oderobot)
    {
      // for wbs
      it->UpdateSensorFeedback();
      const RobotStateClass& irobot = it->getRobotState();

      if (!it->IsFinishHoming) { // for ODE homing
        const std::vector<double>& OffsetAng = it->RobotPara().HOMING_POS();
        for (int i = 0; i < 31; i++) {
          it->joint_angle_ref[i ] = DEGTORAD(OffsetAng[i]);
        }
        double homing_time = 0.2;
        if (sim_t >= homing_time) {
          it->IsFinishHoming = true;
          // it->RTControl.HomingInit(it->joint_angle_ref);
        }
        it->controlMotors("position");
      }
      else {
        for (int i = 0; i < RobotParaClass::JOINT_NUM(); i++) {
          if (i >= RIGHT_HIP_PITCH &&
              i <= LEFT_FOOT_PITCH) {
            it->controlMotor("position", i);
          }
          else {
            it->controlMotor("position", i);
          }
        }
      }

      if (it->IsInitRTctrl())
      {
        it->controlMotors("position");
        // it->controlMotors("torque");
      }
      // it->controlMotors("impedance");
      // it->controlMotors("admittance");
    }

  }
}

void ComanSimClass::UpdateODE()
{
  if (IsStart == 1) {
    dWorldStep(world, dt);
    dJointGroupEmpty(contactgroup);
//     if (command('c'))
//     {
//     loop++;      
//     }
      
    loop++;
  }

  for (auto it : oderobot) {
    it->drawRobot();
    if (!objs.empty() && objs[0].id) {
      for (auto it_obj : objs) {
        it->drawLink(it_obj);
      }
    }
  }

  drawToughTerrain();

  if (!balls.empty() && balls[0].id) {
    for (auto it : balls) {
      dsSetColor(it.color_r, it.color_g, it.color_b);
      dsDrawSphere(dBodyGetPosition(it.id), dBodyGetRotation(it.id), it.r);
    }
  }
}



/////push 

double interval = 10;
double z_turn = 0.0;
void ComanSimClass::command(int cmd)
{
  for (auto it : oderobot) {
    it->KeyBoardControl(cmd);
  }

//   int target = TORSO;
  // interval = 150;
  int target = PELVIS;
  // int target = 100; // >=100 means applying force on objs

  switch (cmd)
  {
  case 'q':
    for (auto it : oderobot) {
      // it->savedata();
    }
    RobotParaClass::robot_push_force.clear();
    dsStop();
    std::cout << "The simulation runs " << loop*dt << " s." << std::endl;
    break;

  // case 'p': case 'P':
  //   IsStart *= -1;
  //   break;

  case 'r': restartODE(); break;

  case '1':
  {
    Eigen::Vector3d amplitude(x_push, y_push, z_push);
    initExternalForce(target, amplitude);
    for (auto it : oderobot) {
      const dReal *pos = dBodyGetPosition(it->rlink[target].id);
      COUT(it->rlink[target].name, pos[2]);
    }
    // initExternalForce((int) LEFT_FOOT, -amplitude);
  }
  break;

  case '2':
  {
    Eigen::Vector3d amplitude(-x_push, -y_push, -z_push);
    initExternalForce(target, amplitude);
    // Eigen::Vector3d amplitude(0.03 * x_push, 0.03 * y_push, 0.03 * z_push);
    // initExternalTorque(target, amplitude);
  }
  break;

  case '3':
  {
    x_push = 0.0;
    y_push = 0.0;
    z_push = 0.0;
  }
  break;


  case '4':
  {
    x_push += interval;
    COUT("x_push", x_push);
  }
  break;

  case '5':
  {
    x_push -= interval;
    COUT("x_push", x_push);
  }
  break;

  case '7':
  {
    y_push += interval;
    COUT("y_push", y_push);
  }
  break;

  case '8':
  {
    y_push -= interval;
    COUT("y_push", y_push);
  }
  break;

  case '9':
  {
    z_turn += 0.5 * interval;
    z_push += interval;
    COUT("z_push", z_push);
  }
  break;

  case '0':
  {
    z_turn -= 0.5 * interval;
    z_push -= interval;
    COUT("z_push", z_push);
  }
  break;

  case 'u':
  {
    // makeObj();
    // makeObj2(); // make table
  }
  break;

  case 't':
    if (!IsDropDownRobot) {
      for (auto it : fixedhead) {
        dJointDestroy(it);
      }
      IsDropDownRobot = true;
    }
    break;
  case 'w':
    applyBallAttack();
    break;

  default: break;
  }
}


void ComanSimClass::nearCallback(void *data, dGeomID o1, dGeomID o2)
{
  dBodyID b1 = dGeomGetBody(o1);
  dBodyID b2 = dGeomGetBody(o2);
  if (b1 && b2 && dAreConnected(b1, b2)) return;

  for (auto it : oderobot)
  {
    if (
      (b1 == it->rlink[RIGHT_FOOT_DUMMY].id) ||
      (b1 == it->rlink[LEFT_FOOT_DUMMY].id) ||
      (b1 == it->rlink[WAIST_DUMMY1].id) ||
      (b1 == it->rlink[RIGHT_HAND_DUMMY1].id) ||
      (b1 == it->rlink[LEFT_HAND_DUMMY1].id) ||
      (b1 == it->rlink[NECK].id) ||
      (b1 == it->rlink[HEAD].id)
    )
      return;

    if (((b1 == it->rlink[WAIST_DUMMY2].id) && (b2 == it->rlink[PELVIS].id)) ||
        ((b2 == it->rlink[WAIST_DUMMY2].id) && (b1 == it->rlink[PELVIS].id))) return;

    if (((b1 == it->rlink[RIGHT_UPPER_ARM_DUMMY2].id) && (b2 == it->rlink[TORSO].id)) ||
        ((b2 == it->rlink[RIGHT_UPPER_ARM_DUMMY2].id) && (b1 == it->rlink[TORSO].id))) return;
    if (((b1 == it->rlink[LEFT_UPPER_ARM_DUMMY2].id) && (b2 == it->rlink[TORSO].id)) ||
        ((b2 == it->rlink[LEFT_UPPER_ARM_DUMMY2].id) && (b1 == it->rlink[TORSO].id))) return;

    if (((b1 == it->rlink[RIGHT_THIGH_DUMMY2].id) && (b2 == it->rlink[PELVIS].id)) ||
        ((b2 == it->rlink[RIGHT_THIGH_DUMMY2].id) && (b1 == it->rlink[PELVIS].id))) return;
    if (((b1 == it->rlink[LEFT_THIGH_DUMMY2].id) && (b2 == it->rlink[PELVIS].id)) ||
        ((b2 == it->rlink[LEFT_THIGH_DUMMY2].id) && (b1 == it->rlink[PELVIS].id))) return;

    // Tue 17 Jan 2017 07:13:41 PM CET
    // There are collisions between the following links, cause internal forces that move the arm outside during pure force control
    if (((b1 == it->rlink[LEFT_UPPER_ARM_DUMMY1].id) && (b2 == it->rlink[LEFT_UPPER_ARM].id)) ||
        ((b2 == it->rlink[LEFT_UPPER_ARM_DUMMY1].id) && (b1 == it->rlink[LEFT_UPPER_ARM].id))) return;
    if (((b1 == it->rlink[RIGHT_UPPER_ARM_DUMMY1].id) && (b2 == it->rlink[RIGHT_UPPER_ARM].id)) ||
        ((b2 == it->rlink[RIGHT_UPPER_ARM_DUMMY1].id) && (b1 == it->rlink[RIGHT_UPPER_ARM].id))) return;


    // Thu 12 Jan 2017 03:28:48 PM CET
    // I found there are collisions between the HAND and FORE_ARM which cause big force reading in the Hands F/T sensors
    if (((b1 == it->rlink[RIGHT_HAND].id) && (b2 == it->rlink[RIGHT_FORE_ARM].id)) ||
        ((b2 == it->rlink[RIGHT_HAND].id) && (b1 == it->rlink[RIGHT_FORE_ARM].id))) return;
    if (((b1 == it->rlink[LEFT_HAND].id) && (b2 == it->rlink[LEFT_FORE_ARM].id)) ||
        ((b2 == it->rlink[LEFT_HAND].id) && (b1 == it->rlink[LEFT_FORE_ARM].id))) return;

    if (!objs.empty() && objs[0].id) {
      for (auto it_obj : objs) {
        if ((b1 == it_obj.id) || (b2 == it_obj.id))
          return;
      }
    }
  }


  const int contact_no = 5;
  dContact contact[contact_no];
  int numc = dCollide(o1, o2, contact_no, &contact[0].geom, sizeof(dContact));

  bool show_contacts = true;
  for (int i = 0; i < numc; i++)
  {
    contact[i].surface.mode     = dContactSoftERP | dContactSoftCFM | dContactApprox1;
    contact[i].surface.soft_erp   = GlobalERP;
    contact[i].surface.soft_cfm   = GlobalCFM;
    contact[i].surface.mu     = dInfinity;
    contact[i].surface.mu2      = 0;
    contact[i].surface.bounce   = 0.0;      // bouncing the objects
    contact[i].surface.bounce_vel = 0.0;      // bouncing velocity

    //dJointID c = dJointCreateContact(world, contactgroup, contact+i);
    //dJointAttach(c, b1, b2);
    dJointID c = dJointCreateContact(world, contactgroup, &contact[i]);
    dJointAttach(c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
    if (show_contacts)
    {
      dMatrix3 RI;
      dRSetIdentity (RI);
      const dReal ss[3] = {0.01, 0.01, 0.01};
      dsSetColorAlpha (0, 0, 1, 0.5);
      dsDrawBox (contact[i].geom.pos, RI, ss);
      dReal *pos  = contact[i].geom.pos;
      dReal depth = contact[i].geom.depth;
      dReal *norm = contact[i].geom.normal;
      dReal endp[3] = {pos[0] + depth*norm[0], pos[1] + depth*norm[1], pos[2] + depth*norm[2]};
      dsSetColorAlpha (1, 1, 1, 1);
      dsDrawLine (contact[i].geom.pos, endp);
    }
  }
}

void ComanSimClass::applyExternalForce()
{
  if (!extforce.empty()) {
    for (auto it : extforce) {
      it.ApplyImpulse(sim_t);
    }
  }
}

void ComanSimClass::initExternalForce(const int& link_name, const Eigen::Vector3d& amp)
{
  // extforce.clear();
  if (link_name < 100) {
    for (auto it : oderobot) {
      Eigen::Vector3d amplitude(amp);
      ExternalForceProfile temp;
      dBodyID link_id = it->rlink[link_name].id;
      ODERobotBaseClass::RotZ(amplitude(0), amplitude(1), it->start_yaw);
      double duration = Dt_imp;
      // double duration = 0.114;
      // // double duration = 0.13;
      temp.init(link_id, "force", sim_t, amplitude, duration);
      extforce.push_back(temp);
    }
  }
  else {
    if (!objs.empty() && objs[0].id) {
      // for (auto it : oderobot) {
      for (auto it_obj : objs) {
        Eigen::Vector3d amplitude(amp);
        ExternalForceProfile temp;
        ODERobotBaseClass::RotZ(amplitude(0), amplitude(1), oderobot[0]->start_yaw);
        double duration = 0.3;
        temp.init(it_obj.id, "force", sim_t, amplitude, duration);
        extforce.push_back(temp);
      }
      // }
    }
  }
}

void ComanSimClass::initExternalTorque(const int& link_name, const Eigen::Vector3d& amp)
{
  // extforce.clear();
  if (link_name < 100) {
    for (auto it : oderobot) {
      Eigen::Vector3d amplitude(amp);
      ExternalForceProfile temp;
      dBodyID link_id = it->rlink[link_name].id;
      ODERobotBaseClass::RotZ(amplitude(0), amplitude(1), it->start_yaw);
      double duration = 0.13;
      temp.init(link_id, "torque", sim_t, amplitude, duration);
      extforce.push_back(temp);
    }
  }
  else {
    if (!objs.empty() && objs[0].id) {
      // for (auto it : oderobot) {
      for (auto it_obj : objs) {
        Eigen::Vector3d amplitude(amp);
        ExternalForceProfile temp;
        ODERobotBaseClass::RotZ(amplitude(0), amplitude(1), oderobot[0]->start_yaw);
        double duration = 0.3;
        temp.init(it_obj.id, "torque", sim_t, amplitude, duration);
        extforce.push_back(temp);
      }
      // }
    }
  }
}

void ComanSimClass::applyBallAttack()
{
  balls.clear();
  double sign_x = -sign(x_push);
  if (x_push == 0.0) sign_x = 0.0;
  double sign_y = -sign(y_push);
  if (y_push == 0.0) sign_y = 0.0;
  Eigen::Vector3d rel_pos(sign_x * 0.4, sign_y * 0.5, 0);
  RobotLink ball;
  ball.m = 1;
  ball.r = 0.1;
  ball.lx = -sign_x * 6; // velocity
  ball.ly = -sign_y * 6; // velocity
  ball.color_r = 1.0;
  ball.color_g = 0.8;
  ball.color_b = 0.0;
  for (auto it : oderobot)
  {
    ODERobotBaseClass::RotZ(rel_pos(0), rel_pos(1), it->start_yaw);
    BallAttack(it->rlink[TORSO].id, rel_pos, ball);
    balls.push_back(ball);
  }
}

void ComanSimClass::restartODE()
{
  // destroy
  dJointGroupDestroy(contactgroup);
  for (auto it : oderobot) {
    for (int i = 0; i < it->rlink.size(); i++) {
      dBodyDestroy(it->rlink[i].id);
      dGeomDestroy(it->rlink[i].gid);
    }
  }

// create again
  contactgroup = dJointGroupCreate(0);
  makeRobot();
//makeSlope();
}

void ComanSimClass::init()
{
  loop = 0;
  IsDropDownRobot = false;
  oderobot.clear();
  oderobot.push_back(boost::shared_ptr<ComanODERobotClass>(new ComanODERobotClass(0.55, -0, 0, 0)));

  COUT("===================== Start loading robot ====================");

  double z = RobotParaClass::FULL_LEG() + 0.02;
  // if (RobotParaClass::RobotName() == "coman") {
  if (oderobot[0]->name == "coman") {
    // oderobot.clear();
    // oderobot.push_back(boost::shared_ptr<ComanODERobotClass>(new ComanODERobotClass(z, -0, 0, 0)));
  }
  // else if (RobotParaClass::RobotName() == "bigman") {
  else if (oderobot[0]->name == "bigman") {
    oderobot.clear();
    oderobot.push_back(boost::shared_ptr<BigmanODERobotClass>(new BigmanODERobotClass(z, -0, 0, 0)));
  }
  else if (oderobot[0]->name == "cogimon") {
    oderobot.clear();
    oderobot.push_back(boost::shared_ptr<CogimonODERobotClass>(new CogimonODERobotClass(z, -0, 0, 0)));
  }

  dt = RobotParaClass::dT();//oderobot[0]->dt;
  g = RobotParaClass::G();//oderobot[0]->g;
  fixedhead.clear();
}

void ComanSimClass::makeRobot()
{
  init();
  if (!oderobot.empty()) {
    for (auto it : oderobot) {
      it->makeRobot(world, space);
      fixedhead.push_back(dJointCreateFixed(world, 0));
      dJointAttach(fixedhead.back(), it->rlink[PELVIS].id, 0);
      dJointSetFixed(fixedhead.back());
    }
  }
  makeToughTerrain();
}

void ComanSimClass::makeObj()
{
  if (!oderobot.empty()) {
    fixedjntID.clear();
    objs.clear();
    objfixedjnts.clear();
    RobotLink obj;
    obj.type = BOX;
    obj.m = 0.0001;
    obj.r = 0.1;
    obj.lx = 0.05;
    obj.ly = 0.5;
    obj.lz = 0.02;
    obj.color_r = 1.0;
    obj.color_g = 0.8;
    obj.color_b = 0.0;
    Eigen::Vector3d rel_pos(0.04 + 0.5 * obj.lx, 0, 0);
    for (auto it : oderobot) {
      ODERobotBaseClass::RotZ(rel_pos(0), rel_pos(1), it->start_yaw);
      const dReal *pos1 = dBodyGetPosition(it->rlink[LEFT_HAND].id);
      const dReal *pos2 = dBodyGetPosition(it->rlink[RIGHT_HAND].id);
      obj.px = 0.5 * (pos1[0] + pos2[0]) + rel_pos[0];
      obj.py = 0.5 * (pos1[1] + pos2[1]) + rel_pos[1];
      obj.pz = 0.5 * (pos1[2] + pos2[2]) + rel_pos[2];

      it->makeLink(obj, world, space);
      objs.push_back(obj);

      fixedjntID.push_back(dJointCreateFixed(world, 0));
      dJointAttach(fixedjntID.back(), obj.id, it->rlink[LEFT_HAND].id);
      dJointSetFixed(fixedjntID.back());
      dJointSetFeedback(fixedjntID.back(), &objs_jnt_fb[0]);

      RobotJoint jnt_l;
      jnt_l.id = fixedjntID.back();
      jnt_l.links[0] = &it->rlink[LEFT_HAND];
      jnt_l.links[1] = &objs.back();
      objfixedjnts.push_back(jnt_l);

      fixedjntID.push_back(dJointCreateFixed(world, 0));
      dJointAttach(fixedjntID.back(), obj.id, it->rlink[RIGHT_HAND].id);
      dJointSetFixed(fixedjntID.back());
      dJointSetFeedback(fixedjntID.back(), &objs_jnt_fb[1]);

      RobotJoint jnt_r;
      jnt_r.id = fixedjntID.back();
      jnt_r.links[0] = &it->rlink[RIGHT_HAND];
      jnt_r.links[1] = &objs.back();
      objfixedjnts.push_back(jnt_r);
    }
  }
}

///make steps and terrain
void ComanSimClass::makeObj2()
{
}

Eigen::Vector3d obs_size, obs_center_pos, obs_center_ori;
#define OBS_NUM 12
dGeomID obstacle[OBS_NUM];
void ComanSimClass::makeToughTerrain()
{
 
  // //////////////////////////////////obstacle;
  ////////////flat ground terrain configuration
  dReal lx = 0.2, ly = 0.5, lz = 0.02;
  dReal x = 10.3 + lx / 2, y = 0,  z = lz / 2;
//   dReal spacing_x = 0.0 + 1 * lx, spacing_y = 0, spacing_z = 0.1 + 0 * lz;
  
/// for 0.02m step walking
//   dReal lx = 0.25, ly = 0.5, lz = 0.02;
//   dReal x = 0.6 + lx / 2, y = 0,  z = lz / 2;
  
/// for 0.05m step walking  
//   dReal lx = 0.3, ly = 0.8, lz = 0.05;
//   dReal x = 0.66 + lx / 2, y = 0,  z = lz / 2;
  
/// for 0.06m step walking  
//  dReal lx = 0.3, ly = 0.5, lz = 0.06;
//  dReal x = 0.65 + lx / 2, y = 0,  z = lz / 2;  
  
/// for 0.07m step walking  
//  dReal lx = 0.3, ly = 0.5, lz = 0.07;
//  dReal x = 0.66 + lx / 2, y = 0,  z = lz / 2;    
  
/// for 0.08m step walking  
//  dReal lx = 0.3, ly = 0.5, lz = 0.08;
//  dReal x = 0.86 + lx / 2, y = 0,  z = lz / 2;   
 
/// for 0.09m step walking  
//  dReal lx = 0.3, ly = 0.5, lz = 0.09;
 // dReal x = 0.86 + lx / 2, y = 0,  z = lz / 2;     
 
/// for 0.10m step walking  
/* dReal lx = 0.3, ly = 0.5, lz = 0.1;
 dReal x = 1.17 + lx / 2, y = 0,  z = lz / 2;  */  
  
  dReal spacing_x = 0.0 + 1 * lx, spacing_y = 0, spacing_z = 0 + 1 * lz;  
  
  
  dReal ax = 0,   ay = 0, az = 1;   // 回転軸ベクトル;
  dReal angle = DEGTORAD(0);      // 回転角;
  dMatrix3 R;
  dRFromAxisAndAngle(R, ax, ay, az, angle);

  for (int i = 0; i < OBS_NUM; i++) {
//     obstacle[i]  = dCreateBox(space, lx, ly, lz); // 直方体ジオメトリの生成;
//     dGeomSetPosition(obstacle[i], x + i * spacing_x, y + i * spacing_y, z + i * spacing_z);  // 位置の設定;
//     dGeomSetRotation(obstacle[i], R);                // 姿勢の設定;
    
    
  
    if (i<5)
    {
    obstacle[i]  = dCreateBox(space, lx, ly, lz); // 直方体ジオメトリの生成;
    dGeomSetPosition(obstacle[i], x + i * spacing_x, y + i * spacing_y, z + i * spacing_z);
    dGeomSetRotation(obstacle[i], R);                
      
    }
    else
    {
/*      obstacle[i]  = dCreateBox(space, lx, ly, lz);
      dGeomSetPosition(obstacle[i], x + i * spacing_x, y + i * spacing_y, z + 4  * spacing_z);  
      dGeomSetRotation(obstacle[i], R); */  

      /// 5cm up&down  
//       if (i<7)
//       {
// 	obstacle[i]  = dCreateBox(space, 0.3, ly, lz);
// 	dGeomSetPosition(obstacle[i], x + i * spacing_x-0.1, y + i * spacing_y, z + 4  * spacing_z); 
// 	dGeomSetRotation(obstacle[i], R);                	
//       }
//       else
//       {
// 	obstacle[i]  = dCreateBox(space, 0.3, ly, lz);
// 	dGeomSetPosition(obstacle[i], x + i * spacing_x-0.1, y + i * spacing_y, z + 4  * spacing_z-(i-7)*spacing_z);  
// 	dGeomSetRotation(obstacle[i], R);                	
//       }

      ///// 10cm up&down  
      if (i<7)
      {
	obstacle[i]  = dCreateBox(space, lx, ly, lz);
	dGeomSetPosition(obstacle[i], x + i * spacing_x-0.1, y + i * spacing_y, z + 4  * spacing_z); 
	dGeomSetRotation(obstacle[i], R);                	
      }
      else
      {
	obstacle[i]  = dCreateBox(space, 0.3, ly, lz);
	dGeomSetPosition(obstacle[i], x + i * spacing_x-0.1, y + i * spacing_y, z + 4  * spacing_z-(i-7)*spacing_z);  
	dGeomSetRotation(obstacle[i], R);                	
      }
    
    }
    // dGeomSetRotation(obstacle[i], R[i]);               
  }


  const dReal *pos2 = dGeomGetPosition(obstacle[0]);
  COUT("obstacle position", pos2[0], pos2[1], pos2[2]);
}

void ComanSimClass::drawToughTerrain()
{
  const dReal *pos, *R;
  // dsSetTexture(DS_WOOD);
  // dsSetColor(0.8, 0.8, 0.8);
  dsSetColor(1.1, 1.1, 1.1);
  for (int i = 0; i < OBS_NUM; i++) {
    pos = dGeomGetPosition(obstacle[i]);
    R   = dGeomGetRotation(obstacle[i]);
    dVector3 sides;
    dGeomBoxGetLengths(obstacle[i], sides);
    dsDrawBox(pos, R, sides);
  }
  dsSetTexture(DS_NONE);
}
