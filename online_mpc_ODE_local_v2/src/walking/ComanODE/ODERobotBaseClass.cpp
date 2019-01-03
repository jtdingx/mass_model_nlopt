/*****************************************************************************
ODERobotBaseClass.cpp

Description:  cpp file of ODERobotBaseClass

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: 2016/03/18
@Update:  2016/03/20
*****************************************************************************/

#include "ComanODE/ODERobotBaseClass.h"

ODERobotBaseClass::ODERobotBaseClass(const double &z, const double &y, const double &x, const double &yaw)
{
  Init(RTControl.RobotPara(), z, y, x, yaw);
}

ODERobotBaseClass::ODERobotBaseClass(const RobotParaClass& robotpara, const double &z, const double &y, const double &x, const double &yaw)
{
  Init(robotpara, z, y, x, yaw);
};

void ODERobotBaseClass::Init(const RobotParaClass& robotpara, const double &z, const double &y, const double &x, const double &yaw)
{
  // RTControl.initOpenSoTIK();

  IsFinishHoming = false;
  JointNUM = RobotParaClass::JOINT_NUM();
  mComanWBS = &RTControl._WBS;
  start_x = x;
  start_y = y;
  start_z = z;
  start_yaw = yaw;
  name = robotpara.name;
  l_foot_FT = std::vector<double>(6, 0.0);
  l_hand_FT = std::vector<double>(6, 0.0);
  r_foot_FT = std::vector<double>(6, 0.0);
  r_hand_FT = std::vector<double>(6, 0.0);
  totalmass = 0.0;
  gcom << 0.0, 0.0, 0.0;
  gdcom << 0.0, 0.0, 0.0;
  gcom_old << 0.0, 0.0, 0.0;
  COMqueue = std::deque<Eigen::Vector3d>(1500, Eigen::Vector3d(x, y, z));
  dt = robotpara.dt;
  g = robotpara.g;
  totalmass = robotpara.totalmass;
  torque_diff_old.resize(JointNUM);
  rlink.resize(RobotParaClass::LINK_NUM());
  rjoint.resize(JointNUM);
  joint_angle_ref.resize(JointNUM);
  joint_torque_ref.resize(JointNUM);
  joint_torque_comp_ref.resize(JointNUM);
  joint_vel_ref.resize(JointNUM);
  jnt_fb.resize(JointNUM);

  torque_fb_filter.resize(JointNUM);
  // COUT("rlink size is ", rlink.size());
  // COUT("rjoint size is ", rjoint.size());
  readJointParam(robotpara);
  readLinkParam(robotpara);

#ifdef PUBLISH_TO_ROS
  data_publisher.reset(new ChartClass(robotpara.name));
  data_publisher->addPublisher("com", "Point");
  data_publisher->addPublisher("l_foot_FT", "Wrench");
  data_publisher->addPublisher("r_foot_FT", "Wrench");
  data_publisher->addPublisher("l_hand_FT", "Wrench");
  data_publisher->addPublisher("r_hand_FT", "Wrench");
#endif

  std::cout << "=============== Finish Initializing the ODERobotBaseClass ====================\n";
}

void ODERobotBaseClass::readLinkParam(const RobotParaClass& robotpara)
{
  for (int i = 0; i < rlink.size(); ++i)
  {
    rlink[i].type = CAPSULE;
    rlink[i].ds_type = CYLINDER;

    if (robotpara.getLink(i)) {
      rlink[i].IsExist = true;
      rlink[i].name = robotpara.getLink(i)->name;
      // COUT("link name is ", robotpara.getLink(i)->name);

      if (robotpara.getLink(i)->inertial) {
        rlink[i].m = robotpara.getLink(i)->inertial->mass;
        rlink[i].cgx = robotpara.getLink(i)->inertial->origin.position.x;
        rlink[i].cgy = robotpara.getLink(i)->inertial->origin.position.y;
        rlink[i].cgz = robotpara.getLink(i)->inertial->origin.position.z;
        rlink[i].I11 = robotpara.getLink(i)->inertial->ixx;
        rlink[i].I22 = robotpara.getLink(i)->inertial->iyy;
        rlink[i].I33 = robotpara.getLink(i)->inertial->izz;
        rlink[i].I12 = robotpara.getLink(i)->inertial->ixy;
        rlink[i].I13 = robotpara.getLink(i)->inertial->ixz;
        rlink[i].I23 = robotpara.getLink(i)->inertial->izz;
      }
      else {
        rlink[i].m = JOINT_MASS;
        rlink[i].I11 = JOINT_MASS;
        rlink[i].I22 = JOINT_MASS;
        rlink[i].I33 = JOINT_MASS;
      }

      COUT(robotpara.getLink(i)->name, rlink[i].m, " kg");

      if (robotpara.getLink(i)->visual) {
        rlink[i].color_r = robotpara.getLink(i)->visual->material->color.r;
        rlink[i].color_g = robotpara.getLink(i)->visual->material->color.g;
        rlink[i].color_b = robotpara.getLink(i)->visual->material->color.b;
        rlink[i].color_a = robotpara.getLink(i)->visual->material->color.a;
      }
      else {
        rlink[i].lx = JOINT_SIZE;
        rlink[i].ly = JOINT_SIZE;
        rlink[i].lz = JOINT_SIZE;
      }
    }
    rlink[i].r = 0.05;
    rlink[i].direction = 3; // z
  }


}

void ODERobotBaseClass::readJointParam(const RobotParaClass& robotpara)
{
  for (int i = 0; i < rjoint.size(); ++i)
  {
    if (robotpara.getJoint(i)) {
      rjoint[i].IsExist = true;
      rjoint[i].name = robotpara.getJoint(i)->name;
      rjoint[i].type = HINGE;
      if (robotpara.getJoint(i)->limits) {
        rjoint[i].fmax = robotpara.getJoint(i)->limits->effort;
        rjoint[i].vmax = robotpara.getJoint(i)->limits->velocity;
        rjoint[i].lo_stop = robotpara.getJoint(i)->limits->lower;
        rjoint[i].hi_stop = robotpara.getJoint(i)->limits->upper;
      }

      // if (robotpara.getJoint(i)->axis) {
      rjoint[i].axis_x = robotpara.getJoint(i)->axis.x;
      rjoint[i].axis_y = robotpara.getJoint(i)->axis.y;
      rjoint[i].axis_z = robotpara.getJoint(i)->axis.z;
      // }

      rjoint[i].link[0] = robotpara.getLinkName(robotpara.getJoint(i)->parent_link_name); // parent link
      rjoint[i].link[1] = robotpara.getLinkName(robotpara.getJoint(i)->child_link_name); // child link

      rjoint[i].links[0] = &rlink[rjoint[i].link[0]]; // parent link
      rjoint[i].links[1] = &rlink[rjoint[i].link[1]]; // child link

      rjoint[i].px = start_x + getJointPosition(robotpara, i).at(0);
      rjoint[i].py = start_y + getJointPosition(robotpara, i).at(1);
      rjoint[i].pz = start_z + getJointPosition(robotpara, i).at(2);

      rjoint[i].fudge_factor = 1; // It should have a value between zero and one (the default value). If the jumping motion is too visible in a joint, the value can be reduced. Making this value too small can prevent the motor from being able to move the joint away from a stop.
      rjoint[i].bounce = 0; // The bouncyness of the stops. This is a restitution parameter in the range 0..1. 0 means the stops are not bouncy at all, 1 means maximum bouncyness.
    }
  }
}

/**
 * @brief      { return joint position in global frame }
 *
 * @param[in]  robotpara   { parameter_description }
 * @param[in]  joint_name  { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
std::vector<double> ODERobotBaseClass::getJointPosition(const RobotParaClass& robotpara, const int& joint_name)
{
  std::vector<double> joint_position(3, 0.0);
  auto it = robotpara.getJoint(joint_name);

  // COUT(joint_name, it->name);

  while (it->parent_link_name != robotpara.getLinkName(PELVIS)) {
    // COUT("it->name", it->name);
    // COUT("it->parent_link_name", it->parent_link_name);
    joint_position[0] += it->parent_to_joint_origin_transform.position.x;
    joint_position[1] += it->parent_to_joint_origin_transform.position.y;
    joint_position[2] += it->parent_to_joint_origin_transform.position.z;
    it = robotpara.getLink(it->parent_link_name)->parent_joint;
    // IOClass::vout(joint_position, "joint_position");
  }

  // COUT("it->parent_link_name", it->parent_link_name);
  joint_position[0] += it->parent_to_joint_origin_transform.position.x;
  joint_position[1] += it->parent_to_joint_origin_transform.position.y;
  joint_position[2] += it->parent_to_joint_origin_transform.position.z;

  // joint_position[0] *= -1;
  // joint_position[1] *= -1;
  // IOClass::vout(joint_position, "joint_position");
  // COUT("");

  return joint_position;
}


void ODERobotBaseClass::makeRobot(const dWorldID &world, const dSpaceID &space)
{
  for (int i = 0; i < rlink.size(); ++i) {
    if (rlink[i].IsExist) {
      makeLink(rlink[i], world, space);
    }
  }

  // Set joint paramter
  for (int i = 0; i < rjoint.size(); ++i) {
    if (rjoint[i].IsExist) {
      rjoint[i].id = dJointCreateHinge(world, 0);
      dJointAttach(rjoint[i].id , rlink[rjoint[i].link[1]].id, rlink[rjoint[i].link[0]].id);
      dJointSetHingeAnchor(rjoint[i].id, rjoint[i].px, rjoint[i].py, rjoint[i].pz);
      dJointSetHingeAxis(rjoint[i].id, rjoint[i].axis_x, rjoint[i].axis_y, rjoint[i].axis_z);
      dJointSetHingeParam(rjoint[i].id, dParamLoStop, rjoint[i].lo_stop);
      dJointSetHingeParam(rjoint[i].id, dParamHiStop, rjoint[i].hi_stop);
      dJointSetHingeParam(rjoint[i].id, dParamFMax,   0 * rjoint[i].fmax);
      // dJointSetHingeParam(rjoint[i].id, dParamFudgeFactor, rjoint[i].fudge_factor);
      dJointSetFeedback(rjoint[i].id, &jnt_fb[i]);
    }
  }
}

void ODERobotBaseClass::makeLink(RobotLink& robotlink, const dWorldID &world, const dSpaceID &space)
{
  robotlink.id = dBodyCreate(world);
  dBodySetPosition(robotlink.id, robotlink.px, robotlink.py, robotlink.pz);

  dMass m;

  if (robotlink.type == SPHERE) {
    dQuaternion q;
    dQFromAxisAndAngle(q, 0, 0, 1, start_yaw);
    dBodySetQuaternion(robotlink.id, q);
    dMassSetSphereTotal(&m, robotlink.m, robotlink.r);
    dMassAdjust(&m, robotlink.m);
    dBodySetMass(robotlink.id, &m);
    robotlink.gid = dCreateSphere(space, robotlink.r);
  }
  else  if (robotlink.type == BOX) {
    dQuaternion q;
    dQFromAxisAndAngle(q, 0, 0, 1, start_yaw);
    dBodySetQuaternion(robotlink.id, q);
    dMassSetBoxTotal(&m, robotlink.m, robotlink.lx, robotlink.ly, robotlink.lz);
    dMassAdjust(&m, robotlink.m);
    dBodySetMass(robotlink.id, &m);
    robotlink.gid = dCreateBox(space, robotlink.lx, robotlink.ly, robotlink.lz);
  }
  else if (robotlink.type == CAPSULE) {
    std::vector<dReal> axis(3, 0.0);
    dReal capsule_length = 0.0;
    dReal angle = 0.5 * M_PI;
    if (robotlink.direction == 1) {
      capsule_length = robotlink.lx - 2 * robotlink.r;
      axis[1] = 1.0;
    }
    else if (robotlink.direction == 2) {
      capsule_length = robotlink.ly - 2 * robotlink.r;
      axis[0] = 1.0;
    }
    else if (robotlink.direction == 3) {
      capsule_length = robotlink.lz - 2 * robotlink.r;
      axis[2] = 1.0;
      angle = 0.0;
    }
    if (capsule_length < JOINT_SIZE) {
      capsule_length = JOINT_SIZE;
      robotlink.r = 0.1 * JOINT_SIZE;
    }
    dQuaternion qc; // local orientation
    dQFromAxisAndAngle(qc, axis[0], axis[1], axis[2], angle);
    dQuaternion qb; // start_yaw
    dQFromAxisAndAngle(qb, 0, 0, 1, start_yaw);
    dQuaternion q;
    dQMultiply0(q, qb, qc);
    dBodySetQuaternion(robotlink.id, q);
    dMassSetCapsuleTotal(&m, robotlink.m, robotlink.direction, robotlink.r, capsule_length);
    dMassAdjust(&m, robotlink.m);
    dBodySetMass(robotlink.id, &m);
    robotlink.gid = dCreateCapsule(space, robotlink.r, capsule_length);
  }
  else  if (robotlink.type == CYLINDER) {
    std::vector<dReal> axis(3, 0.0);
    dReal cylinder_length = 0.0;
    dReal angle = 0.5 * M_PI;
    if (robotlink.direction == 1) {
      cylinder_length = robotlink.lx;
      axis[1] = 1.0;
    }
    else if (robotlink.direction == 2) {
      cylinder_length = robotlink.ly;
      axis[0] = 1.0;
    }
    else if (robotlink.direction == 3) {
      cylinder_length = robotlink.lz;
      axis[2] = 1.0;
      angle = 0.0;
    }
    dQuaternion qc; // local orientation
    dQFromAxisAndAngle(qc, axis[0], axis[1], axis[2], angle);
    dQuaternion qb; // start_yaw
    dQFromAxisAndAngle(qb, 0, 0, 1, start_yaw);
    dQuaternion q;
    dQMultiply0(q, qb, qc);
    dBodySetQuaternion(robotlink.id, q);
    dMassSetCylinderTotal(&m, robotlink.m, robotlink.direction, robotlink.r, cylinder_length);
    dMassAdjust(&m, robotlink.m);
    dBodySetMass(robotlink.id, &m);
    robotlink.gid = dCreateCylinder(space, robotlink.r, cylinder_length);
  }

  dGeomSetBody(robotlink.gid, robotlink.id);
}


void ODERobotBaseClass::drawLink(const RobotLink& robotlink)
{
  dReal sides[3], sides_torso1[3];
  const dReal *sides_torso2;

  sides[0] = robotlink.lx;
  sides[1] = robotlink.ly;
  sides[2] = robotlink.lz;

  dsSetColor(robotlink.color_r, robotlink.color_g, robotlink.color_b);

  if (robotlink.ds_type == SPHERE) {
    dsDrawSphere(dBodyGetPosition(robotlink.id), dBodyGetRotation(robotlink.id), robotlink.r);
  }
  else if (robotlink.ds_type == BOX) {
    dsDrawBox(dBodyGetPosition(robotlink.id), dBodyGetRotation(robotlink.id), sides);
  }
  else if (robotlink.ds_type == CAPSULE) {
    dReal capsule_length = - 2 * robotlink.r;
    if (robotlink.direction == 1)  capsule_length += robotlink.lx;
    if (robotlink.direction == 2)  capsule_length += robotlink.ly;
    if (robotlink.direction == 3)  capsule_length += robotlink.lz;
    if (capsule_length < JOINT_SIZE) {
      capsule_length = JOINT_SIZE;
    }
    dsDrawCapsule(dBodyGetPosition(robotlink.id), dBodyGetRotation(robotlink.id), capsule_length, robotlink.r);
  }
  else if (robotlink.ds_type == CYLINDER) {
    double cylinder_length = 0.0;
    if (robotlink.direction == 1)  cylinder_length = robotlink.lx;
    if (robotlink.direction == 2)  cylinder_length = robotlink.ly;
    if (robotlink.direction == 3)  cylinder_length = robotlink.lz;
    dsDrawCylinder(dBodyGetPosition(robotlink.id), dBodyGetRotation(robotlink.id), cylinder_length, robotlink.r);
  }
}

void ODERobotBaseClass::controlMotors(std::string method)
{
  if (method != "position" && method != "torque" && method != "impedance" && method != "velocity" && method != "admittance") {
    std::cerr << "Error!!! controlMotors(std::string method) method is not correctly defined. Could be \"position\" or \"impedance\" or \"torque\" or \"velocity\" or \"admittance\". Now runing in default position control mode." << std::endl;
    method = "position";
  }
  for (int i = 0; i < rjoint.size(); i++) {
    controlMotor(method, i);
  }
}

void ODERobotBaseClass::controlMotor(std::string method, const int& joint_name)
{
  if (method != "position" && method != "torque" && method != "impedance" && method != "velocity" && method != "admittance") {
    std::cerr << "Error!!! controlMotor(std::string method) method is not correctly defined. Could be \"position\" or \"impedance\" or \"torque\" or \"velocity\" or \"admittance\". Now runing in default position control mode." << std::endl;
    method = "position";
  }

  if (rjoint[joint_name].IsExist) {
    if (method == "position") {
      PositionMotorControl(joint_name, joint_angle_ref[joint_name]);
    }
    else if (method == "impedance") {
      ImpedanceMotorControl(joint_name, joint_angle_ref[joint_name], joint_vel_ref[joint_name], joint_torque_ref[joint_name]);
    }
    else if (method == "torque") {
      TorqueMotorControl(joint_name, joint_torque_ref[joint_name]);
    }
    else if (method == "velocity") {
      VelocityMotorControl(joint_name, joint_vel_ref[joint_name]);
    }
    else if (method == "admittance") {
      AdmittanceMotorControl(joint_name, joint_angle_ref[joint_name]);
    }
  }
}

void ODERobotBaseClass::PositionMotorControl(const int& joint_name, dReal joint_angle)
{
  dReal kp = RTControl.RobotPara().kp_pos, kd = RTControl.RobotPara().kd_pos, ki = 2, u, diff;

  
  clamp(joint_angle, rjoint[joint_name].lo_stop, rjoint[joint_name].hi_stop);
  diff = joint_angle - dJointGetHingeAngle(rjoint[joint_name].id);
  u = kp * diff - kd * dJointGetHingeAngleRate(rjoint[joint_name].id);// + ki*idiff;

  dJointSetHingeParam(rjoint[joint_name].id, dParamVel, u);
  dJointSetHingeParam(rjoint[joint_name].id, dParamFMax, rjoint[joint_name].fmax);
  dJointSetHingeParam(rjoint[joint_name].id, dParamFudgeFactor, 0.1);
}

void ODERobotBaseClass::AdmittanceMotorControl(const int& joint_name, dReal joint_angle)
{
  dReal kp = RTControl.RobotPara().kp_adm, kd = RTControl.RobotPara().kd_adm, diff, u;

  double torque_fb;
  fb = dJointGetFeedback(rjoint[joint_name].id);
  if (rjoint[joint_name].axis_x) {
    torque_fb = fb->t2[0];
  }
  if (rjoint[joint_name].axis_y) {
    torque_fb = fb->t2[1];
  }
  if (rjoint[joint_name].axis_z) {
    torque_fb = fb->t2[2];
  }

  double FilterCutOff = 6;
  double N_ButWth = 4;
  diff = joint_torque_comp_ref[joint_name] - torque_fb;
  double A = - diff;
  u = (dt / (kp * dt + kd)) * A + kd / (kp * dt + kd) * torque_diff_old[joint_name] ;

  torque_diff_old[joint_name] = u;

  clamp(u, -DEGTORAD(15), DEGTORAD(15));
  PositionMotorControl(joint_name, joint_angle + u);
}

void ODERobotBaseClass::TorqueMotorControl(const int& joint_name, const dReal& torque)
{
  // the gravity compensation part is necessary for impedance mode
  dReal u = 0 * joint_torque_comp_ref[joint_name] + torque ;

  if (dJointGetHingeAngleRate(rjoint[joint_name].id) >= rjoint[joint_name].vmax) {
    u = u > 0 ? 0 : u;
  }
  else if (dJointGetHingeAngleRate(rjoint[joint_name].id) <= -rjoint[joint_name].vmax) {
    u = u < 0 ? 0 : u;
  }

  clamp(u, -(rjoint[joint_name].fmax), rjoint[joint_name].fmax);

  dJointAddHingeTorque(rjoint[joint_name].id, u);// for torque control
  // Set the joint friction
  dJointSetHingeParam(rjoint[joint_name].id, dParamVel, 0.0);
  dJointSetHingeParam(rjoint[joint_name].id, dParamFMax, 0.1); // torque friction no less than 0.7

}

void ODERobotBaseClass::ImpedanceMotorControl(const int& joint_name, dReal joint_angle, dReal joint_vel, const dReal& tau_ff)//target in radian
{
  dReal kp = RTControl.RobotPara().kp_imp, kd = RTControl.RobotPara().kd_imp, u, diff, vel_err;

  clamp(joint_angle, rjoint[joint_name].lo_stop, rjoint[joint_name].hi_stop);
  clamp(joint_vel, -rjoint[joint_name].vmax, rjoint[joint_name].vmax);
  diff = joint_angle - dJointGetHingeAngle(rjoint[joint_name].id);
  vel_err = joint_vel - dJointGetHingeAngleRate(rjoint[joint_name].id);
  u    = tau_ff + kp * diff + kd * vel_err;

  TorqueMotorControl(joint_name, u);
}

void ODERobotBaseClass::VelocityMotorControl(const int& joint_name, const dReal& joint_angular_velocity)
{

  dJointSetHingeParam(rjoint[joint_name].id, dParamVel, joint_angular_velocity);
  dJointSetHingeParam(rjoint[joint_name].id, dParamFMax, rjoint[joint_name].fmax);
}

void ODERobotBaseClass::RotZ(double& px, double& py, const double& yaw)
{
  double x = px, y = py;
  px = std::cos(yaw) * x - std::sin(yaw) * y;
  py = std::sin(yaw) * x + std::cos(yaw) * y;
}

void ODERobotBaseClass::updateRobotYaw()
{
  for (int i = 0; i < rlink.size(); ++i) {
    rlink[i].px -= start_x;
    rlink[i].py -= start_y;
    RotZ(rlink[i].px, rlink[i].py, start_yaw);
    rlink[i].px += start_x;
    rlink[i].py += start_y;
  }

  for (int i = 0; i < rjoint.size(); ++i) {
    rjoint[i].px -= start_x;
    rjoint[i].py -= start_y;
    RotZ(rjoint[i].px, rjoint[i].py, start_yaw);
    rjoint[i].px += start_x;
    rjoint[i].py += start_y;
    RotZ(rjoint[i].axis_x, rjoint[i].axis_y, start_yaw);
  }
}

std::vector<double> ODERobotBaseClass::getLocalJointForceTorque(const RobotJoint& robotjoint)
{
  std::vector<double> ft;
  // in Foot frame
  Eigen::Vector3d force(0, 0, 0), torque(0, 0, 0);
  Eigen::Matrix3d R;

  // const dReal *rot = dBodyGetRotation(rlink[robotjoint.link[1]].id);
  const dReal *rot = dBodyGetRotation(robotjoint.links[1]->id);
  R << rot[0], rot[1], rot[2],
  rot[4], rot[5], rot[6],
  rot[8], rot[9], rot[10]; // ポインタ的表記では　r11 = *(rot+0); r12 = *(rot+1); r13 = *(rot+2)

  dJointFeedback *fb = dJointGetFeedback(robotjoint.id);
  force(0) = fb->f2[0];
  force(1) = fb->f2[1];
  force(2) = fb->f2[2];
  torque(0) = fb->t2[0];
  torque(1) = fb->t2[1];
  torque(2) = fb->t2[2];

// the joint feedback is in global frame therefore I transfer it to the local frame
  force = R.transpose() * force;
  torque = R.transpose() * torque;

  ft.push_back(force(0) );
  ft.push_back(force(1) );
  ft.push_back(force(2) );
  ft.push_back(torque(0));
  ft.push_back(torque(1));
  ft.push_back(torque(2));

  return ft;
}

void ODERobotBaseClass::ReadForeTorqueSensor()
{
  l_foot_FT = getLocalJointForceTorque(rjoint[LEFT_FOOT_PITCH]);
  r_foot_FT = getLocalJointForceTorque(rjoint[RIGHT_FOOT_PITCH]);
  if (!RTControl.RobotPara().IsOnlyLowerBody) {
    l_hand_FT = getLocalJointForceTorque(rjoint[LEFT_WRIST_ROLL]);
    r_hand_FT = getLocalJointForceTorque(rjoint[RIGHT_WRIST_ROLL]);
  }

}

void ODERobotBaseClass::CalcGlobalCOM()
{
  Eigen::Vector3d TotalCOM(0.0, 0.0, 0.0);

  const dReal *pos;
  for (auto it : rlink) {
    if (it.m > 0.0001) {
      pos = dBodyGetPosition(it.id);
      TotalCOM[0] += it.m * pos[0];
      TotalCOM[1] += it.m * pos[1];
      TotalCOM[2] += it.m * pos[2];
    }
  };

  gcom = TotalCOM / totalmass;
  gdcom = (gcom - gcom_old) / dt;
  gcom_old = gcom;

  COMqueue.push_back(gcom);
  COMqueue.pop_front();
}

void ODERobotBaseClass::PublishToRos()
{
#ifdef PUBLISH_TO_ROS
  data_publisher->send("com", gcom);
  data_publisher->send("l_foot_FT", l_foot_FT);
  data_publisher->send("r_foot_FT", r_foot_FT);
  data_publisher->send("l_hand_FT", l_hand_FT);
  data_publisher->send("r_hand_FT", r_hand_FT);
#endif
}

void ODERobotBaseClass::UpdateWBS()
{

  Eigen::Vector3d Ghip, pelvisOFF;
  Eigen::Vector3d Glft, Grft, Gcom;
  Eigen::Matrix3d Rpelvis;

  const dReal *rot = dBodyGetRotation(rlink[PELVIS].id);
  Rpelvis << rot[0], rot[1], rot[2],
          rot[4], rot[5], rot[6],
          rot[8], rot[9], rot[10]; // ポインタ的表記では　r11 = *(rot+0); r12 = *(rot+1); r13 = *(rot+2)

  const dReal *pos;
  pos = dBodyGetPosition(rlink[LEFT_FOOT_DUMMY].id);
  Glft(0) = pos[0];
  Glft(1) = pos[1];
  Glft(2) = pos[2];
  pos = dBodyGetPosition(rlink[RIGHT_FOOT_DUMMY].id);
  Grft(0) = pos[0];
  Grft(1) = pos[1];
  Grft(2) = pos[2];
  pos = dBodyGetPosition(rlink[PELVIS].id);
  Ghip(0) = pos[0];
  Ghip(1) = pos[1];
  Ghip(2) = pos[2];
  Gcom(0) = gcom[0];
  Gcom(1) = gcom[1];
  Gcom(2) = gcom[2];
  const dReal *linear_vel_raw, *angular_vel_raw;
  Eigen::Vector3d linear_accel(0, 0, 0), angular_accel(0, 0, 0);
  Eigen::Vector3d linear_vel(0, 0, 0), angular_vel(0, 0, 0);
  static std::vector<double> linear_vel_old(3, 0), angular_vel_old(3, 0);
  linear_vel_raw  = dBodyGetLinearVel(rlink[PELVIS].id);
  angular_vel_raw = dBodyGetAngularVel(rlink[PELVIS].id);

  linear_vel[0] = (linear_vel_raw[0]);
  linear_vel[1] = (linear_vel_raw[1]);
  linear_vel[2] = (linear_vel_raw[2]);

  angular_vel[0] = (angular_vel_raw[0]);
  angular_vel[1] = (angular_vel_raw[1]);
  angular_vel[2] = (angular_vel_raw[2]);

  angular_vel = Rpelvis.transpose() * angular_vel; // need to check

  for (int i = 0; i < 3; i++) {
    linear_accel[i] = (linear_vel[i] - linear_vel_old[i])  / dt;
    angular_accel[i] = (angular_vel[i] - angular_vel_old[i]) / dt;

    linear_vel_old[i] = linear_vel[i];
    angular_vel_old[i] = angular_vel[i];
  }
  linear_accel = Rpelvis.transpose() * linear_accel;

  // dReal qall[JOINT_NUM];
  std::vector<double> qall(JointNUM, 0);
  std::vector<double> tall(JointNUM, 0);
  for (int i = 0; i < JointNUM; i++) {
    if (rjoint[i].IsExist) {
      qall[i] = dJointGetHingeAngle(rjoint[i].id);
      double torque_fb = 0.0;
      fb = dJointGetFeedback(rjoint[i].id);
      if (rjoint[i].axis_x) {
        torque_fb = fb->t2[0];
      }
      if (rjoint[i].axis_y) {
        torque_fb = fb->t2[1];
      }
      if (rjoint[i].axis_z) {
        torque_fb = fb->t2[2];
      }
      tall[i] = torque_fb;
    }
  }

  Eigen::Vector3d EulerAng(0, 0, 0);
  EulerAng[0] = atan2(Rpelvis(2, 1), Rpelvis(2, 2));   //phi
  EulerAng[1] = asin(-Rpelvis(2, 0));           //theta
  EulerAng[2] = atan2(Rpelvis(1, 0), Rpelvis(0, 0));    //greek Y

  std::vector<double> FTSensor = r_foot_FT;
  FTSensor.insert(FTSensor.end(), l_foot_FT.begin(), l_foot_FT.end());
  std::vector<double> HandFTSensor = r_hand_FT;
  HandFTSensor.insert(HandFTSensor.end(), l_hand_FT.begin(), l_hand_FT.end());
  //******************** implementation of WBS class ************************
  const RobotStateClass& irobot = mComanWBS->getRobotState();
  RTControl.UpdateWBS(dt, qall, EulerAng, linear_accel, angular_vel, FTSensor, HandFTSensor);
  RTControl.UpdateJointTorqueFB(tall);
  //**************************************************************************

  //-------- Put your own controller here ------------------
  std::vector<float> OffsetAng(JointNUM, 0);
  RTControl.Run(dt, joint_angle_ref, OffsetAng);

  // send to the robot
  for (int i = 0; i < JointNUM; ++i) {
    joint_angle_ref[i] = irobot.SendToRobot->q_ref[i];
    joint_torque_ref[i] = irobot.SendToRobot->jtau[i];
    joint_vel_ref[i] = irobot.SendToRobot->qdot[i];
  }

}

void ODERobotBaseClass::UpdateSensorFeedback()
{
  ReadForeTorqueSensor();
  CalcGlobalCOM();
  PublishToRos();

  if (IsFinishHoming) {
    UpdateWBS();
  }
}

void ODERobotBaseClass::drawCOPLine()
{
  const RobotStateClass& irobotDS = mComanWBS->getRobotState();
  dVector3 pos1, pos2;
  dReal k1 = 0.0007;
  pos1[0] = irobotDS.gcop(0) + start_x;
  pos1[1] = irobotDS.gcop(1) + start_y;
  pos1[2] = 0; //irobotDS.gcop(2);
  pos2[0] = pos1[0];
  pos2[1] = pos1[1];
  pos2[2] = 0.1 + pos1[2];

  dsSetColor(0.2, 1, 0.2);
  dsDrawLine(pos1, pos2);
  dsDrawSphere(pos1, dBodyGetRotation(rlink[WAIST_YAW].id), 0.008);
}
