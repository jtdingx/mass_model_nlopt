/*****************************************************************************
ComanODERobotClass.cpp

Description:  cpp file of ComanODERobotClass

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: 2016/03/18
@Update:  2016/03/20
*****************************************************************************/

#include "ComanODE/ComanODERobotClass.h"

ComanODERobotClass::ComanODERobotClass(const double &z, const double &y, const double &x, const double &yaw)
  : ODERobotBaseClass(z, y, x, yaw)
{
  std::cout << "Start to construct Coman robot \n";
  updateJointParam(RTControl.RobotPara());
  updateLinkParam(RTControl.RobotPara());
  updateRobotYaw();

  std::cout << "No." << RobotParaClass::ROBOT_COUNT() << " Coman is ready..." << std::endl;
};

ComanODERobotClass::~ComanODERobotClass()
{

};

/**
 * @brief      used to update link positions and dimensions
 *
 * @param[in]  robotpara  { parameter_description }
 */
void ComanODERobotClass::updateLinkParam(const RobotParaClass& robotpara)
{
  for (int i = 0; i < rlink.size(); ++i) {
    if (robotpara.getLink(i)) {
      bool HasChild = !robotpara.getLink(i)->child_joints.empty();
      if (HasChild)
      {
        if (i == PELVIS) { // PELVIS has three children
          int child = WAIST_ROLL;
          rlink[i].ds_type = CAPSULE;
          rlink[i].lx = std::abs(rjoint[child].px - start_x);
          rlink[i].ly = std::abs(rjoint[child].py - start_y);
          rlink[i].lz = std::abs(rjoint[child].pz - start_z) + 0.05; // + offset for visualization
          rlink[i].px = rjoint[child].px - 0.5 * rlink[i].lx;
          rlink[i].py = rjoint[child].py - 0.5 * rlink[i].ly;
          rlink[i].pz = rjoint[child].pz - 0.5 * rlink[i].lz;
        }
        else if (i != RIGHT_FOOT && i != LEFT_FOOT && i != RIGHT_HAND && i != LEFT_HAND) {
          int parent = robotpara.getJointName(robotpara.getLink(i)->parent_joint->name);
          int child = robotpara.getJointName(robotpara.getLink(i)->child_joints[0]->name);
          if (i == TORSO) { // TORSO actually has three children: NECK, RIGHT_UPPER_ARM_DUMMY1, LEFT_UPPER_ARM_DUMMY1
            child = HEAD_PITCH;
          }

          rlink[i].px = rjoint[parent].px + 0.5 * (rjoint[child].px - rjoint[parent].px);
          rlink[i].py = rjoint[parent].py + 0.5 * (rjoint[child].py - rjoint[parent].py);
          rlink[i].pz = rjoint[parent].pz + 0.5 * (rjoint[child].pz - rjoint[parent].pz);

          rlink[i].lx = std::abs(rjoint[child].px - rjoint[parent].px);
          rlink[i].ly = std::abs(rjoint[child].py - rjoint[parent].py);
          rlink[i].lz = std::abs(rjoint[child].pz - rjoint[parent].pz);
        }

        if ( (rlink[i].lx + rlink[i].ly + rlink[i].lz) < JOINT_SIZE) {
          rlink[i].r = 0.1 * JOINT_SIZE;
          rlink[i].lx = JOINT_SIZE;
          rlink[i].ly = JOINT_SIZE;
          rlink[i].lz = JOINT_SIZE;
        }

        // COUT(rjoint[parent].name);
        // COUT(rjoint[parent].px, rjoint[parent].py, rjoint[parent].pz);

        // COUT(rlink[i].name, rlink[i].m);
        // COUT(rlink[i].px, rlink[i].py, rlink[i].pz);
        // COUT(rlink[i].lx, rlink[i].ly, rlink[i].lz);

        // COUT(rjoint[child].name);
        // COUT(rjoint[child].px, rjoint[child].py, rjoint[child].pz);
        // COUT(" ");
      }
    }
  }


  rlink[TORSO].m = robotpara.TORSO_MASS;
  rlink[HEAD].m = robotpara.HEAD_MASS;
  rlink[HEAD].type = SPHERE;
  rlink[HEAD].ds_type = SPHERE;
  rlink[HEAD].r = JOINT_SIZE * HEAD_RADIUS;
  rlink[HEAD].px = rlink[NECK].px;
  rlink[HEAD].py = rlink[NECK].py;
  rlink[HEAD].pz = rlink[NECK].pz + rlink[HEAD].r;
  rlink[HEAD].lx = HEAD_RADIUS;
  rlink[HEAD].ly = HEAD_RADIUS;
  rlink[HEAD].lz = HEAD_RADIUS;
  rlink[HEAD].color_r = 1;
  rlink[HEAD].color_g = 1;
  rlink[HEAD].color_b = 1;

  rlink[RIGHT_THIGH_DUMMY1].type = CYLINDER;
  rlink[RIGHT_THIGH_DUMMY1].direction = 2; // y

  rlink[LEFT_THIGH_DUMMY1].type = CYLINDER;
  rlink[LEFT_THIGH_DUMMY1].direction = 2; // y

  rlink[RIGHT_FOOT].type = BOX;
  rlink[RIGHT_FOOT].ds_type = BOX;
  rlink[RIGHT_FOOT].lx = RobotParaClass::FOOT_LENGTH();
  rlink[RIGHT_FOOT].ly = RobotParaClass::FOOT_WIDTH();
  rlink[RIGHT_FOOT].lz = RobotParaClass::FOOT_HEIGHT();
  rlink[RIGHT_FOOT].px = rlink[RIGHT_FOOT_DUMMY].px + RobotParaClass::ANKLE_X_OFFSET();
  rlink[RIGHT_FOOT].py = rlink[RIGHT_FOOT_DUMMY].py;
  rlink[RIGHT_FOOT].pz = rlink[RIGHT_FOOT_DUMMY].pz + 0.5 * (rlink[RIGHT_FOOT].lz) - RobotParaClass::ANKLE_HEIGHT();

  rlink[LEFT_FOOT].type = BOX;
  rlink[LEFT_FOOT].ds_type = BOX;
  rlink[LEFT_FOOT].lx = RobotParaClass::FOOT_LENGTH();
  rlink[LEFT_FOOT].ly = RobotParaClass::FOOT_WIDTH();
  rlink[LEFT_FOOT].lz = RobotParaClass::FOOT_HEIGHT();
  rlink[LEFT_FOOT].px = rlink[LEFT_FOOT_DUMMY].px + RobotParaClass::ANKLE_X_OFFSET();
  rlink[LEFT_FOOT].py = rlink[LEFT_FOOT_DUMMY].py;
  rlink[LEFT_FOOT].pz = rlink[LEFT_FOOT_DUMMY].pz +  0.5 * (rlink[LEFT_FOOT].lz) - RobotParaClass::ANKLE_HEIGHT();

  rlink[WAIST_DUMMY2].ds_type = CAPSULE;
  rlink[WAIST_DUMMY2].direction = 2;
  rlink[WAIST_DUMMY2].r = 0.045;
  rlink[WAIST_DUMMY2].ly = 0.13;
  // rlink[WAIST_DUMMY2].color_r = 0.1; //1.2 * rlink[TORSO].color_r;
  // rlink[WAIST_DUMMY2].color_g = 0.1; //1.2 * rlink[TORSO].color_g;
  // rlink[WAIST_DUMMY2].color_b = 0.1; //1.2 * rlink[TORSO].color_b;

  rlink[TORSO].ds_type = CAPSULE;
  rlink[TORSO].r = 0.1;

  rlink[RIGHT_UPPER_ARM_DUMMY1].type = CYLINDER;
  rlink[RIGHT_UPPER_ARM_DUMMY1].direction = 2; // y
  rlink[RIGHT_UPPER_ARM_DUMMY2].type = CYLINDER;
  rlink[RIGHT_UPPER_ARM_DUMMY2].r = 0.04;
  rlink[RIGHT_UPPER_ARM].ds_type = CAPSULE;
  rlink[RIGHT_ELBOW_FORE_ARM].type = CAPSULE;
  rlink[RIGHT_ELBOW_FORE_ARM].r = 0.04;
  rlink[RIGHT_FORE_ARM].r = 0.04;

  rlink[LEFT_UPPER_ARM_DUMMY1].type = CYLINDER;
  rlink[LEFT_UPPER_ARM_DUMMY1].direction = 2; // y
  rlink[LEFT_UPPER_ARM_DUMMY2].type = CYLINDER;
  rlink[LEFT_UPPER_ARM_DUMMY2].r = 0.04;
  rlink[LEFT_UPPER_ARM].ds_type = CAPSULE;
  rlink[LEFT_ELBOW_FORE_ARM].type = CAPSULE;
  rlink[LEFT_ELBOW_FORE_ARM].r = 0.04;
  rlink[LEFT_FORE_ARM].r = 0.04;

  // ==== box hands
  rlink[RIGHT_HAND].type = BOX;
  rlink[RIGHT_HAND].ds_type = BOX;
  rlink[RIGHT_HAND].lx = HAND_LENGTH;
  rlink[RIGHT_HAND].ly = HAND_WIDTH;
  rlink[RIGHT_HAND].lz = HAND_HEIGHT;
  rlink[RIGHT_HAND].px = rlink[RIGHT_HAND_DUMMY1].px;
  rlink[RIGHT_HAND].py = rlink[RIGHT_HAND_DUMMY1].py;
  rlink[RIGHT_HAND].pz = rlink[RIGHT_HAND_DUMMY1].pz - 0.5 * rlink[RIGHT_HAND].lz;

  rlink[LEFT_HAND].type = BOX;
  rlink[LEFT_HAND].ds_type = BOX;
  rlink[LEFT_HAND].lx = HAND_LENGTH;
  rlink[LEFT_HAND].ly = HAND_WIDTH;
  rlink[LEFT_HAND].lz = HAND_HEIGHT;
  rlink[LEFT_HAND].px = rlink[LEFT_HAND_DUMMY1].px;
  rlink[LEFT_HAND].py = rlink[LEFT_HAND_DUMMY1].py;
  rlink[LEFT_HAND].pz = rlink[LEFT_HAND_DUMMY1].pz - 0.5 * rlink[LEFT_HAND].lz;

  // === sphere hands
  // rlink[RIGHT_HAND].type = SPHERE;
  // rlink[RIGHT_HAND].ds_type = SPHERE;
  // rlink[RIGHT_HAND].r = 0.03;
  rlink[RIGHT_HAND].color_r = 1.0;
  rlink[RIGHT_HAND].color_g = 1.0;
  rlink[RIGHT_HAND].color_b = 1.0;
  // rlink[LEFT_HAND].type = SPHERE;
  // rlink[LEFT_HAND].ds_type = SPHERE;
  // rlink[LEFT_HAND].r = 0.03;
  rlink[LEFT_HAND].color_r = 1.0;
  rlink[LEFT_HAND].color_g = 1.0;
  rlink[LEFT_HAND].color_b = 1.0;
}

void ComanODERobotClass::updateJointParam(const RobotParaClass& robotpara)
{
  // for neck and head pitch joints
  rjoint[NECK_PITCH].fmax    = rjoint[LEFT_WRIST_PITCH].fmax    ;
  rjoint[NECK_PITCH].vmax    = rjoint[LEFT_WRIST_PITCH].vmax    ;
  rjoint[NECK_PITCH].lo_stop = rjoint[LEFT_WRIST_PITCH].lo_stop ;
  rjoint[NECK_PITCH].hi_stop = rjoint[LEFT_WRIST_PITCH].hi_stop ;
  rjoint[NECK_PITCH].axis_x  = rjoint[LEFT_WRIST_PITCH].axis_x  ;
  rjoint[NECK_PITCH].axis_y  = rjoint[LEFT_WRIST_PITCH].axis_y  ;
  rjoint[NECK_PITCH].axis_z  = rjoint[LEFT_WRIST_PITCH].axis_z  ;

  rjoint[HEAD_PITCH].fmax    = rjoint[LEFT_WRIST_PITCH].fmax    ;
  rjoint[HEAD_PITCH].vmax    = rjoint[LEFT_WRIST_PITCH].vmax    ;
  rjoint[HEAD_PITCH].lo_stop = rjoint[LEFT_WRIST_PITCH].lo_stop ;
  rjoint[HEAD_PITCH].hi_stop = rjoint[LEFT_WRIST_PITCH].hi_stop ;
  rjoint[HEAD_PITCH].axis_x  = rjoint[LEFT_WRIST_PITCH].axis_x  ;
  rjoint[HEAD_PITCH].axis_y  = rjoint[LEFT_WRIST_PITCH].axis_y  ;
  rjoint[HEAD_PITCH].axis_z  = rjoint[LEFT_WRIST_PITCH].axis_z  ;

  rjoint[RIGHT_SHOULDER_PITCH].lo_stop = DEGTORAD(-170);
  rjoint[LEFT_SHOULDER_PITCH].lo_stop = DEGTORAD(-170);

  rjoint[RIGHT_SHOULDER_ROLL].hi_stop = DEGTORAD(5);
  rjoint[RIGHT_SHOULDER_ROLL].lo_stop = DEGTORAD(-180);
  rjoint[LEFT_SHOULDER_ROLL].hi_stop = DEGTORAD(180);
  rjoint[LEFT_SHOULDER_ROLL].lo_stop = DEGTORAD(-5);

  // rjoint[RIGHT_ELBOW_PITCH].hi_stop = DEGTORAD(120);
  // rjoint[LEFT_ELBOW_PITCH].hi_stop = DEGTORAD(120);
}


void ComanODERobotClass::drawRobot()
{
  dReal sides[3], sides_torso1[3];
  const dReal *sides_torso2;

  for (int i = 0; i < rlink.size(); ++i) {
    if (rlink[i].IsExist) {
      if (i == RIGHT_ELBOW_FORE_ARM || i == RIGHT_FORE_ARM) {
        const dReal *pos1, *pos2;
        dReal pos[3] = {0.0, 0.0, 0.0};
        pos1 = dBodyGetPosition(rlink[RIGHT_ELBOW_FORE_ARM].id);
        pos2 = dBodyGetPosition(rlink[RIGHT_FORE_ARM].id);
        pos[0] = 0.5 * (pos1[0] + pos2[0]);
        pos[1] = 0.5 * (pos1[1] + pos2[1]);
        pos[2] = 0.5 * (pos1[2] + pos2[2]);
        dReal capsule_length = rlink[RIGHT_ELBOW_FORE_ARM].lz + rlink[RIGHT_FORE_ARM].lz - 2 * rlink[RIGHT_FORE_ARM].r;
        dsDrawCapsule(pos, dBodyGetRotation(rlink[RIGHT_FORE_ARM].id), capsule_length, rlink[RIGHT_FORE_ARM].r);
      }
      else if (i == LEFT_ELBOW_FORE_ARM || i == LEFT_FORE_ARM) {
        const dReal *pos1, *pos2;
        dReal pos[3] = {0.0, 0.0, 0.0};
        pos1 = dBodyGetPosition(rlink[LEFT_ELBOW_FORE_ARM].id);
        pos2 = dBodyGetPosition(rlink[LEFT_FORE_ARM].id);
        pos[0] = 0.5 * (pos1[0] + pos2[0]);
        pos[1] = 0.5 * (pos1[1] + pos2[1]);
        pos[2] = 0.5 * (pos1[2] + pos2[2]);
        dReal capsule_length = rlink[LEFT_ELBOW_FORE_ARM].lz + rlink[LEFT_FORE_ARM].lz - 2 * rlink[LEFT_FORE_ARM].r;
        dsDrawCapsule(pos, dBodyGetRotation(rlink[LEFT_FORE_ARM].id), capsule_length, rlink[RIGHT_FORE_ARM].r);
      }
      else {
        drawLink(rlink[i]);
      }
    }
  }

  dsSetColor(0, 0, 0);
  for (int i = 0; i < rjoint.size(); ++i) {
    if (rjoint[i].IsExist) {
      if (i == RIGHT_HIP_ROLL || i == LEFT_HIP_ROLL ||
          i == RIGHT_KNEE_PITCH || i == LEFT_KNEE_PITCH ||
          i == RIGHT_FOOT_PITCH || i == LEFT_FOOT_PITCH ||
          i == RIGHT_SHOULDER_ROLL || i == LEFT_SHOULDER_ROLL)
      {
        dVector3 pos;
        dJointGetHingeAnchor(rjoint[i].id, pos);
        dsDrawSphere(pos, dBodyGetRotation(rlink[0].id), JOINT_RADIUS);
      }
      else if (i == RIGHT_ELBOW_PITCH || i == LEFT_ELBOW_PITCH) {
        dQuaternion q;
        dMatrix3 R;
        dQFromAxisAndAngle(q, 1, 0, 0,  0.5 * M_PI);
        dQtoR(q, R);
        dVector3 pos;
        dJointGetHingeAnchor(rjoint[i].id, pos);
        dsDrawCylinder(pos, R, 0.06, 0.03);
      }
    }
  }

  // draw real COM line
  dsSetColor(2, 0.4, 1);
  int interval = 10;
  for (int i = 0; i < COMqueue.size() - interval ; i = i + interval) {
    dReal startp[3] = {COMqueue[i](0), COMqueue[i](1), COMqueue[i](2)};
    dReal endp[3] = {COMqueue[i + interval](0), COMqueue[i + interval](1), COMqueue[i + interval](2)};
    dsDrawLine (startp, endp);
  }

  // drawCOPLine();

}



