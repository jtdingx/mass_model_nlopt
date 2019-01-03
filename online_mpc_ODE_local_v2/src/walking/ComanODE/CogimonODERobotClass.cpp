/*****************************************************************************
CogimonODERobotClass.cpp

Description:  cpp file of CogimonODERobotClass

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: Thu 30 Nov 2017 08:30:19 PM CET
@Update:  Thu 30 Nov 2017 08:30:15 PM CET
*****************************************************************************/

#include "ComanODE/CogimonODERobotClass.h"


#define JOINT_RADIUS_B   0.08

CogimonODERobotClass::CogimonODERobotClass(const double &z, const double &y, const double &x, const double &yaw)
  : ODERobotBaseClass(z, y, x, yaw)
{
  std::cout << "Start to construct Cogimon robot \n";
  updateJointParam(RTControl.RobotPara());
  updateLinkParam(RTControl.RobotPara());
  updateRobotYaw();

  std::cout << "No." << RobotParaClass::ROBOT_COUNT() << " Cogimon is ready..." << std::endl;
};

CogimonODERobotClass::~CogimonODERobotClass()
{

};

/**
 * @brief      used to update link positions and dimensions
 *
 * @param[in]  robotpara  { parameter_description }
 */
void CogimonODERobotClass::updateLinkParam(const RobotParaClass& robotpara)
{
  for (int i = 0; i < rlink.size(); ++i) {
    if (robotpara.getLink(i)) {
      bool HasChild = !robotpara.getLink(i)->child_joints.empty();
      COUT("\n", i, rlink[i].name, HasChild);
      if (HasChild && robotpara.getLink(i)->child_joints[0]->type == 1)
      {
        if (i == PELVIS) {
          rlink[i].type = CYLINDER;
          rlink[i].ds_type = CYLINDER;
          rlink[i].lx = 0.1;
          rlink[i].ly = 0.1;
          rlink[i].lz = 0.1; // + offset for visualization
          rlink[i].px = start_x;
          rlink[i].py = start_y;
          rlink[i].pz = start_z;
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

        // // COUT(rjoint[parent].name);
        // // COUT(rjoint[parent].px, rjoint[parent].py, rjoint[parent].pz);

        // COUT(rlink[i].name, rlink[i].m);
        // COUT(rlink[i].px, rlink[i].py, rlink[i].pz);
        // COUT(rlink[i].lx, rlink[i].ly, rlink[i].lz);

        // // COUT(rjoint[child].name);
        // // COUT(rjoint[child].px, rjoint[child].py, rjoint[child].pz);
        // COUT(" ");
      }
    }
  }
  rlink[RIGHT_THIGH_DUMMY1].type = CYLINDER;
  rlink[RIGHT_THIGH_DUMMY1].direction = 1;

  rlink[LEFT_THIGH_DUMMY1].type = CYLINDER;
  rlink[LEFT_THIGH_DUMMY1].direction = 1;

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

  rlink[TORSO].type = CYLINDER;
  rlink[TORSO].ds_type = CYLINDER;
  rlink[TORSO].lx = 0.1;
  rlink[TORSO].ly = 0.1;
  rlink[TORSO].lz = 0.1;
  rlink[TORSO].px = rlink[PELVIS].px;
  rlink[TORSO].py = rlink[PELVIS].py;
  rlink[TORSO].pz = rlink[PELVIS].pz +  0.3;

}

void CogimonODERobotClass::updateJointParam(const RobotParaClass& robotpara)
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

}


void CogimonODERobotClass::drawRobot()
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

  drawCOPLine();
}



