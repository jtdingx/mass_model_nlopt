/*****************************************************************************
BigmanODERobotClass.cpp

Description:  cpp file of BigmanODERobotClass

@Version: 1.0
@Author:  Chengxu Zhou (zhouchengxu@gmail.com)
@Release: Wed 31 Aug 2016 06:40:49 PM CEST
@Update:  Wed 31 Aug 2016 06:40:52 PM CEST
*****************************************************************************/

#include "ComanODE/BigmanODERobotClass.h"


#define JOINT_RADIUS_B   0.08

BigmanODERobotClass::BigmanODERobotClass(const double &z, const double &y, const double &x, const double &yaw)
  : ODERobotBaseClass(z, y, x, yaw)
{
  updateJointParam(RTControl.RobotPara());
  updateLinkParam(RTControl.RobotPara());
  updateRobotYaw();

  std::cout << "No." << RobotParaClass::ROBOT_COUNT() << " Bigman is ready..." << std::endl;
};

BigmanODERobotClass::~BigmanODERobotClass()
{

};

/**
 * @brief      used to update link positions and dimensions
 *
 * @param[in]  robotpara  { parameter_description }
 */
void BigmanODERobotClass::updateLinkParam(const RobotParaClass& robotpara)
{
  for (int i = 0; i < rlink.size(); ++i) {
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


  rlink[TORSO].m = robotpara.TORSO_MASS;
  rlink[HEAD].m = robotpara.HEAD_MASS;
  rlink[HEAD].type = SPHERE;
  rlink[HEAD].ds_type = SPHERE;
  rlink[HEAD].r = 0.12;
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

  rlink[TORSO].ds_type = BOX;
  rlink[TORSO].r = 0.1;
  rlink[TORSO].lx = 0.15;
  rlink[TORSO].ly = robotpara.HALF_SHOULDER_WIDTH;
  rlink[TORSO].lz = robotpara.TORSO_HEIGHT;
  rlink[TORSO].color_r = 0.0;
  rlink[TORSO].color_g = 0.0;
  rlink[TORSO].color_b = 1.0;

  rlink[RIGHT_UPPER_ARM_DUMMY1].type = CYLINDER;
  rlink[RIGHT_UPPER_ARM_DUMMY1].direction = 2; // y
  rlink[RIGHT_UPPER_ARM_DUMMY2].type = CYLINDER;
  rlink[RIGHT_UPPER_ARM_DUMMY2].r = 0.04;
  rlink[RIGHT_UPPER_ARM].ds_type = CAPSULE;
  rlink[RIGHT_ELBOW_FORE_ARM].type = CYLINDER;
  rlink[RIGHT_ELBOW_FORE_ARM].r = 0.04;
  rlink[RIGHT_FORE_ARM].r = 0.04;

  rlink[LEFT_UPPER_ARM_DUMMY1].type = CYLINDER;
  rlink[LEFT_UPPER_ARM_DUMMY1].direction = 2; // y
  rlink[LEFT_UPPER_ARM_DUMMY2].type = CYLINDER;
  rlink[LEFT_UPPER_ARM_DUMMY2].r = 0.04;
  rlink[LEFT_UPPER_ARM].ds_type = CAPSULE;
  rlink[LEFT_ELBOW_FORE_ARM].type = CYLINDER;
  rlink[LEFT_ELBOW_FORE_ARM].r = 0.04;
  rlink[LEFT_FORE_ARM].r = 0.04;

  rlink[RIGHT_THIGH].color_r = 0.65;
  rlink[RIGHT_THIGH].color_g = 0.65;
  rlink[RIGHT_THIGH].color_b = 0.65;
  rlink[RIGHT_THIGH].r = 0.068;

  rlink[LEFT_THIGH].color_r = 0.65;
  rlink[LEFT_THIGH].color_g = 0.65;
  rlink[LEFT_THIGH].color_b = 0.65;
  rlink[LEFT_THIGH].r = 0.068;

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


  rlink[RIGHT_HAND].type = SPHERE;
  rlink[RIGHT_HAND].ds_type = SPHERE;
  rlink[RIGHT_HAND].r = 0.05;
  rlink[RIGHT_HAND].color_r = 0.0;
  rlink[RIGHT_HAND].color_g = 0.0;
  rlink[RIGHT_HAND].color_b = 1.0;
  rlink[LEFT_HAND].type = SPHERE;
  rlink[LEFT_HAND].ds_type = SPHERE;
  rlink[LEFT_HAND].r = 0.05;
  rlink[LEFT_HAND].color_r = 0.0;
  rlink[LEFT_HAND].color_g = 0.0;
  rlink[LEFT_HAND].color_b = 1.0;


}

void BigmanODERobotClass::updateJointParam(const RobotParaClass& robotpara)
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


void BigmanODERobotClass::drawRobot()
{
  dReal sides[3], sides_torso1[3];
  const dReal *sides_torso2;

  for (int i = 0; i < rlink.size(); ++i) {
    if (i == RIGHT_ELBOW_FORE_ARM || i == RIGHT_FORE_ARM) {
      const dReal *pos1, *pos2;
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos1 = dBodyGetPosition(rlink[RIGHT_ELBOW_FORE_ARM].id);
      pos2 = dBodyGetPosition(rlink[RIGHT_FORE_ARM].id);
      pos[0] = pos2[0];
      pos[1] = 0.5 * (pos1[1] + pos2[1]);
      pos[2] = 0.5 * (pos1[2] + pos2[2]);
      double cylinder_length = 0.22;
      dsDrawCylinder(pos, dBodyGetRotation(rlink[i].id), cylinder_length, 0.04);
    }
    else if (i == LEFT_ELBOW_FORE_ARM || i == LEFT_FORE_ARM) {
      const dReal *pos1, *pos2;
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos1 = dBodyGetPosition(rlink[LEFT_ELBOW_FORE_ARM].id);
      pos2 = dBodyGetPosition(rlink[LEFT_FORE_ARM].id);
      pos[0] = pos2[0];
      pos[1] = 0.5 * (pos1[1] + pos2[1]);
      pos[2] = 0.5 * (pos1[2] + pos2[2]);
      double cylinder_length = 0.22;
      dsDrawCylinder(pos, dBodyGetRotation(rlink[i].id), cylinder_length, 0.04);
    }
    else if (i == RIGHT_UPPER_ARM || i == RIGHT_UPPER_ARM_DUMMY2) {
      const dReal *pos1, *pos2;
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos1 = dBodyGetPosition(rlink[RIGHT_UPPER_ARM].id);
      pos2 = dBodyGetPosition(rlink[RIGHT_UPPER_ARM_DUMMY2].id);
      pos[0] = 0.5 * (pos1[0] + pos2[0]);
      pos[1] = 0.5 * (pos1[1] + pos2[1]);
      pos[2] = 0.5 * (pos1[2] + pos2[2]);
      double cylinder_length = 0.25;
      dsDrawCylinder(pos, dBodyGetRotation(rlink[i].id), cylinder_length, 0.06);
    }
    else if (i == LEFT_UPPER_ARM || i == LEFT_UPPER_ARM_DUMMY2) {
      const dReal *pos1, *pos2;
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos1 = dBodyGetPosition(rlink[LEFT_UPPER_ARM].id);
      pos2 = dBodyGetPosition(rlink[LEFT_UPPER_ARM_DUMMY2].id);
      pos[0] = 0.5 * (pos1[0] + pos2[0]);
      pos[1] = 0.5 * (pos1[1] + pos2[1]);
      pos[2] = 0.5 * (pos1[2] + pos2[2]);
      double cylinder_length = 0.25;
      dsDrawCylinder(pos, dBodyGetRotation(rlink[i].id), cylinder_length, 0.06);
    }
    else if (i == TORSO) {
      dsSetColor(rlink[i].color_r, rlink[i].color_g, rlink[i].color_b);
      dReal sides[3];
      sides[0] = rlink[i].lx;
      sides[1] = rlink[i].ly + 0.10;
      sides[2] = rlink[i].lz + 0.03;
      const dReal *pos1;
      pos1 = dBodyGetPosition(rlink[i].id);
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos[0] = pos1[0];
      pos[1] = pos1[1];
      pos[2] = pos1[2] - 0.06;
      dsDrawBox(pos, dBodyGetRotation(rlink[i].id), sides);
    }
    else if (i == PELVIS) {
      dsSetColor(1, 1, 1);
      dReal sides[3];
      sides[0] = 0.15;
      sides[1] = 0.19;
      sides[2] = 0.25;
      const dReal *pos1;
      pos1 = dBodyGetPosition(rlink[i].id);
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos[0] = pos1[0];
      pos[1] = pos1[1];
      pos[2] = pos1[2] - 0.06;
      dsDrawBox(pos, dBodyGetRotation(rlink[i].id), sides);
    }
    else if (i == RIGHT_THIGH_DUMMY1 || i == LEFT_THIGH_DUMMY1) {
      dsSetColor(rlink[i].color_r, rlink[i].color_g, rlink[i].color_b);
      double cylinder_length = rlink[i].ly + 0.08;
      dsDrawCylinder(dBodyGetPosition(rlink[i].id), dBodyGetRotation(rlink[i].id), cylinder_length, rlink[i].r);
    }
    else if (i == HEAD) {
      dsSetColor(0.3, 0.3, 0.3);
      const dReal *pos1;
      pos1 = dBodyGetPosition(rlink[i].id);
      dReal pos[3] = {0.0, 0.0, 0.0};
      pos[0] = pos1[0];
      pos[1] = pos1[1];
      pos[2] = pos1[2] + 0.1;
      dsDrawSphere(pos, dBodyGetRotation(rlink[i].id), rlink[i].r);
    }
    else {
      drawLink(rlink[i]);
    }
  }

  dsSetColor(0.1, 0.1, 0.1);
  for (int i = 0; i < rjoint.size(); ++i) {
    // if (i == RIGHT_FOOT_PITCH || i == LEFT_FOOT_PITCH )
    // {
    //   dQuaternion q;
    //   dMatrix3 R;
    //   dQFromAxisAndAngle(q, 0, 0, 0,  0.5 * M_PI);
    //   dQtoR(q, R);
    //   dVector3 pos;
    //   dJointGetHingeAnchor(rjoint[i].id, pos);
    //   dsDrawCylinder(pos, R, 0.12, JOINT_RADIUS_B);
    // }
    // else
    if ( i == RIGHT_HIP_PITCH || i == LEFT_HIP_PITCH ||
         i == RIGHT_ELBOW_PITCH || i == LEFT_ELBOW_PITCH ||
         i == RIGHT_KNEE_PITCH || i == LEFT_KNEE_PITCH || i == WAIST_PITCH) {
      dQuaternion q;
      dMatrix3 R;
      dQFromAxisAndAngle(q, 1, 0, 0,  0.5 * M_PI);
      dQtoR(q, R);
      dVector3 pos;
      dJointGetHingeAnchor(rjoint[i].id, pos);
      if (i == RIGHT_ELBOW_PITCH || i == LEFT_ELBOW_PITCH) {
        pos[0] -= 0.05;
      }
      double len = 0.17;
      if (i == WAIST_PITCH) len = 0.23;
      dsDrawCylinder(pos, R, len, 0.07);
    }
    else if (i == RIGHT_FOOT_ROLL || i == LEFT_FOOT_ROLL ||
             i == RIGHT_SHOULDER_ROLL || i == LEFT_SHOULDER_ROLL) {
      dQuaternion q;
      dMatrix3 R;
      dQFromAxisAndAngle(q, 0, 1, 0,  0.5 * M_PI);
      dQtoR(q, R);
      dVector3 pos;
      dJointGetHingeAnchor(rjoint[i].id, pos);
      dsDrawCylinder(pos, R, 0.12, JOINT_RADIUS_B);
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



