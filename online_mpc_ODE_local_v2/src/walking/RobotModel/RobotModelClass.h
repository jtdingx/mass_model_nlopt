/*****************************************************************************
RobotModelClass.h

Description:	Header file of RobotModelClass

@Version:	1.0
@Author:	Chengxu Zhou (zhouchengxu@gmail.com)
@Release:	2015/10/01
@Update:	Tue 12 Apr 2016 12:03:22 PM CEST
*****************************************************************************/
#pragma once

#ifdef USE_KDL
#include <boost/shared_ptr.hpp>
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/frames_io.hpp>
#endif

#include "rbdl/rbdl.h"
namespace RBDL = RigidBodyDynamics;

#include "utils/Eigen_utils.hpp"
#include "utils/QuaternionClass.h"

#include "RobotPara/RobotParaClass.h"

#ifdef USE_KDL
class RobotChainClass
{
public:
	RobotChainClass(const RobotParaClass& robotpara, const KDL::Tree& kdl_tree, const int& link_name, const double& steptime = 0.005)
		: nDOF(0)
		, startID2(0)
		, nDOF2(0)
	{
		dt = steptime;

		EndLinkName = link_name;
		std::string chain_root = robotpara.getLinkName(PELVIS);
		std::string chain_tip = robotpara.getLinkName(EndLinkName);

		switch (EndLinkName) {
		case RIGHT_THIGH:
			startID = RIGHT_HIP_PITCH;
			nDOF = 3;
			break;
		case RIGHT_CALF:
			startID = RIGHT_HIP_PITCH;
			nDOF = 4;
			break;
		case RIGHT_FOOT:
			startID = RIGHT_HIP_PITCH;
			nDOF = 6;
			break;
		case LEFT_THIGH:
			startID = LEFT_HIP_PITCH;
			nDOF = 3;
			break;
		case LEFT_CALF:
			startID = LEFT_HIP_PITCH;
			nDOF = 4;
			break;
		case LEFT_FOOT:
			startID = LEFT_HIP_PITCH;
			nDOF = 6;
			break;
		case LEFT_HAND:
			// chain_root = robotpara.getLinkName(TORSO);
			startID = WAIST_ROLL;
			nDOF = 10;
			startID2 = LEFT_SHOULDER_PITCH;
			nDOF2 = 7;
			break;
		case RIGHT_HAND:
			// chain_root = robotpara.getLinkName(TORSO);
			startID = WAIST_ROLL;
			nDOF = 10;
			startID2 = RIGHT_SHOULDER_PITCH;
			nDOF2 = 7;
			break;
		default: std::cout << "Counld not find the link in RobotChainClass: " << chain_tip << std::endl;  break;
		}
		kdl_tree.getChain(chain_root, chain_tip, _kdl_chain);
		_kdl_jdot_solver.reset(new KDL::ChainJntToJacDotSolver(_kdl_chain));
		_kdl_jac_solver.reset(new KDL::ChainJntToJacSolver(_kdl_chain));
		_kdl_fk_solver.reset(new KDL::ChainFkSolverPos_recursive(_kdl_chain));
		_q_in.resize(_kdl_chain.getNrOfJoints());
		_qdot_in.resize(_kdl_chain.getNrOfJoints());
		nSeg = _kdl_chain.getNrOfSegments();
		EndSegName = _kdl_chain.getSegment(nSeg - 1).getName();
		assert(chain_tip == EndSegName);

	};
	~RobotChainClass() {};

	void Update(const std::vector<double>& qall, const std::vector<double>& dqall = std::vector<double>(31, 0.0))
	{
		if (startID != WAIST_ROLL) {
			for (int j = 0; j < nDOF; j++) {
				// _qdot_in(j) = dqall[j + startID];
				_qdot_in(j) = (qall[j + startID] - _q_in(j)) / dt;
				_q_in(j) = qall[j + startID];
			}
		}
		else {
			for (int j = 0; j < 3; j++) {
				_qdot_in(j) = (qall[j] - _q_in(j)) / dt;
				_q_in(j) = qall[j];
			}
			for (int j = 0; j < nDOF2; j++) {
				// _qdot_in(j) = dqall[j + startID];
				_qdot_in(j + 3) = (qall[j + startID2] - _q_in(j + 3)) / dt;
				_q_in(j + 3) = qall[j + startID2];
			}
		}
	};

	Eigen::MatrixXd Jdot()
	{
		KDL::Jacobian jdot_kdl(nDOF);
		_kdl_jdot_solver->JntToJacDot(KDL::JntArrayVel(_q_in, _qdot_in), jdot_kdl);
		return jdot_kdl.data;
	};

	Eigen::Vector6d Jdotqdot()
	{
		KDL::Twist jac_dot_q_dot;
		_kdl_jdot_solver->JntToJacDot(KDL::JntArrayVel(_q_in, _qdot_in), jac_dot_q_dot);
		return vKDLtoEigen(jac_dot_q_dot);
	};

	Eigen::MatrixXd Jacobian()
	{
		// std::vector<double> jntAngle(qall.begin()+startID, qall.begin()+startID+nDOF);

		KDL::Jacobian jacobian_kdl(nDOF);
		_kdl_jac_solver->JntToJac(_q_in, jacobian_kdl);
		// COUT(EndSegName, jacobian_kdl.data.rows(),jacobian_kdl.data.cols());
		return jacobian_kdl.data;
	};

	KDL::Frame ChainFK()
	{
		_kdl_fk_solver->JntToCart(_q_in, EndEffector);
		return EndEffector;
	};

	// KDL::Frame EndSegment()
	// {
	// 	_kdl_fk_solver->JntToCart(_q_in, EndEffector);
	// 	return EndEffector;
	// };
	//
	void PrintSelf()
	{
		std::cout << "+++++++++++++++++  " << EndSegName << " \t+++++++++++++++++" << std::endl;
		std::cout << "EndEffector Position: " << EndEffector.p << std::endl;

		std::cout << "It has " << nSeg << " segments, which are: " << std::endl;
		for (int i = 0; i < nSeg; ++i) {
			std::cout << "\t" << i << "\t" << _kdl_chain.getSegment(i).getName() << "\t" << _kdl_chain.getSegment(i).getJoint().getName() << std::endl;
		}
		std::cout << std::endl;
	};

	int startID, nDOF;
	int startID2, nDOF2;
	KDL::Chain _kdl_chain;

	int nSeg;
	std::string EndSegName;
	int EndLinkName;


private:
	boost::shared_ptr<KDL::ChainJntToJacDotSolver> _kdl_jdot_solver;
	boost::shared_ptr<KDL::ChainJntToJacSolver> _kdl_jac_solver;
	boost::shared_ptr<KDL::ChainFkSolverPos_recursive> _kdl_fk_solver;

	KDL::Frame EndEffector;
	KDL::JntArray _q_in;
	KDL::JntArray _qdot_in;

	double dt;

};




class LinkLine
{
public:
	LinkLine(): P0(0, 0, 0), P1(0, 0, 0), radius(0.0) {};
	~LinkLine() {};

	void setPara(const std::string& linename, const std::string& segname, const double& r) {name = linename; EndSegName = segname; radius = r;};

	void setPara(const std::string& linename, const int& linkname, const double& r) {name = linename; EndLinkName = linkname; radius = r;};

	void Update(Eigen::Vector3d L0, Eigen::Vector3d L1) {P0 = L0; P1 = L1;};

	inline void UpdateJacobian(Eigen::MatrixXd& jaco) {Jacobian = jaco;};

	void PrintSelf()
	{
		std::cout << "+++++++++++++++++  " << name << " \t+++++++++++++++++" << std::endl;
		std::cout << "P0: " << P0.transpose() << std::endl;
		std::cout << "P1: " << P1.transpose() << std::endl;

		std::cout << std::endl;
	};

	Eigen::MatrixXd Jacobian;
	Eigen::Vector3d P0, P1;
	double radius;
	std::string name;
	std::string EndSegName;
	int EndLinkName;
};

class LinkPair
{
public:
	LinkPair(LinkLine& s1, LinkLine& s2)
		: S1(s1)
		, S2(s2)
		, CP1(0, 0, 0)
		, CP2(0, 0, 0)
		, d_i(0.05)
		, d_s(0.005)
		, damping(0.3)
		, nDOF(12)
		, _B_row(0.0)
		, minDis(0.0)
	{
		_A_row = Eigen::RowVectorXd::Zero(nDOF);
	};
	~LinkPair() {};

	void PrintSelf()
	{
		std::cout << "+++++++++++++++++  " << S1.name << "   " << S2.name << " \t+++++++++++++++++" << std::endl;
		std::cout << "Link Pair MIN Distance: " << dist3D_Capsule_to_Capsule() << std::endl;
		std::cout << "A_row: " << _A_row << std::endl;
		std::cout << "B_row: " << _B_row << std::endl;
		std::cout << "closepoint_dir: " << closepoint_dir.transpose() << std::endl;
		std::cout << "CP1: " << CP1.transpose() << std::endl;
		std::cout << "CP2: " << CP2.transpose() << std::endl;

		std::cout << std::endl;
	};

	inline void Damping(const double& value) {damping = value;};
	inline void InfluenceDistance(const double& value) {d_i = value;};
	inline void SecurityDistance(const double& value) {d_s = value;};
	inline double Damping() {return damping;};
	inline double InfluenceDistance() {return d_i;};
	inline double SecurityDistance() {return d_s;};
	inline bool IsTooClose() {return (dist3D_Capsule_to_Capsule() < d_i);};

	double dist3D_Capsule_to_Capsule()
	{
		Eigen::Vector3d cp1, cp2, cp12;
		dist3D_Segment_to_Segment(S1, S2, cp1, cp2);
		cp12 = cp2 - cp1;
		CP1 = cp1 + S1.radius * (cp12 / cp12.norm());
		CP2 = cp2 + S2.radius * (-cp12 / cp12.norm());
		Eigen::Vector3d  dP = CP2 - CP1;
		closepoint_dir = dP / dP.norm();
		minDis = dP.norm();
		return minDis;
	};

	void Update()
	{
		Eigen::MatrixXd Link1_CP_Jaco(6, nDOF), Link2_CP_Jaco(6, nDOF);
		Link1_CP_Jaco.setZero();
		Link2_CP_Jaco.setZero();
		Link1_CP_Jaco.block(0, 0, 6, S1.Jacobian.cols()) = S1.Jacobian;
		Link2_CP_Jaco.block(0, 6, 6, S2.Jacobian.cols()) = S2.Jacobian;

		Link1_CP_Jaco = skewSymmetricOperator(CP1 - S1.P0) * Link1_CP_Jaco;
		Link2_CP_Jaco = skewSymmetricOperator(CP2 - S2.P0) * Link2_CP_Jaco;

		_A_row = closepoint_dir.transpose() * (-Link1_CP_Jaco + Link2_CP_Jaco);
		_B_row = -damping * (minDis - d_s) / (d_i - d_s);
	};

	inline void setDistance(double& di, double& ds) {d_i = di; d_s = ds;};
	inline void setDamping(double& damp) {damping = damp;};
	inline Eigen::RowVectorXd A_row() {return _A_row;};
	inline double B_row() {return _B_row;};

	LinkLine &S1, &S2;
	double d_i, d_s;
	double damping;
	Eigen::Vector3d CP1, CP2;
	Eigen::Vector3d closepoint_dir;
	Eigen::RowVectorXd _A_row;
	double _B_row;
	double minDis;
	const int nDOF;


private:
	Eigen::MatrixXd skewSymmetricOperator(const Eigen::Vector3d& r_cp)
	{
		Eigen::MatrixXd J_transform(3, 6);
		Eigen::Matrix3d left_part, right_part;
		left_part = Eigen::MatrixXd::Identity(3, 3);
		right_part << 0, r_cp(2), -r_cp(1),
		           -r_cp(2), 0, r_cp(0),
		           r_cp(1), -r_cp(0), 0;
		J_transform.block(0, 0, 3, 3) = left_part;
		J_transform.block(0, 3, 3, 3) = right_part;

		return J_transform;
	}

	double dist3D_Segment_to_Segment(const LinkLine& S1, const LinkLine& S2, Eigen::Vector3d & CP1, Eigen::Vector3d & CP2)
	{
		double SMALL_NUM = 1e-10;
		Eigen::Vector3d   u = S1.P1 - S1.P0;
		Eigen::Vector3d   v = S2.P1 - S2.P0;
		Eigen::Vector3d   w = S1.P0 - S2.P0;
		double    a = u.dot(u);         // always >= 0
		double    b = u.dot(v);
		double    c = v.dot(v);         // always >= 0
		double    d = u.dot(w);
		double    e = v.dot(w);
		double    D = a * c - b * b;    // always >= 0
		double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
		double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0


		// compute the line parameters of the two closest points
		if (D < SMALL_NUM) { // the lines are almost parallel
			sN = 0.0;         // force using point P0 on segment S1
			sD = 1.0;         // to prevent possible division by 0.0 later
			tN = e;
			tD = c;
		}
		else {                 // get the closest points on the infinite lines
			sN = (b * e - c * d);
			tN = (a * e - b * d);
			if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
				sN = 0.0;
				tN = e;
				tD = c;
			}
			else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}

		if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
			tN = 0.0;
			// recompute sc for this edge
			if (-d < 0.0)
				sN = 0.0;
			else if (-d > a)
				sN = sD;
			else {
				sN = -d;
				sD = a;
			}
		}
		else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
			tN = tD;
			// recompute sc for this edge
			if ((-d + b) < 0.0)
				sN = 0;
			else if ((-d + b) > a)
				sN = sD;
			else {
				sN = (-d +  b);
				sD = a;
			}
		}
		// finally do the division to get sc and tc
		sc = (std::abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
		tc = (std::abs(tN) < SMALL_NUM ? 0.0 : tN / tD);
		// std::cout<<"abs(sN): "<<std::abs(sN)<<"\t SMALL_NUM: "<<SMALL_NUM<<std::endl;
		// std::cout<<"D: "<<D<<"\t abs(sN): "<<std::abs(sN)<<"\t sD: "<<sD<<"\t sN: "<<sN<<std::endl;

		// std::cout<<"sc: "<<sc<<"\t"<<tc<<std::endl;

		CP1 = S1.P0 + sc * u;
		CP2 = S2.P0 + tc * v;

		// get the difference of the two closest points
		Eigen::Vector3d   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

		return dP.norm();   // return the closest distance
	};

};

//====================================================================
class KDLModelClass
{
public:
	KDLModelClass();
	virtual ~KDLModelClass() {};

	void InitKDL(const RobotParaClass& robotpara);
	void UpdateKDL(const std::vector<double>& qall_msr);
	void UpdateCollisionInEq(Eigen::MatrixXd& CI, Eigen::Vector3d& ci);
	void WhichFoot(const int& which, const double& percent);

	inline const Eigen::MatrixXd& CollisionInEqA() {return _A_LBn_eigen;};
	inline const Eigen::VectorXd& CollisionInEqB() {return _lbA_LBn_eigen;};
	inline int nCollisionInEq() {return nC_self_pair + nC_Obs_pair;};
	inline double getColWeight() {return We;};

protected:
	KDL::Tree _kdl_tree;

	Eigen::MatrixXd getJacobian(const int& body_id);
	Eigen::MatrixXd getJdot(const int& body_id);
	Eigen::Vector6d getJdotQdot(const int& body_id);

private:
	void SelfCollisionCheck();
	double CalcManipulability(const int& body_id);
	double CalcManipulability(const Eigen::MatrixXd& jaco);
	double FadeOutRatio(double whichfoot, double percent, double BeginSS, double EndSS, double Interval);
	double Poly3(double ds0, double s0, double t0, double dsf, double sf, double tf, double time);
	inline double setWeight(const double &scalar, const double &Lgain, const double &Hgain) {return Hgain - scalar * (Hgain - Lgain);};


	double scalar;
	double scalar2;

	boost::shared_ptr<RobotChainClass> RHandChain, LHandChain;

	boost::shared_ptr<RobotChainClass> RThighChain, RCalfChain, RFootChain, LThighChain, LCalfChain, LFootChain;

	KDL::Frame RThighFrame, RCalfFrame, RFootFrame, LThighFrame, LCalfFrame, LFootFrame;

	KDL::Frame FootP1, FootP2, FootP3, FootP4;
	KDL::Frame AnkleProjection;

	LinkLine RThighLowLine, RCalfLine, LThighLowLine, LCalfLine;
	LinkLine RFootLine1, RFootLine2, RFootLine3, RFootLine4;
	LinkLine LFootLine1, LFootLine2, LFootLine3, LFootLine4;
	LinkLine RAnkleLine, LAnkleLine;

	std::vector<LinkPair> SelfCollisionPair;

	Eigen::MatrixXd _A_LBn_eigen;
	Eigen::VectorXd _lbA_LBn_eigen;

	int nDOF, nC_self_pair, nC_Obs_pair;
	double IsLftSwing, IsRftSwing, We;

	double dt;
};
#endif

class RBDLModelClass
{
public:
	RBDLModelClass() {};
	virtual ~RBDLModelClass() {};

	void InitRBDL(const RobotParaClass& robotpara);
	void UpdateRBDL(const std::vector<double>& qall_msr);
	void UpdateRBDL(const Eigen::VectorXd& qall_msr);

	int getJointNameFromRBDLJntID(const int& rbdlID)
	{
		for (auto it : _map_JointName_to_rbdl) {
			if (it.second == rbdlID) return it.first;
		};
	};

	inline const int& getRBDLJointID(const int& joint_name) {return _map_JointName_to_rbdl[joint_name];};
	inline const int& getRBDLBodyID(const int& link_name) {return _map_LinkName_to_rbdl[link_name];};

	RBDL::Math::SpatialTransform getLocalBodyFrame(const int& link_name);
	Eigen::Vector3d PositionFromGlobal2Hip(const Eigen::Vector3d& from);
	Eigen::Vector3d VelocityFromGlobal2Hip(const Eigen::Vector3d& vel_from);
	Eigen::Matrix3d getGloballBodyFrame(const int& link_name);
	Eigen::Vector3d getGloballBodyPosition(const int& link_name);
	Eigen::Vector3d getGloballBodyVelocity(const int& link_name);

	const RBDL::Math::Matrix3d& getOrientationFromWorldToBody(const int& link_name);
	const RBDL::Math::MatrixNd& InertiaMatrix();
	const RBDL::Math::VectorNd& NonlinearEffects();
	// Eigen::VectorXd sortVecFromRBDLtoJointName(const Eigen::VectorXd& from);

	template <typename T>
	void getRequiredTorques(T& torque_all) // if it is the template func used by other cpp file, the implementation should be in the header file, otherwise will cause " undefined reference to" error because other cpp file cannot see the implementation
	{
		InverseDynamics();
		vRBDLToJointName(tau_all_floating, torque_all, 0, torque_all.size() - 1);
	};

	void CalcCenterOfMass(const bool& IsUpdateKinematics = true);
	void InverseDynamics();
	Eigen::VectorXd InverseDynamicsContacts(RBDL::ConstraintSet& CS);
	void ForwardDynamicsContacts(RBDL::ConstraintSet& CS);
	// Eigen::MatrixXd CMM() {return centroidal_momentum_matrix();};
	Eigen::MatrixXd CalcCMM();

	/** \brief Returns the spatial jacobian of given frame with respect to required reference frame
	 *
	 * \param tar_body_id id of the target body to which the target frame are expressed
	 * \param tar_body_Rotation_tar_frame the orientation of target frame with respect to target body frame
	 * \param tar_body_Translation_tar_frame the position of target frame expressed in body frame
	 * \param ref_body_id id of the body to which the spatial velocity are expressed (default: 0, which means the reference frame is base frame).
	 * \param ref_body_Rotation_ref_frame the orientation of reference frame with respect to reference body frame (defual: I, no rotation)
	 * \param ref_body_Translation_ref_frame the position of reference frame expressed in reference body frame (defual: 0, no translation)
	 * \param IsUpdateKinematics whether UpdateKinematics() should be called or not
	 * (default: false).
	 *
	 * \returns a jacobian matrix of the target frame wrt reference frame (rotational part above translational part)
	 */
	Eigen::MatrixXd CalcFrameJacobian ( unsigned int tar_body_id,
	                                    RBDL::Math::Matrix3d tar_body_Rotation_tar_frame,
	                                    RBDL::Math::Vector3d tar_body_Translation_tar_frame,
	                                    unsigned int ref_body_id = 0,
	                                    RBDL::Math::Matrix3d ref_body_Rotation_ref_frame = RBDL::Math::Matrix3d::Identity(),
	                                    RBDL::Math::Vector3d ref_body_Translation_ref_frame = RBDL::Math::Vector3d::Zero(),
	                                    bool IsUpdateKinematics = false );


	/** \brief Returns the spatial transform of given frame with respect to required reference frame
	 *
	 * \param model the rigid body model
	 * \param Q the curent genereralized positions
	 * \param tar_body_id id of the target body to which the target frame are expressed
	 * \param tar_body_Rotation_tar_frame the orientation of target frame with respect to target body frame
	 * \param tar_body_Translation_tar_frame the position of target frame expressed in body frame
	 * \param ref_body_id id of the body to which the spatial velocity are expressed (default: 0, which means the reference frame is base frame).
	 * \param ref_body_Rotation_ref_frame the orientation of reference frame with respect to reference body frame (defual: I, no rotation)
	 * \param ref_body_Translation_ref_frame the position of reference frame expressed in reference body frame (defual: 0, no translation)
	 *
	 * \returns a spatial transform of the target frame wrt reference frame(spactial transform: E->rotation matrix, r->translation vector)
	 */
	RBDL::Math::SpatialTransform CalcFrameTransform (
	    unsigned int tar_body_id,
	    RBDL::Math::Matrix3d tar_body_Rotation_tar_frame,
	    RBDL::Math::Vector3d tar_body_Translation_tar_frame,
	    unsigned int ref_body_id = 0,
	    RBDL::Math::Matrix3d ref_body_Rotation_ref_frame = RBDL::Math::Matrix3d::Identity(),
	    RBDL::Math::Vector3d ref_body_Translation_ref_frame = RBDL::Math::Vector3d::Zero());

	/**
	 * @brief      transfer a vector with total dof length from predefined joint name order to rbdl order (6 floating dof + joint num)
	 *
	 * @param[in]  from     The from
	 * @param      to_rbdl  To rbdl
	 *
	 * @tparam     T1       std::vector or Eigen::VectorXd
	 * @tparam     T2       std::vector or Eigen::VectorXd
	 */
	template <typename T1, typename T2>
	void vJointNameToRBDL(const T1& from, T2& to_rbdl)
	{
		for (auto it : _map_JointName_to_rbdl) {
			if (RBDL_API_VERSION == rbdl250) {
				to_rbdl[it.second + 4] = from[it.first]; // for rbdl 2.5.0
			}
			else if (RBDL_API_VERSION == rbdl240) {
				to_rbdl[it.second] = from[it.first]; // for rbdl 2.4.0
			}
		}
	};

	/**
	 * @brief      sort vector "from_rbdl" to "to" from rbdl order to link_name order
	 *
	 * @param[in]  from_rbdl  The from rbdl
	 * @param      to         resulted vector in link_name order
	 * @param[in]  start      The start link_name id
	 * @param[in]  end        The end link_name id
	 *
	 * @tparam     T1         std::vector or Eigen::VectorXd
	 * @tparam     T2         std::vector or Eigen::VectorXd
	 */
	template <typename T1, typename T2>
	void vRBDLToJointName(const T1& from_rbdl, T2& to, const int& start = 0, const int& end = RobotParaClass::JOINT_NUM() - 1)
	{
		for (auto it : _map_JointName_to_rbdl) {
			if ((it.first >= start) && (it.first <= end)) {
				assert(it.second < from_rbdl.size());
				assert(it.first < to.size());
				// to[it.first] = from_rbdl[it.second]; // for rbdl 2.4.0, has 6 dof at the start
				to[it.first] = from_rbdl[it.second]; // for rbdl 2.5.0
			}
		}
	};

	RBDL::ConstraintSet constraint_set_lft, constraint_set_rft, constraint_set_bothft;

	RBDL::Math::MatrixNd Jcom;
	RBDL::Math::MatrixNd dJcom;
	RBDL::Math::MatrixNd _CMM;
	RBDL::Math::MatrixNd _dCMM;
	RBDL::Math::VectorNd q_all_floating;
	RBDL::Math::VectorNd dq_all_floating;
	RBDL::Math::VectorNd ddq_all_floating;
	RBDL::Math::VectorNd tau_all_floating;
	RBDL::Math::Vector3d com_floating, com_ft;
	RBDL::Math::Vector3d com_hip, dcom_hip;
	RBDL::Math::Vector3d com_velocity, angular_momentum;
	RBDL::Math::VectorNd q_all_fixed;
	RBDL::Math::VectorNd dq_all_fixed;
	double Ep, Ek; // potential and kinetic energy of the robot

	double mass;

	RBDL::Model* getRBDLModel() {return &_rbdl_model;};

	int rbdl240, rbdl250;

	Eigen::Matrix3d Rot_World_to_Pelvis;
	RBDL::Math::SpatialTransform local_origin_PelvisProj;
	RBDL::Math::SpatialTransform lwrist, rwrist;
	RBDL::Math::SpatialTransform lsole, rsole;
	RBDL::Math::SpatialTransform lankle, rankle;
	Eigen::Vector3d base_pos;

	void UpdateKinematicsOnce();

protected:
	RBDL::Model _rbdl_model;

	Eigen::MatrixXd getJacobian(const int& body_id, const bool& IsFullJacobian = false, const bool& IsUpdateKinematics = false);

	Eigen::MatrixXd getJdot(const int& body_id, const bool& IsFullJacobian = false);

	Eigen::VectorXd getJdotQdot(const int& body_id, const bool& IsFullJacobian = false);

private:
	double dt;
	double nDoF_floating;
	double nLink_floating;

	Eigen::Vector3d _ankle_offset;

	RBDL::Math::MatrixNd _H; //! inertial matrix
	RBDL::Math::VectorNd _C; //! gravity + centrifugal and coriolis forces

	RBDL::Math::VectorNd q_all_floating_old, dq_all_floating_old;

	std::map<int, int> _map_LinkName_to_rbdl, _map_JointName_to_rbdl;

	RBDL::Math::MatrixNd _CMM_old;
	Eigen::MatrixXd compute_P();
	Eigen::MatrixXd compute_J();
	Eigen::MatrixXd centroidal_momentum_matrix();
};

//===============================================================

#ifdef USE_KDL
class RobotModelClass: public KDLModelClass, public RBDLModelClass
#else
class RobotModelClass: public RBDLModelClass
#endif
{
public:
	RobotModelClass() {};
	virtual ~RobotModelClass() {};

	int Init(const RobotParaClass& robotpara)
	{
		InitRBDL(robotpara);

#ifdef USE_KDL
		InitKDL(robotpara);
#endif
		// COUT("KDLModelClass::dt:",KDLModelClass::dt);
		// COUT("RBDLModelClass::dt:",RBDLModelClass::dt);
	};

	int Update(const std::vector<double>& qall_msr)
	{
		UpdateRBDL(qall_msr);
#ifdef USE_KDL
		UpdateKDL(qall_msr);
#endif
	};

	int Update(const Eigen::VectorXd& qall_msr)
	{
		UpdateRBDL(qall_msr);
#ifdef USE_KDL
		UpdateKDL(qall_msr);
#endif
	};

	Eigen::MatrixXd getJacobian(const int& body_id, const std::string& lib_name = "kdl", const bool& IsFullJacobian = false)
	{
		if (lib_name == "rbdl") {
			return RBDLModelClass::getJacobian(body_id, IsFullJacobian);
		}
#ifdef USE_KDL
		else if (lib_name == "kdl") {
			assert(!IsFullJacobian);
			return KDLModelClass::getJacobian(body_id);
		}
#endif
		else {
			assert(!"not defined lib_name to getJacobian, should be kdl or rbdl");
		}
	};

	Eigen::MatrixXd getJdot(const int& body_id, const std::string& lib_name = "kdl", const bool& IsFullJacobian = false)
	{
		if (lib_name == "rbdl") {
			return RBDLModelClass::getJdot(body_id, IsFullJacobian);
		}
#ifdef USE_KDL
		else if (lib_name == "kdl") {
			assert(!IsFullJacobian);
			return KDLModelClass::getJdot(body_id);
		}
#endif
		else {
			assert(!"not defined lib_name to getJdot, should be kdl or rbdl");
		}
	};

	Eigen::MatrixXd getJdotQdot(const int& body_id, const std::string& lib_name = "kdl", const bool& IsFullJacobian = false)
	{
		if (lib_name == "rbdl") {
			return RBDLModelClass::getJdotQdot(body_id, IsFullJacobian);
		}
#ifdef USE_KDL
		else if (lib_name == "kdl") {
			assert(!IsFullJacobian);
			return KDLModelClass::getJdotQdot(body_id);
		}
#endif
		else {
			assert(!"not defined lib_name to getJdotQdot, should be kdl or rbdl");
		}
	};

private:


};


