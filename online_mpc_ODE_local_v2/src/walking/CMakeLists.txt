SET( CMAKE_BUILD_TYPE Release )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
MESSAGE(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")


ENABLE_TESTING()
FIND_PACKAGE (Threads)
FIND_PACKAGE(GTest REQUIRED)

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
MESSAGE(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")


# set(Eigen_INCLUDE_DIRS /usr/include/eigen3)
set(AnkleCtrl_DIRS ${SRC_PATH}/AnkleCtrl)
set(Filters_DIRS ${SRC_PATH}/Filters)
set(IK_DIRS ${SRC_PATH}/IK)
set(LIPM_DIRS ${SRC_PATH}/LIPM)
set(LQR_DIRS ${SRC_PATH}/LQR)
# set(MatrixVector_DIRS ${SRC_PATH}/MatrixVector)
set(PreviewCtrl_DIRS ${SRC_PATH}/PreviewCtrl)
set(RobotPara_DIRS ${SRC_PATH}/RobotPara)
set(ReactiveStep_DIRS ${SRC_PATH}/ReactiveStep)
set(RTControl_DIRS ${SRC_PATH}/RTControl)
set(Stabilizer_DIRS ${SRC_PATH}/Stabilizer)
set(WBS_DIRS ${SRC_PATH}/WBS)
set(ZMPWalk_DIRS ${SRC_PATH}/ZMPWalk)
set(ComanODE_DIRS ${SRC_PATH}/ComanODE)
set(RobotModel_DIRS ${SRC_PATH}/RobotModel)
set(Utils_DIRS ${SRC_PATH}/utils)
set(QP_DIRS ${SRC_PATH}/QP)
set(ML_DIRS ${SRC_PATH}/MachineLearning)
set(MPC_DIRS ${SRC_PATH}/MPC)


set(Ext_DIRS ${SRC_PATH}/utils/External)
INCLUDE_DIRECTORIES(${Ext_DIRS}/include)
LINK_DIRECTORIES(${Ext_DIRS}/lib)

INCLUDE_DIRECTORIES(
	${Utils_DIRS}
	${SRC_PATH}
	/usr/local/include
)

FILE(GLOB SOURCES
	${SRC_PATH}/main.cpp
	${Filters_DIRS}/*.cpp
	${IK_DIRS}/*.cpp
	${Stabilizer_DIRS}/*.cpp
	${RTControl_DIRS}/RTControlBaseClass.cpp
	${RTControl_DIRS}/MpcRTControlClass.cpp
	${RobotModel_DIRS}/*.cpp
	${RobotPara_DIRS}/*.cpp
	${QP_DIRS}/*.cpp
	${WBS_DIRS}/*.cpp
	${ComanODE_DIRS}/*.cpp
	${MPC_DIRS}/*.cpp
	${Utils_DIRS}/EiQuadProg/*.cpp
	${Utils_DIRS}/XBotInterface/*.cpp
)

# yaml file path definition, will pass the file path ROBOT_PARA_YAML_PATH to the source code
ADD_DEFINITIONS( -DROBOT_PARA_YAML_PATH="${SRC_PATH}/RobotPara/" )

# ADD_DEFINITIONS( -D_GLIBCXX_USE_CXX11_ABI=0)
# ADD_DEFINITIONS( -DUSE_MRDPLOT_LOGGER)
ADD_DEFINITIONS( -DUSE_XBOT_LOGGER)
# ADD_DEFINITIONS( -DUSE_KDL)
# ADD_DEFINITIONS( -DUSE_OPENSOT)

# make own shared library and put in ZHOU_LIB_DIRS
SET(ZHOU_LIB_DIRS ${SRC_PATH}/utils/lib)
ADD_SUBDIRECTORY(${SRC_PATH}/utils)
LINK_DIRECTORIES(${ZHOU_LIB_DIRS})
# -----------------------------------------------------------

MESSAGE(STATUS "Chengxu's Library Directory: ${ZHOU_LIB_DIRS}")
MESSAGE(STATUS "External Library Directory: ${Ext_DIRS}/lib")
MESSAGE(STATUS "Project Source Directory: ${SRC_PATH}")

include_directories(include)

add_executable(${PROJ_NAME} ${SOURCES})
target_link_libraries(${PROJ_NAME} gtest gtest_main)
target_link_libraries(${PROJ_NAME} m ode drawstuff GLU GL X11 -lpthread)
target_link_libraries(${PROJ_NAME} ${catkin_LIBRARIES})
target_link_libraries(${PROJ_NAME} ${XCM_LIBRARIES} ${OpenSoT_LIBRARIES})
target_link_libraries(${PROJ_NAME} yaml-cpp urdfdom_model urdfdom_world)
target_link_libraries(${PROJ_NAME} rbdl rbdl_urdfreader)
target_link_libraries(${PROJ_NAME} matio) # for independent XBotLogger only 
target_link_libraries(${PROJ_NAME} ${CMAKE_THREAD_LIBS_INIT})
add_test(NAME unittest COMMAND unittest)
