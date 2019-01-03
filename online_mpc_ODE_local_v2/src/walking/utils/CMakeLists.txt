# # ================================================================================
# ADD_LIBRARY(ReactiveStep SHARED 
# 	${ReactiveStep_DIRS}/ReactiveStepClass.cpp
# 	${Filters_DIRS}/FilterClass.cpp
# 	${LIPM_DIRS}/LIPMDynClass.cpp
# 	${WBS_DIRS}/RobotStateClass.cpp
# 	)
# SET_TARGET_PROPERTIES(ReactiveStep PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# ADD_LIBRARY(stabilizer SHARED 
# 	${Stabilizer_DIRS}/StabilizerClass.cpp
# 	${Stabilizer_DIRS}/HandStabilizerClass.cpp
# 	${Filters_DIRS}/FilterClass.cpp
# 	${Filters_DIRS}/MeanFilterClass.cpp
# 	${WBS_DIRS}/RobotStateClass.cpp
# 	)
# SET_TARGET_PROPERTIES(stabilizer PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# ADD_LIBRARY(wbs SHARED 
# 	${RobotModel_DIRS}/RobotModelClass.cpp
# 	${WBS_DIRS}/RobotStateClass.cpp
# 	${WBS_DIRS}/WholeBodySensingClass.cpp
# 	${Filters_DIRS}/FilterClass.cpp
# 	${Filters_DIRS}/MeanFilterClass.cpp
# 	${Filters_DIRS}/KalmanFilterClass.cpp
# 	)
# SET_TARGET_PROPERTIES(wbs PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# ADD_LIBRARY(ZMPWalk SHARED 
# 	${ZMPWalk_DIRS}/ZMPGeneratorClass.cpp
# 	${ZMPWalk_DIRS}/FixedFootHoldsClass.cpp
# 	${ZMPWalk_DIRS}/swingfootcartclass.cpp
# 	)
# SET_TARGET_PROPERTIES(ZMPWalk PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# ADD_LIBRARY(PreviewCtrl SHARED 
# 	${PreviewCtrl_DIRS}/PreviewClass.cpp
# 	${LQR_DIRS}/LQRClass.cpp
# 	)
# SET_TARGET_PROPERTIES(PreviewCtrl PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")


# # ================================================================================
# ADD_LIBRARY(svm SHARED ${Utils_DIRS}/libsvm/svm.cpp)
# SET_TARGET_PROPERTIES(svm PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# ADD_LIBRARY(svmpp SHARED ${Utils_DIRS}/libsvm/svm.cpp ${Utils_DIRS}/libsvm/svmpp.cpp)
# SET_TARGET_PROPERTIES(svmpp PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# # ADD_LIBRARY(mrdplot SHARED ${Utils_DIRS}/mrdplot/mrdplot.cpp)
# # SET_TARGET_PROPERTIES(mrdplot PROPERTIES 
# #    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# # ADD_LIBRARY(EiQuadProg SHARED ${Utils_DIRS}/EiQuadProg/EiQuadProg.cpp)
# # SET_TARGET_PROPERTIES(EiQuadProg PROPERTIES 
# #    LIBRARY_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")

# # ================================================================================
# set(PSOPT_DIRS ${SRC_PATH}/utils/PSOPT)
# set(PSOPT_LIB_DIRS ${SRC_PATH}/utils/PSOPT/lib)
# INCLUDE_DIRECTORIES(${PSOPT_DIRS} ${SRC_PATH}/utils/IPOPT)
# LINK_DIRECTORIES(${PSOPT_LIB_DIRS})

# # MESSAGE(STATUS "Building lusol library")
# # ADD_DEFINITIONS( -DYZHANG -DNDEBUG )
# # INCLUDE_DIRECTORIES(${PSOPT_DIRS}/lusol)
# # FILE(GLOB LUSOL_SRC ${PSOPT_DIRS}/lusol/*.c)
# # ADD_LIBRARY(lusol STATIC ${LUSOL_SRC})
# # TARGET_LINK_LIBRARIES(lusol -Wall -fexceptions -pthread -ansi -fPIC -O3 -I -c)
# # SET_TARGET_PROPERTIES(lusol PROPERTIES 
# # 	ARCHIVE_OUTPUT_DIRECTORY "${PSOPT_LIB_DIRS}"
# # 	# LINKER_LANGUAGE C
# # )

# # MESSAGE(STATUS "Building lusol library")
# # ADD_DEFINITIONS( -DYZHANG -DNDEBUG )
# # INCLUDE_DIRECTORIES(${PSOPT_DIRS}/lusol)
# # FILE(GLOB LUSOL_SRC ${PSOPT_DIRS}/lusol/*.c)
# # ADD_LIBRARY(lusol SHARED ${LUSOL_SRC})
# # target_link_libraries(lusol ${CMAKE_DL_LIBS}) # solve undefined reference to `dlopen' and dlclose dlsym ....
# # TARGET_LINK_LIBRARIES(lusol m X11 -Wall -fexceptions -pthread -ansi -fPIC -O3 -I X11)
# # SET_TARGET_PROPERTIES(lusol PROPERTIES 
# # 	LIBRARY_OUTPUT_DIRECTORY "${PSOPT_LIB_DIRS}"
# # 	# LINKER_LANGUAGE C
# # )


# # MESSAGE(STATUS "Building dmatrix library")
# # ADD_DEFINITIONS( -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -DNDEBUG -DHAVE_MALLOC )
# # ADD_LIBRARY(dmatrix STATIC ${Utils_DIRS}/PSOPT/dmatrixv.cxx)
# # TARGET_LINK_LIBRARIES(dmatrix lusol)
# # TARGET_LINK_LIBRARIES(dmatrix -O0 -g -Wno-long-long -fomit-frame-pointer -pipe -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC )
# # SET_TARGET_PROPERTIES(dmatrix PROPERTIES 
# #    ARCHIVE_OUTPUT_DIRECTORY "${PSOPT_LIB_DIRS}")


# # MESSAGE(STATUS "Building dmatrix library")
# # ADD_DEFINITIONS( -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -DNDEBUG -DHAVE_MALLOC )
# # ADD_LIBRARY(dmatrix SHARED ${Utils_DIRS}/PSOPT/dmatrixv.cxx)
# # TARGET_LINK_LIBRARIES(dmatrix lusol cxsparse)
# # TARGET_LINK_LIBRARIES(dmatrix -O0 -g -Wno-long-long -fomit-frame-pointer -pipe -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC )
# # SET_TARGET_PROPERTIES(dmatrix PROPERTIES 
# #    LIBRARY_OUTPUT_DIRECTORY "${PSOPT_LIB_DIRS}")


# # MESSAGE(STATUS "Building ModernPsoptInterface library")
# # ADD_DEFINITIONS( -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -DNDEBUG -DHAVE_MALLOC )
# # ADD_LIBRARY(ModernPsoptInterface STATIC ${PSOPT_DIRS}/ModernPsoptInterface.cxx)
# # TARGET_LINK_LIBRARIES(ModernPsoptInterface X11 dmatrix psopt)
# # SET_TARGET_PROPERTIES(ModernPsoptInterface PROPERTIES 
# #    ARCHIVE_OUTPUT_DIRECTORY "${ZHOU_LIB_DIRS}")


# MESSAGE(STATUS "Building ModernPsoptInterface library")
# ADD_DEFINITIONS( -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -DNDEBUG -DHAVE_MALLOC )
# ADD_LIBRARY(ModernPsoptInterface SHARED ${PSOPT_DIRS}/ModernPsoptInterface.cxx)
# TARGET_LINK_LIBRARIES(ModernPsoptInterface X11 dmatrix psopt adolc ColPack )
# SET_TARGET_PROPERTIES(ModernPsoptInterface PROPERTIES 
#    LIBRARY_OUTPUT_DIRECTORY "${PSOPT_LIB_DIRS}")



