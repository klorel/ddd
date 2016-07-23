CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(qp)


MESSAGE("CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}")
MESSAGE("CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")

SET(EXECUTABLE_OUTPUT_PATH	${CMAKE_BINARY_DIR})
SET(   LIBRARY_OUTPUT_PATH	${CMAKE_BINARY_DIR})

IF( WIN32 )

	SET( CMAKE_CONFIGURATION_TYPES "Debug;Release;" CACHE STRING "limited configs" FORCE)
	SET( CMAKE_CXX_FLAGS_DEBUG          "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1" )
	SET( CMAKE_CXX_FLAGS_RELEASE        "/MT /O2 /Ob2 /D NDEBUG")
ELSE()
	SET ( CMAKE_BUILD_TYPE "RELEASE"     )
#	SET ( CMAKE_BUILD_TYPE "DEBUG"     )
#	IF( WITH_DEBUG )
#		SET ( CMAKE_BUILD_TYPE "DEBUG"     )
#	ELSE( WITH_DEBUG )GENDIOR
#		SET ( CMAKE_BUILD_TYPE "RELEASE"     )
#	ENDIF( WITH_DEBUG )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DNDEBUG")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -pedantic -fmessage-length=0 -fPIC" )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-conversion -Wno-sign-compare")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-strict-aliasing -Wno-unused-parameter")
ENDIF( WIN32 )

SET(EIGEN_DIR $ENV{EIGENDIR})
MESSAGE("EIGEN_DIR is ${EIGEN_DIR}")

SET(XPRESS_DIR $ENV{XPRESSDIR})
MESSAGE("XPRESS_DIR is ${XPRESS_DIR}")

SET(MOSEK_DIR $ENV{MOSEKDIR})
MESSAGE("MOSEK_DIR is ${MOSEK_DIR}")

INCLUDE_DIRECTORIES(	
	cpp_pop_modeler
	cpp_pooling_modeler 
	cpp_sdp_modeler
	# cpp_opf_modeler 
	exe 
	${EIGEN_DIR}
	${XPRESS_DIR}/include
	${MOSEK_DIR}
	)

ADD_SUBDIRECTORY(cpp_pop_modeler)
ADD_SUBDIRECTORY(cpp_pooling_modeler)

ADD_SUBDIRECTORY(cpp_sdp_modeler)

#ADD_EXECUTABLE(pop exe/main.cc)
#TARGET_LINK_LIBRARIES(pop pop_modeler)

ADD_EXECUTABLE(pooling exe/pooling.cc)
TARGET_LINK_LIBRARIES(pooling pop_modeler pooling_modeler)

ADD_EXECUTABLE(sdp exe/sdp.cc)
TARGET_LINK_LIBRARIES(sdp pop_modeler sdp_modeler pooling_modeler )
#TARGET_LINK_LIBRARIES(sdp pop_modeler sdp_modeler pooling_modeler ${MOSEK_DIR}/../bin/mosek64_7_1.lib ${XPRESS_DIR}/lib/xprs.lib)
 


