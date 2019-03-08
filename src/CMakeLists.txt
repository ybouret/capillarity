CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(capillarity)

IF(WIN32)
	SET(Y_ROOT /e/program/yocto/sdk)
	SET(Y_ROOT "e:/program/yocto/sdk")
ELSE()
	SET(Y_ROOT $ENV{HOME}/program/yocto/sdk)
ENDIF()

INCLUDE(${Y_ROOT}/share/yocto/yocto.cmake)

INCLUDE_DIRECTORIES(${Y_ROOT}/include)
LINK_DIRECTORIES(${Y_ROOT}/lib)

MESSAGE( STATUS "Y_ROOT=${Y_ROOT}" )
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/..)

#ADD_EXECUTABLE(contact_angle contact_angle.cpp)
#TARGET_LINK_YOCTO(contact_angle gfx tiff png lang lua++)

#ADD_EXECUTABLE(build_time build_time.cpp)
#TARGET_LINK_LIBRARIES(build_time)
#TARGET_LINK_YOCTO(build_time mk)

IF(FALSE)
ADD_EXECUTABLE(gather_angles gather_angles.cpp)
TARGET_LINK_Y(gather_angles)

ADD_EXECUTABLE(parabole parabole.cpp)
TARGET_LINK_Y(parabole mk)

ADD_EXECUTABLE(water water.cpp)
TARGET_LINK_Y(water mk)

INSTALL(TARGETS #contact_angle
gather_angles parabole water RUNTIME DESTINATION bin )
ENDIF()

ADD_SUBDIRECTORY(cbridge4)



