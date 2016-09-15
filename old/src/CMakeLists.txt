CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(capillarity)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_EXECUTABLE(extract_lens extract_lens.cpp)
TARGET_LINK_YOCTO(extract_lens mk gfx png tiff)

ADD_EXECUTABLE(polar polar.cpp)
TARGET_LINK_YOCTO(polar mk)

#ADD_EXECUTABLE(contact_angle contact_angle.cpp)
#TARGET_LINK_YOCTO(contact_angle mk graphics png tiff)

#ADD_EXECUTABLE(contact_angle2 contact_angle2.cpp)
#TARGET_LINK_YOCTO(contact_angle2 mk gfx png tiff)

#ADD_EXECUTABLE(contact_angle3 contact_angle3.cpp)
#TARGET_LINK_YOCTO(contact_angle3 mk gfx png tiff)

ADD_EXECUTABLE(contact_angle_v4 contact_angle_v4.cpp)
TARGET_LINK_YOCTO(contact_angle_v4 mk gfx png tiff)

INSTALL( TARGETS polar extract_lens contact_angle_v4 RUNTIME DESTINATION bin )

ADD_SUBDIRECTORY(cbridge10)

