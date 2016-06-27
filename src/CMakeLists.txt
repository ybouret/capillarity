CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(capillarity)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_EXECUTABLE(extract_lens extract_lens.cpp)
TARGET_LINK_YOCTO(extract_lens mk graphics png tiff)

#ADD_EXECUTABLE(contact_angle contact_angle.cpp)
#TARGET_LINK_YOCTO(contact_angle mk graphics png tiff)

ADD_EXECUTABLE(contact_angle2 contact_angle2.cpp)
TARGET_LINK_YOCTO(contact_angle2 mk graphics png tiff)


INSTALL( TARGETS extract_lens contact_angle2 RUNTIME DESTINATION bin )

ADD_SUBDIRECTORY(cbridge10)


