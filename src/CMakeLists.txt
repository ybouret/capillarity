CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(capillarity)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_EXECUTABLE(extract_lens extract_lens.cpp)
TARGET_LINK_YOCTO(extract_lens mk graphics png z jpeg tiff lzma jbig)

INSTALL( TARGETS extract_lens RUNTIME DESTINATION bin )

ADD_SUBDIRECTORY(cbridge10)


