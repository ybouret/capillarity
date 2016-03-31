ADD_EXECUTABLE(profile profile.cpp)
TARGET_LINK_YOCTO(profile mk)

INSTALL( TARGETS profile RUNTIME DESTINATION bin )

ADD_EXECUTABLE(bridge bridge.hpp bridge.cpp main.cpp)
TARGET_LINK_YOCTO(bridge mk)


