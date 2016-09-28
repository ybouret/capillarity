ADD_LIBRARY( bridge_v2 STATIC bridge.hpp bridge.cpp bridge-ode.cpp bridge-fcn.cpp)

ADD_EXECUTABLE(cbridge_v2 main.cpp)
TARGET_LINK_LIBRARIES(cbridge_v2 bridge_v2)
TARGET_LINK_YOCTO(cbridge_v2 mk lua++ seem)

INSTALL( TARGETS cbridge_v2 RUNTIME DESTINATION bin )

