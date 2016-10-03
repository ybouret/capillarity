ADD_LIBRARY( bridge_v2 STATIC 
bridge.hpp 
bridge.cpp 
bridge-ode.cpp 
bridge-alpha.cpp
bridge-theta.cpp)

ADD_EXECUTABLE(cbridge_v2 main.cpp)
TARGET_LINK_LIBRARIES(cbridge_v2 bridge_v2)
TARGET_LINK_YOCTO(cbridge_v2 mk lua++)

ADD_EXECUTABLE(abacus_v2 abacus_v2.cpp)
TARGET_LINK_LIBRARIES(abacus_v2 bridge_v2)
TARGET_LINK_YOCTO(abacus_v2 mk lua++)

ADD_EXECUTABLE(h_of_xi h_of_xi.cpp)
TARGET_LINK_LIBRARIES(h_of_xi bridge_v2)
TARGET_LINK_YOCTO(h_of_xi mk lua++)

INSTALL( TARGETS cbridge_v2 abacus_v2 h_of_xi RUNTIME DESTINATION bin )

