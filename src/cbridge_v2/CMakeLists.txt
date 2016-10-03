ADD_LIBRARY( bridge_v2 STATIC 
bridge.hpp 
bridge.cpp 
bridge-ode.cpp 
bridge-alpha.cpp
bridge-theta.cpp
par_bridge.hpp
par_bridge.cpp
)

ADD_EXECUTABLE(cbridge_v2 main.cpp)
TARGET_LINK_LIBRARIES(cbridge_v2 bridge_v2)
TARGET_LINK_YOCTO(cbridge_v2 mk lua++)

ADD_EXECUTABLE(abacus_v2 abacus_v2.cpp)
TARGET_LINK_LIBRARIES(abacus_v2 bridge_v2)
TARGET_LINK_YOCTO(abacus_v2 mk lua++)

ADD_EXECUTABLE(extract_v2 extract_v2.cpp)
TARGET_LINK_LIBRARIES(extract_v2 bridge_v2)
TARGET_LINK_YOCTO(extract_v2 mk lua++)

INSTALL( TARGETS cbridge_v2 abacus_v2 extract_v2 RUNTIME DESTINATION bin )

