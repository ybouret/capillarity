ADD_LIBRARY( bridge_v2 STATIC 
bridge.hpp 
bridge.cpp 
bridge-ode.cpp 
bridge-alpha.cpp
bridge-theta.cpp
par_bridge.hpp
par_bridge.cpp
)

#ADD_EXECUTABLE(cbridge_v2 main.cpp)
#TARGET_LINK_LIBRARIES(cbridge_v2 bridge_v2)
#TARGET_LINK_YOCTO(cbridge_v2 mk lua++)

ADD_EXECUTABLE(find_theta find_theta.cpp)
TARGET_LINK_LIBRARIES(find_theta bridge_v2)
TARGET_LINK_YOCTO(find_theta mk lua++)

ADD_EXECUTABLE(find_theta_v2 find_theta_v2.cpp)
TARGET_LINK_LIBRARIES(find_theta_v2 bridge_v2)
TARGET_LINK_YOCTO(find_theta_v2 mk lua++)

ADD_EXECUTABLE(find_alpha find_alpha.cpp)
TARGET_LINK_LIBRARIES(find_alpha bridge_v2)
TARGET_LINK_YOCTO(find_alpha mk lua++)


ADD_EXECUTABLE(abacus_v2 abacus_v2.cpp)
TARGET_LINK_LIBRARIES(abacus_v2 bridge_v2)
TARGET_LINK_YOCTO(abacus_v2 mk lua++)

ADD_EXECUTABLE(extract_v2 extract_v2.cpp)
TARGET_LINK_LIBRARIES(extract_v2 bridge_v2)
TARGET_LINK_YOCTO(extract_v2 mk lua++)

IF(OFF)
INSTALL( TARGETS 
abacus_v2 
extract_v2 
find_alpha 
find_theta 
find_theta_v2
RUNTIME DESTINATION bin )
ENDIF()


