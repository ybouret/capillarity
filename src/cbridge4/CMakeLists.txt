ADD_LIBRARY(bridge4 STATIC
bridge.hpp
bridge.cpp
bridge-ode.cpp)

ADD_EXECUTABLE(cbridge4 main.cpp)
TARGET_LINK_LIBRARIES(cbridge4 bridge4)
TARGET_LINK_YOCTO(cbridge4 mk lua++)

ADD_EXECUTABLE(find-theta4 find-theta.cpp)
TARGET_LINK_LIBRARIES(find-theta4 bridge4)
TARGET_LINK_YOCTO(find-theta4 mk lua++)

INSTALL( TARGETS 
cbridge4
find-theta4 
RUNTIME DESTINATION bin )
