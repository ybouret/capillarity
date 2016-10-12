ADD_LIBRARY(bridge4 STATIC
bridge.hpp
bridge.cpp
bridge-ode.cpp)

ADD_EXECUTABLE(cbridge4 main.cpp)
TARGET_LINK_LIBRARIES(cbridge4 bridge4)
TARGET_LINK_YOCTO(cbridge4 mk lua++)

INSTALL( TARGETS cbridge4 RUNTIME DESTINATION bin )
