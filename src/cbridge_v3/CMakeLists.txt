ADD_LIBRARY(bridge_v3 STATIC
cbridge.hpp
cbridge.cpp
cbridge-ode.cpp
)

ADD_EXECUTABLE(cbridge_v3 main.cpp)
TARGET_LINK_LIBRARIES(cbridge_v3 bridge_v3)
TARGET_LINK_YOCTO(cbridge_v3 mk lua++)

INSTALL( TARGETS cbridge_v3 RUNTIME DESTINATION bin )

