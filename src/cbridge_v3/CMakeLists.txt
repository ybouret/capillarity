ADD_LIBRARY(bridge_v3 STATIC
cbridge.hpp
cbridge.cpp
cbridge-ode.cpp
cbridge-alpha.cpp
cbridge-theta.cpp
par-bridge.hpp
par-bridge.cpp
)

ADD_EXECUTABLE(cbridge_v3 main.cpp)
TARGET_LINK_LIBRARIES(cbridge_v3 bridge_v3)
TARGET_LINK_YOCTO(cbridge_v3 mk lua++)

ADD_EXECUTABLE(find-alpha find-alpha.cpp)
TARGET_LINK_LIBRARIES(find-alpha bridge_v3)
TARGET_LINK_YOCTO(find-alpha mk lua++)

ADD_EXECUTABLE(find-theta find-theta.cpp)
TARGET_LINK_LIBRARIES(find-theta bridge_v3)
TARGET_LINK_YOCTO(find-theta mk lua++)

ADD_EXECUTABLE(abacus abacus.cpp)
TARGET_LINK_LIBRARIES(abacus bridge_v3)
TARGET_LINK_YOCTO(abacus mk lua++)

INSTALL( TARGETS 
cbridge_v3 
find-alpha
find-theta
abacus
RUNTIME DESTINATION bin )

