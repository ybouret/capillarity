ADD_LIBRARY(bridge4 STATIC
bridge.hpp
bridge.cpp
bridge-ode.cpp
bridge-theta.cpp
bridge-params.cpp
bridge-alpha.cpp
application.hpp
application.cpp
application_v2.cpp
)

ADD_EXECUTABLE(cbridge4 main.cpp)
TARGET_LINK_LIBRARIES(cbridge4 bridge4)
TARGET_LINK_YOCTO(cbridge4 mk lua++ lang)

ADD_EXECUTABLE(find-theta4 find-theta.cpp)
TARGET_LINK_LIBRARIES(find-theta4 bridge4)
TARGET_LINK_YOCTO(find-theta4 mk lua++)

ADD_EXECUTABLE(find-alpha4 find-alpha.cpp)
TARGET_LINK_LIBRARIES(find-alpha4 bridge4)
TARGET_LINK_YOCTO(find-alpha4 mk lua++)

ADD_EXECUTABLE(abacus4 abacus.cpp)
TARGET_LINK_LIBRARIES(abacus4 bridge4)
TARGET_LINK_YOCTO(abacus4 mk lua++)


ADD_EXECUTABLE(theta-of-alpha theta-of-alpha.cpp)
TARGET_LINK_LIBRARIES(theta-of-alpha bridge4)
TARGET_LINK_YOCTO(theta-of-alpha mk lua++)

INSTALL( TARGETS 
cbridge4
find-theta4 
find-alpha4
abacus4
theta-of-alpha
RUNTIME DESTINATION bin )
