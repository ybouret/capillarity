ADD_LIBRARY(bridge STATIC 
bridge.hpp bridge.cpp 
application.hpp application.cpp
)

IF(FALSE)
ADD_EXECUTABLE(cbridge main.cpp)
TARGET_LINK_YOCTO(cbridge mk)
TARGET_LINK_LIBRARIES(cbridge bridge)

ADD_EXECUTABLE(abacus abacus.cpp setup.hpp setup.cpp)
TARGET_LINK_YOCTO(abacus mk)
TARGET_LINK_LIBRARIES(abacus bridge)

ADD_EXECUTABLE(heuristic heuristic.cpp)
TARGET_LINK_YOCTO(heuristic mk)
TARGET_LINK_LIBRARIES(heuristic bridge)

ADD_EXECUTABLE(extheta extheta.cpp setup.hpp setup.cpp)
TARGET_LINK_YOCTO(extheta mk)
TARGET_LINK_LIBRARIES(extheta bridge)

ADD_EXECUTABLE(extheta_raw extheta_raw.cpp setup.hpp setup.cpp)
TARGET_LINK_YOCTO(extheta_raw mk seem)
TARGET_LINK_LIBRARIES(extheta_raw bridge)

INSTALL( TARGETS cbridge abacus heuristic extheta extheta_raw RUNTIME DESTINATION bin )
ENDIF()


ADD_EXECUTABLE(compute_theta_v1 compute_theta_v1.cpp)
TARGET_LINK_YOCTO(compute_theta_v1 mk seem)
TARGET_LINK_LIBRARIES(compute_theta_v1 bridge)

INSTALL( TARGETS compute_theta_v1 RUNTIME DESTINATION bin )
