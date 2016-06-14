ADD_LIBRARY(bridge STATIC bridge.hpp bridge.cpp)

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


INSTALL( TARGETS cbridge abacus heuristic extheta RUNTIME DESTINATION bin )

