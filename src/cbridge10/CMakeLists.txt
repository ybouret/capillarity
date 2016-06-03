ADD_EXECUTABLE(cbridge main.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(cbridge mk)

ADD_EXECUTABLE(abacus abacus.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(abacus mk)

ADD_EXECUTABLE(limits limits.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(limits mk)

INSTALL( TARGETS cbridge abacus limits RUNTIME DESTINATION bin )

