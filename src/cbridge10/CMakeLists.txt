ADD_EXECUTABLE(cbridge main.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(cbridge mk)

INSTALL( TARGETS cbridge RUNTIME DESTINATION bin )

