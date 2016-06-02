ADD_EXECUTABLE(cbridge main.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(cbridge mk)

ADD_EXECUTABLE(abacus abacus.cpp bridge.hpp bridge.cpp)
TARGET_LINK_YOCTO(abacus mk)

INSTALL( TARGETS cbridge abacus RUNTIME DESTINATION bin )

