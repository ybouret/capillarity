ADD_LIBRARY(cbridge STATIC lens.hpp lens.cpp bridge.hpp bridge.cpp)
ADD_EXECUTABLE(cbridge7 cbridge7.cpp)
TARGET_LINK_LIBRARIES(cbridge7 cbridge)
TARGET_LINK_YOCTO(cbridge7 mk)
INSTALL( TARGETS cbridge7 RUNTIME DESTINATION bin )


