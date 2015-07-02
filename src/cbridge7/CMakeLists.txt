ADD_LIBRARY(cbridge STATIC 
lens.hpp 
lens.cpp 
bridge.hpp 
bridge.cpp 
bridge-ode.cpp
bridge-alpha.cpp
bridge-hmax.cpp
bridge-theta.cpp
bridge-data.cpp
bridge-opt.cpp
datafile.hpp
datafile.cpp
window.hpp
window.cpp
bridge-mt.cpp
)

ADD_EXECUTABLE(cbridge7 cbridge7.cpp)
TARGET_LINK_LIBRARIES(cbridge7 cbridge)

TARGET_LINK_YOCTO(cbridge7 mk)
INSTALL( TARGETS cbridge7 RUNTIME DESTINATION bin )

ADD_EXECUTABLE(abacus abacus.cpp)
TARGET_LINK_LIBRARIES(abacus cbridge)
TARGET_LINK_YOCTO(abacus mk)

INSTALL( TARGETS abacus RUNTIME DESTINATION bin )

