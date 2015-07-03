ADD_LIBRARY(cbridge8lib STATIC
lens.hpp
lens.cpp
vlens.hpp
vlens.cpp

bridge.hpp
bridge.cpp
bridge-ode.cpp
bridge-alpha.cpp
bridge-hmax.cpp
bridge-theta.cpp
)

ADD_EXECUTABLE(cbridge8 cbridge8.cpp)
TARGET_LINK_LIBRARIES(cbridge8 cbridge8lib)
TARGET_LINK_YOCTO(cbridge8 mk)

INSTALL( TARGETS cbridge8 RUNTIME DESTINATION bin )

ADD_EXECUTABLE(abacus8 abacus8.cpp)
TARGET_LINK_LIBRARIES(abacus8 cbridge8lib)
TARGET_LINK_YOCTO(abacus8 mk)

INSTALL( TARGETS abacus8 RUNTIME DESTINATION bin )


