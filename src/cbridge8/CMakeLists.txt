ADD_LIBRARY(cbridge8lib STATIC
lens.hpp
lens.cpp
vlens.hpp
vlens.cpp
)

ADD_EXECUTABLE(cbridge8 cbridge8.cpp)
TARGET_LINK_LIBRARIES(cbridge8 cbridge8lib)
TARGET_LINK_YOCTO(cbridge8 mk)

INSTALL( TARGETS cbridge8 RUNTIME DESTINATION bin )


