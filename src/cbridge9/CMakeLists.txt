ADD_EXECUTABLE(profile profile.cpp)
TARGET_LINK_YOCTO(profile mk)

ADD_EXECUTABLE(bridge bridge.hpp bridge.cpp lens.cpp lens.hpp main.cpp)
TARGET_LINK_YOCTO(bridge mk)

ADD_EXECUTABLE(fit_lens fit_lens.cpp)
TARGET_LINK_YOCTO(fit_lens mk)

ADD_EXECUTABLE(abaccus bridge.hpp bridge.cpp lens.cpp lens.hpp abaccus.cpp)
TARGET_LINK_YOCTO(abaccus mk)

ADD_EXECUTABLE(geom bridge.hpp bridge.cpp lens.cpp lens.hpp geom.cpp)
TARGET_LINK_YOCTO(geom mk)

INSTALL( TARGETS profile bridge fit_lens abaccus geom RUNTIME DESTINATION bin )
