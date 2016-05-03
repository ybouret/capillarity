ADD_EXECUTABLE(profile profile.cpp)
TARGET_LINK_YOCTO(profile mk)

ADD_EXECUTABLE(bridge bridge.hpp bridge.cpp lens.cpp lens.hpp main.cpp)
TARGET_LINK_YOCTO(bridge mk)

ADD_EXECUTABLE(find_lens find_lens.cpp)
TARGET_LINK_YOCTO(find_lens mk)

ADD_EXECUTABLE(real_lens real_lens.cpp)
TARGET_LINK_YOCTO(real_lens mk)

ADD_EXECUTABLE(fit_lens fit_lens.cpp)
TARGET_LINK_YOCTO(fit_lens mk)

INSTALL( TARGETS profile bridge fit_lens RUNTIME DESTINATION bin )
