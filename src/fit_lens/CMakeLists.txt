ADD_EXECUTABLE(fit_lens fit_lens.cpp)
TARGET_LINK_YOCTO(fit_lens mk)

ADD_EXECUTABLE(sim_lens sim_lens.cpp)
TARGET_LINK_YOCTO(sim_lens mk)


