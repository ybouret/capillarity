ADD_EXECUTABLE(fit_lens fit_lens.cpp)
TARGET_LINK_YOCTO(fit_lens mk)
INSTALL( TARGETS fit_lens RUNTIME DESTINATION bin )


#ADD_EXECUTABLE(sim_lens sim_lens.cpp)
#TARGET_LINK_YOCTO(sim_lens mk)


