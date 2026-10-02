if(PCMS_ENABLE_XGC AND PCMS_ENABLE_OMEGA_H)
  add_executable(osh2XgcRC XgcRCfromOsh.cpp)
  target_link_libraries(osh2XgcRC pcms::core Omega_h::omega_h)
  add_executable(cpn_gen CpnFromOsh.cpp)
  target_link_libraries(cpn_gen pcms::core Omega_h::omega_h)
endif()

add_executable(distributed_field_eval_server distributed_field_eval_server.cpp)
target_link_libraries(distributed_field_eval_server pcms::coupler)

