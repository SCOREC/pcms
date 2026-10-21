if(PCMS_ENABLE_XGC AND PCMS_ENABLE_OMEGA_H)
  add_executable(osh2XgcRC XgcRCfromOsh.cpp)
  target_link_libraries(osh2XgcRC pcms::core Omega_h::omega_h)
  add_executable(cpn_gen CpnFromOsh.cpp)
  target_link_libraries(cpn_gen pcms::core Omega_h::omega_h)
endif()

if(PCMS_ENABLE_GNET)
  add_executable(gnet2gtc gnet/gnet2gtc.cpp)
  target_link_libraries(gnet2gtc pcms::gnet)
  # let the installed tool find the pcms and dependency libraries
  set_target_properties(gnet2gtc PROPERTIES
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE)
  install(TARGETS gnet2gtc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
