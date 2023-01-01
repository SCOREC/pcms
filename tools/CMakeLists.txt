if(WDMCPL_ENABLE_XGC AND WDMCPL_ENABLE_OMEGA_H)
    add_executable(osh2XgcRC XgcRCfromOsh.cpp)
    target_link_libraries(osh2XgcRC wdmcpl::wdmcpl Omega_h::omega_h)
endif()