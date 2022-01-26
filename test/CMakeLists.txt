add_executable(passthrough passthrough.cc)
target_link_libraries(passthrough coupler)

add_executable(test_init test_init.cc)
target_link_libraries(test_init coupler)

add_executable(cpl cpl.cc)
add_executable(gemxgc gemxgc.cc)

target_link_libraries(cpl coupler)
target_link_libraries(gemxgc coupler)

set(drivers passthrough test_init cpl gemxgc)
install(TARGETS ${drivers})

if(BUILD_TESTING)
  include(testing.cmake)
endif()
