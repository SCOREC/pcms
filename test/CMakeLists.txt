add_executable(passthrough passthrough.cc)
target_link_libraries(passthrough coupler)

add_executable(test_init test_init.cc)
target_link_libraries(test_init coupler)

add_executable(cpl cpl.cc)
target_link_libraries(cpl coupler)

set(drivers passthrough test_init cpl)
install(TARGETS ${drivers})

#if(BUILD_TESTING)
#  include(testing.cmake)
#endif()
