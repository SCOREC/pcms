add_executable(cpl cpl.cc)
target_link_libraries(cpl coupler)

add_executable(test_init test_init.cc)
target_link_libraries(test_init coupler)

add_executable(rinter rinter.cc)
target_link_libraries(rinter coupler)

set(sources cpl test_init)
install(TARGETS ${sources} coupler)
