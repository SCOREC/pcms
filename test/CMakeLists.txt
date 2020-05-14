add_executable(cpl cpl.cc)
target_link_libraries(cpl coupler)

add_executable(test_init test_init.cc)
target_link_libraries(test_init coupler)

add_executable(driver driver.cc)
target_link_libraries(driver coupler)

set(drivers cpl test_init driver)
install(TARGETS ${drivers} DESTINATION ${CMAKE_BINARY_DIR}/bin)

if(BUILD_TESTING)
  include(testing.cmake)
endif()
