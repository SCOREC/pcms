add_executable(driver driver.cc)
target_link_libraries(driver coupler)

install(TARGETS driver coupler)
