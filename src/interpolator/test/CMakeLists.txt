cmake_minimum_required(VERSION 3.20)

project(test_1 LANGUAGES CXX)

find_program(BASH_EXECUTABLE NAMES bash REQUIRED)

add_executable(test_polybasis test_polybasis.cpp)

target_link_libraries(test_polybasis PUBLIC pcms::core)

target_include_directories(test_polybasis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()

add_test(
    NAME test_polybasis
    COMMAND $<TARGET_FILE:cpp_test>
    )
