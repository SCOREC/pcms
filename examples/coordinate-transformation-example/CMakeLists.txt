add_executable(coordinate_transformations coordinate_transformations.cpp)
target_link_libraries(coordinate_transformations PRIVATE pcms::field
                                                         pcms::transfer)
