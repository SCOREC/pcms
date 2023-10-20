find_package(pybind11)

pybind11_add_module(pcms_client client_pybind.cpp)
target_link_libraries()