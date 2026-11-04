# Omega_h + MeshFields bridge: the MeshFields field adapter, part of
# pcms::field. Built only when both backends are enabled. No public headers.
target_sources(pcms_field PRIVATE layout.cpp layout.hpp evaluator_backend.hpp
                                  evaluator.hpp field_data.hpp)
target_link_libraries(pcms_field PUBLIC meshfields::meshfields)
