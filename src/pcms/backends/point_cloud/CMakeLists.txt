# Point-cloud backend (no external dependency). No public headers.

# layout and field data: part of pcms::field
target_sources(pcms_field PRIVATE layout.cpp field_data.cpp layout.hpp
                                  field_data.hpp)

# MLS evaluation and the geometry use Omega_h arrays (as before)
if(PCMS_ENABLE_OMEGA_H)
  target_sources(pcms_field PRIVATE mls_evaluator.hpp evaluator.hpp
                                    geometry.hpp)
endif()
