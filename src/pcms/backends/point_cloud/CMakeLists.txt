# Point-cloud backend (no external dependency). No public headers.

# discretization: part of pcms::discretization
target_sources(pcms_discretization PRIVATE discretization.cpp
                                           discretization.hpp)

# layout and field data: part of pcms::field
target_sources(pcms_field PRIVATE layout.cpp field_data.cpp layout.hpp
                                  field_data.hpp)

# MLS evaluation and localization use Omega_h arrays (as before)
if(PCMS_ENABLE_OMEGA_H)
  target_sources(pcms_field PRIVATE localization.hpp mls_evaluator.hpp
                                    evaluator.hpp)
endif()
