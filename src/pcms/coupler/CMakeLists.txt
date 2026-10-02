set(PCMS_COUPLER_HEADERS
        coupler_types.h
        field_communicator.hpp
        coupler.hpp
        distributed_evaluation_channel.h
        distributed_evaluation_partition.h
        distributed_point_evaluator.h
        field_serializer.h
        field_layout_communicator.h
        field_exchange_planner.h
        partition.h
        overlap_mask.h
        global_communicator.h
)


set(PCMS_COUPLER_SOURCES
        coupler.cpp
        field_layout_communicator.cpp
        field_exchange_planner.cpp
)

if(PCMS_ENABLE_XGC)
    list(APPEND PCMS_COUPLER_HEADERS
            serializer/xgc.h)
endif()

add_library(pcms_coupler ${PCMS_COUPLER_SOURCES})
set_target_properties(pcms_coupler PROPERTIES
        OUTPUT_NAME pcmscoupler
        EXPORT_NAME coupler)
target_sources(pcms_coupler PUBLIC
        FILE_SET coupler
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_COUPLER_HEADERS})
add_library(pcms::coupler ALIAS pcms_coupler)
target_compile_features(pcms_coupler PUBLIC cxx_std_20)

target_link_libraries(pcms_coupler PUBLIC
        redev::redev
        MPI::MPI_CXX
        Kokkos::kokkos
        perfstubs
        pcms::utility
        pcms::field)

if(PCMS_ENABLE_OMEGA_H)
    target_link_libraries(pcms_coupler PUBLIC Omega_h::omega_h)
endif()

if (PCMS_ENABLE_MESHFIELDS)
    target_link_libraries(pcms_coupler PUBLIC meshfields::meshfields)
endif ()

target_include_directories(pcms_coupler INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
        TARGETS pcms_coupler
        EXPORT pcms_coupler-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/coupler
        FILE_SET coupler DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(
        EXPORT pcms_coupler-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
