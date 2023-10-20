#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../capi/client.h"
#include "pcms.h"
#include "pcms/xgc_field_adapter.h"
#include <variant>
#include <redev_variant_tools.h>
// #ifdef PCMS_HAS_OMEGA_H
//   #include "pcms/omega_h_field.h"
// #endif
#include <fstream>
#include "pcms/xgc_reverse_classification.h"
#include "pcms/dummy_field_adapter.h"
#include <type_traits>

template <typename T>
struct VALID_PCMS_TYPE : std::false_type {};
template<>
struct VALID_PCMS_TYPE<double> : std::true_type {};
template<>
struct VALID_PCMS_TYPE<float> : std::true_type {};
template<>
struct VALID_PCMS_TYPE<int> : std::true_type {};
template<>
struct VALID_PCMS_TYPE<long int> : std::true_type {};



namespace py = pybind11;


template <typename T>
PcmsFieldAdapterHandle xgc_field_adapter_numpy(const char* name, MPI_Comm comm, py::array_t<T> data, 
const PcmsReverseClassificationHandle rc, in_overlap_function in_overlap) {
    //Check T typing (SFINAE)
    static_assert(VALID_PCMS_TYPE<T>::value, "Passed type is invalid");

    auto buffer_info = data.request(); // Get information about the NumPy array

    int size = buffer_info.size; // Get the size of the NumPy array

    auto field_adapter = new pcms::FieldAdapterVariant{};
    PCMS_ALWAYS_ASSERT(rc.pointer != nullptr);

    auto reverse_classification = reinterpret_cast<const pcms::ReverseClassificationVertex*>(rc.pointer);
    PCMS_ALWAYS_ASSERT(reverse_classification != nullptr);

    pcms_create_xgc_field_adapter_t<T>(name, comm, buffer_info.ptr, size, 
    *reverse_classification, in_overlap, *field_adapter);

    return PcmsFieldAdapterHandle{reinterpret_cast<void*>(field_adapter)};
}


PYBIND11_MODULE(pcms_client, m)
{
    m.doc() = "Python bindings for C PCMS client functions";

    py::enum_<PcmsAdapterType>(m, "PcmsAdapterType")
        .value("PCMS_ADAPTER_XGC", PCMS_ADAPTER_XGC)
        .value("PCMS_ADAPTER_OMEGAH", PCMS_ADAPTER_OMEGAH)
        .value("PCMS_ADAPTER_GENE", PCMS_ADAPTER_GENE)
        .value("PCMS_ADAPTER_GEM", PCMS_ADAPTER_GEM);

    py::enum_<PcmsType>(m, "PcmsType")
        .value("PCMS_FLOAT", PCMS_FLOAT)
        .value("PCMS_DOUBLE", PCMS_DOUBLE)
        .value("PCMS_INT", PCMS_INT)
        .value("PCMS_LONG_INT", PCMS_LONG_INT);


    m.def("pcms_create_client", &pcms_create_client, "Create a PCMS client");
    m.def("pcms_destroy_client", &pcms_destroy_client, "Destroy a PCMS client");
    m.def("pcms_load_reverse_classification", &pcms_load_reverse_classification, "Load reverse classification data");
    m.def("pcms_destroy_reverse_classification", &pcms_destroy_reverse_classification, "Destroy reverse classification data");
    m.def("pcms_reverse_classification_count_verts", &pcms_reverse_classification_count_verts, "Count vertices in reverse classification");

    m.def("pcms_create_xgc_field_adapter",&xgc_field_adapter_numpy<double>, "Create an XGC field adapter for double");
    m.def("pcms_create_xgc_field_adapter",&xgc_field_adapter_numpy<float>, "Create an XGC field adapter for float");
    m.def("pcms_create_xgc_field_adapter",&xgc_field_adapter_numpy<int>, "Create an XGC field adapter for int");
    m.def("pcms_create_xgc_field_adapter",&xgc_field_adapter_numpy<long int>, "Create an XGC field adapter for long int");

    m.def("pcms_create_xgc_field_adapter", &pcms_create_xgc_field_adapter, "Create an XGC field adapter");
    m.def("pcms_create_dummy_field_adapter", &pcms_create_dummy_field_adapter, "Create a dummy field adapter");
    m.def("pcms_destroy_field_adapter", &pcms_destroy_field_adapter, "Destroy a field adapter");
    m.def("pcms_add_field", &pcms_add_field, "Add a field");
    m.def("pcms_send_field_name", &pcms_send_field_name, "Send field name");
    m.def("pcms_receive_field_name", &pcms_receive_field_name, "Receive field name");
    m.def("pcms_send_field", &pcms_send_field, "Send field");
    m.def("pcms_receive_field", &pcms_receive_field, "Receive field");
    m.def("pcms_begin_send_phase", &pcms_begin_send_phase, "Begin send phase");
    m.def("pcms_end_send_phase", &pcms_end_send_phase, "End send phase");
    m.def("pcms_begin_receive_phase", &pcms_begin_receive_phase, "Begin receive phase");
    m.def("pcms_end_receive_phase", &pcms_end_receive_phase, "End receive phase");
}


