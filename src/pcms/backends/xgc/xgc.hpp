#ifndef PCMS_BACKENDS_XGC_XGC_HPP
#define PCMS_BACKENDS_XGC_XGC_HPP

// Public entry points of the XGC backend: wrap XGC's own field arrays.

#include "pcms/backends/xgc/reverse_classification.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_factory.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

#include <mpi.h>

#include <memory>
#include <string>

namespace pcms::backend::xgc
{

/// Field factory for XGC fields on the plane described by
/// reverse_classification with num_plane_nodes nodes. layout_name identifies
/// the layout in the coupler.
std::shared_ptr<FieldFactory> CreateFieldFactory(
  const ReverseClassificationVertex& reverse_classification, LO num_plane_nodes,
  std::string layout_name = "");

/// Field named name whose values live in XGC's array data (not copied).
/// factory must come from CreateFieldFactory. plane_comm, if given, is the
/// communicator of the ranks sharing the plane. T is one of int8_t, int32_t,
/// int64_t, float, double.
template <typename T>
Field<T> CreateField(const FieldFactory& factory, std::string name,
                     Rank1View<T, HostMemorySpace> data,
                     MPI_Comm plane_comm = MPI_COMM_NULL,
                     ValueBasis basis = {});

} // namespace pcms::backend::xgc

#endif // PCMS_BACKENDS_XGC_XGC_HPP
