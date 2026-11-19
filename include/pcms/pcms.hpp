#ifndef PCMS_PCMS_HPP
#define PCMS_PCMS_HPP

// Umbrella header: includes the public PCMS API that is enabled in this
// build. Backend adapters (pcms/backends/...) are not included; include the
// backend header you use explicitly. pcms/version.hpp is not included (it
// defines a non-inline function).

#include "pcms/configuration.hpp"

// utility
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/entity_types.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/print.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

// geometry
#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/supports.hpp"

// field
#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/differential_operator.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_factory.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/field/gradient_source.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/query_points.hpp"
#include "pcms/field/value_view.hpp"

// transfer
#include "pcms/transfer/copy.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/projection.hpp"
#include "pcms/transfer/transfer_operator.hpp"

// coupler
#ifdef PCMS_ENABLE_COUPLER
#include "pcms/coupler/coupler.hpp"
#include "pcms/coupler/coupler_types.hpp"
#include "pcms/coupler/overlap_mask.hpp"
#endif

// transient
#include "pcms/transient/timestepper.hpp"

#endif // PCMS_PCMS_HPP
