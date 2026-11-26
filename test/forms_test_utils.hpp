#ifndef PCMS_TEST_FORMS_TEST_UTILS_HPP
#define PCMS_TEST_FORMS_TEST_UTILS_HPP

#include <Kokkos_Core.hpp>
#include <catch2/catch_test_macros.hpp>
#include <pcms/field/point_rows.hpp>
#include <pcms/transfer/forms.hpp>
#include <map>
#include <utility>

namespace pcms::test
{

using Matrix = std::map<std::pair<GO, GO>, Real>;

/// @return E^T W E of the form's quadrature from the test evaluator's rows,
///         keyed by the global ids of the DOF holders
inline Matrix AssembleMass(const transfer::detail::BilinearForm& form)
{
  const auto* source = detail::AsRowsSource(form.TestEvaluator());
  REQUIRE(source != nullptr);
  const auto& rows = source->Rows();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.offsets);
  const auto holders =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.holders);
  const auto values =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.values);
  const auto device_gids = source->HolderGids();
  Kokkos::View<GO*, DeviceMemorySpace> gids_view("gids", device_gids.size());
  Kokkos::parallel_for(
    "gids",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<LO>(gids_view.extent(0))),
    KOKKOS_LAMBDA(const LO i) { gids_view(i) = device_gids[i]; });
  const auto gids =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, gids_view);
  const auto weights = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace{}, form.GetQuadrature().weights);
  Matrix matrix;
  for (std::size_t p = 0; p + 1 < offsets.extent(0); ++p) {
    for (LO i = offsets(p); i < offsets(p + 1); ++i) {
      for (LO j = offsets(p); j < offsets(p + 1); ++j) {
        matrix[{gids(holders(i)), gids(holders(j))}] +=
          weights(p) * values(i) * values(j);
      }
    }
  }
  return matrix;
}

} // namespace pcms::test

#endif // PCMS_TEST_FORMS_TEST_UTILS_HPP
