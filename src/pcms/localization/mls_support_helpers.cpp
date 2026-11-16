#include "pcms/localization/mls_support_helpers.h"
#include "pcms/localization/adj_search.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/print.h"

#include <Kokkos_Core.hpp>

namespace pcms
{

// replace with Kokkos::minmax_element when out of experimental
// https://kokkos.org/kokkos-core-wiki/API/algorithms/std-algorithms/all/StdMinMaxElement.html
void minmax(Omega_h::Read<Omega_h::LO> num_supports,
            unsigned& min_supports_found, unsigned& max_supports_found)
{
  using minMaxReducerType = Kokkos::MinMax<unsigned, Kokkos::HostSpace>;
  using minMaxValueType = minMaxReducerType::value_type;
  minMaxValueType minmax_val;
  Kokkos::parallel_reduce(
    num_supports.size(),
    KOKKOS_LAMBDA(int i, minMaxValueType& update) {
      if (static_cast<unsigned>(num_supports[i]) < update.min_val)
        update.min_val = static_cast<unsigned>(num_supports[i]);
      if (static_cast<unsigned>(num_supports[i]) > update.max_val)
        update.max_val = static_cast<unsigned>(num_supports[i]);
    },
    minMaxReducerType(minmax_val));
  Kokkos::fence();
  min_supports_found = minmax_val.min_val;
  max_supports_found = minmax_val.max_val;
}

void adapt_radii(unsigned min_req_supports, unsigned max_allowed_supports,
                 Omega_h::LO n_targets, Omega_h::Write<Omega_h::Real> radii2_l,
                 Omega_h::Write<Omega_h::LO> num_supports)
{
  Omega_h::parallel_for(
    "increase radius", n_targets, OMEGA_H_LAMBDA(const int& i) {
      Omega_h::LO nsupports = num_supports[i];
      if (nsupports < static_cast<Omega_h::LO>(min_req_supports)) {
        double factor =
          Omega_h::Real(min_req_supports) / Omega_h::Real(nsupports);
        OMEGA_H_CHECK_PRINTF(factor > 1.0,
                             "Factor should be more than 1.0: %f\n", factor);
        factor = (nsupports == 0 || factor > 1.5) ? 1.5 : factor;
        radii2_l[i] *= factor;
      } else if (nsupports > static_cast<Omega_h::LO>(max_allowed_supports)) {
        double factor =
          Omega_h::Real(min_req_supports) / Omega_h::Real(nsupports);
        OMEGA_H_CHECK_PRINTF(factor < 1.0,
                             "Factor should be less than 1.0: %f\n", factor);
        factor = (factor < 0.1) ? 0.33 : factor;
        radii2_l[i] *= factor;
      }
      num_supports[i] = 0;
    });
  Kokkos::fence();
}

// ---- point-cloud N² support search helpers --------------------------------

namespace
{

KOKKOS_INLINE_FUNCTION
Omega_h::Vector<3> load_point(const Omega_h::Reals& coords, int id, int dim)
{
  Omega_h::Vector<3> p{0, 0, 0};
  for (int d = 0; d < dim; ++d)
    p[d] = coords[id * dim + d];
  return p;
}

struct NSquareCountFunctor
{
  const int dim;
  const Omega_h::LO n_sources;
  const Omega_h::Reals target_coords;
  const Omega_h::Reals source_coords;
  const Omega_h::Write<Omega_h::Real> radii2;
  const Omega_h::Write<Omega_h::LO> num_supports;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int target_id) const
  {
    auto tgt = load_point(target_coords, target_id, dim);
    auto r2 = radii2[target_id];
    int count = 0;
    for (int s = 0; s < n_sources; ++s) {
      auto src = load_point(source_coords, s, dim);
      double dist2 = 0;
      for (int d = 0; d < dim; ++d) {
        double diff = tgt[d] - src[d];
        dist2 += diff * diff;
      }
      if (dist2 <= r2)
        ++count;
    }
    num_supports[target_id] = count;
  }
};

struct ScanSupportPtrFunctor
{
  Omega_h::Write<Omega_h::LO> support_ptr;
  Omega_h::Write<Omega_h::LO> num_supports;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, unsigned& update, const bool final) const
  {
    update += num_supports[i];
    if (final)
      support_ptr[i + 1] = update;
  }
};

struct FillSupportIdxFunctor
{
  const int dim;
  const Omega_h::LO n_sources;
  const Omega_h::Write<Omega_h::LO> support_ptr;
  const Omega_h::Write<Omega_h::LO> supports_idx;
  const Omega_h::Reals target_coords;
  const Omega_h::Reals source_coords;
  const Omega_h::Write<Omega_h::Real> radii2;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int target_id) const
  {
    auto tgt = load_point(target_coords, target_id, dim);
    auto r2 = radii2[target_id];
    auto pos = support_ptr[target_id];
    const auto end = support_ptr[target_id + 1];
    for (int s = 0; s < n_sources; ++s) {
      auto src = load_point(source_coords, s, dim);
      double dist2 = 0;
      for (int d = 0; d < dim; ++d) {
        double diff = tgt[d] - src[d];
        dist2 += diff * diff;
      }
      if (dist2 <= r2) {
        OMEGA_H_CHECK_PRINTF(pos < end,
                             "Support index out of bounds: pos %d end %d "
                             "target_id %d\n",
                             pos, end, target_id);
        supports_idx[pos++] = s;
      }
    }
  }
};

} // anonymous namespace

SupportResults BuildPointCloudSupports(const Omega_h::Reals& source_coords,
                                       const Omega_h::Reals& target_coords,
                                       int dim, double radius,
                                       unsigned min_req_supports,
                                       bool adapt_radius_flag,
                                       unsigned max_iterations)
{
  PCMS_ALWAYS_ASSERT(radius >= 0.0);
  const Omega_h::LO n_targets = target_coords.size() / dim;
  const Omega_h::LO n_sources = source_coords.size() / dim;
  const unsigned max_allowed = 3 * min_req_supports;
  const auto radius_sq = static_cast<Omega_h::Real>(radius * radius);

  auto radii2 = Omega_h::Write<Omega_h::Real>(n_targets, radius_sq);
  auto num_supports = Omega_h::Write<Omega_h::LO>(n_targets, 0);

  unsigned min_found = 0, max_found = 0;
  unsigned loop_count = 0;

  while (!within_number_of_support_range(min_found, max_found, min_req_supports,
                                         max_allowed)) {
    Kokkos::parallel_for("n2_count", n_targets,
                         NSquareCountFunctor{dim, n_sources, target_coords,
                                             source_coords, radii2,
                                             num_supports});
    Kokkos::fence();

    ++loop_count;
    if (!adapt_radius_flag)
      break;
    if (loop_count > max_iterations) {
      pcms::printError(
        "BuildPointCloudSupports: radius adjustment did not converge after "
        "%d iterations.\n",
        max_iterations);
      break;
    }

    minmax(Omega_h::read(num_supports), min_found, max_found);

    if (!within_number_of_support_range(min_found, max_found, min_req_supports,
                                        max_allowed)) {
      pcms::printInfo("BuildPointCloudSupports: adjusting radius iter %d "
                      "(min=%d max=%d req=%d allowed=%d)\n",
                      loop_count, min_found, max_found, min_req_supports,
                      max_allowed);
      adapt_radii(min_req_supports, max_allowed, n_targets, radii2,
                  num_supports);
    }
  }

  pcms::printInfo(
    "BuildPointCloudSupports: support search done after %d iterations "
    "(min=%d max=%d)\n",
    loop_count, min_found, max_found);

  // Build CSR offset array
  auto support_ptr = Omega_h::Write<Omega_h::LO>(n_targets + 1, 0);
  unsigned total_supports = 0;
  Kokkos::parallel_scan("scan_support_ptr", n_targets,
                        ScanSupportPtrFunctor{support_ptr, num_supports},
                        total_supports);
  Kokkos::fence();

  auto supports_idx = Omega_h::Write<Omega_h::LO>(total_supports, 0);
  Kokkos::parallel_for("fill_support_idx", n_targets,
                       FillSupportIdxFunctor{dim, n_sources, support_ptr,
                                             supports_idx, target_coords,
                                             source_coords, radii2});
  Kokkos::fence();

  return SupportResults{Omega_h::LOs(support_ptr), Omega_h::LOs(supports_idx),
                        radii2};
}

} // namespace pcms
