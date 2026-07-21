#ifndef PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H
#define PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H

#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/localization/point_search.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"

#include <Omega_h_mesh.hpp>
#include <Kokkos_Core.hpp>
#include <memory>
#include <stdexcept>
#include <variant>

namespace pcms
{

// ---------------------------------------------------------------------------
// Localization hint — computed once per query point set, reused across Evaluate
// ---------------------------------------------------------------------------
struct OmegaHLagrangeLocHint
{
  Kokkos::View<LO*, DeviceMemorySpace>
    elem_ids; // containing element (valid pts)
  Kokkos::View<Real**, DeviceMemorySpace> bary;      // [n_valid x (dim+1)]
  Kokkos::View<LO*, DeviceMemorySpace> orig_indices; // original query index
  Kokkos::View<LO*, DeviceMemorySpace> missing_indices;
  OutOfBoundsMode mode;
};

namespace detail
{

template <int Dim>
struct CopyCoordsFunctor
{
  using DefaultLayout =
    detail::default_layout_for_memory_space_t<DeviceMemorySpace>;

  Kokkos::View<Real**, DeviceMemorySpace> coords_d;
  Rank2View<const Real, DeviceMemorySpace, DefaultLayout> raw_coords;

  CopyCoordsFunctor(
    Kokkos::View<Real**, DeviceMemorySpace> coords_d_,
    Rank2View<const Real, DeviceMemorySpace, DefaultLayout> raw_coords_)
    : coords_d(coords_d_), raw_coords(raw_coords_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i) const
  {
    for (int d = 0; d < Dim; ++d) {
      coords_d(i, d) = raw_coords(i, d);
    }
  }
};

template <int Dim>
OmegaHLagrangeLocHint BuildLagrangeLocHint(
  Omega_h::Mesh& mesh, int mesh_dim,
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results,
  Kokkos::View<Real**, DeviceMemorySpace> coords_d,
  Kokkos::View<LO*, DeviceMemorySpace> owning_ids, OutOfBoundsMode mode)
{
  LO n = static_cast<LO>(results.size());

  // First pass: count valid and missing
  Kokkos::View<int*, DeviceMemorySpace> is_valid("is_valid", n);
  Kokkos::parallel_for(
    "CheckValidity",
    Kokkos::RangePolicy<typename DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(LO i) {
      bool out = (owning_ids(i) < 0) || (results(i).element_id < 0);
      is_valid(i) = out ? 0 : 1;
    });

  // Exclusive scan to get compaction indices
  Kokkos::View<LO*, DeviceMemorySpace> valid_offsets("valid_offsets", n);
  LO nv;
  Kokkos::parallel_scan(
    "ScanValid",
    Kokkos::RangePolicy<typename DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i, LO& update, const bool final) {
      const LO val = is_valid(i);
      if (final && val)
        valid_offsets(i) = update;
      update += val;
    },
    nv);

  LO nm = n - nv;

  // Allocate output arrays
  Kokkos::View<LO*, DeviceMemorySpace> elem_ids("elem_ids", nv);
  Kokkos::View<Real**, DeviceMemorySpace> bary("bary", nv, mesh_dim + 1);
  Kokkos::View<LO*, DeviceMemorySpace> orig_indices("orig_indices", nv);
  Kokkos::View<LO*, DeviceMemorySpace> missing_indices("missing_indices", nm);

  // Compact valid and missing separately
  Kokkos::View<LO*, DeviceMemorySpace> missing_offsets("missing_offsets", n);
  Kokkos::parallel_scan(
    "ScanMissing",
    Kokkos::RangePolicy<typename DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i, LO& update, const bool final) {
      const LO val = 1 - is_valid(i);
      if (final && val)
        missing_offsets(i) = update;
      update += val;
    });

  auto elem_verts = mesh.ask_elem_verts();
  auto mesh_coords = mesh.coords();

  Kokkos::parallel_for(
    "CompactData",
    Kokkos::RangePolicy<typename DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(LO i) {
      if (is_valid(i)) {
        LO k = valid_offsets(i);
        LO owner_elem = owning_ids(i);
        elem_ids(k) = owner_elem;
        orig_indices(k) = i;

        // Recompute barycentric coordinates in the owning element
        // using global coordinates
        int nvpe = mesh_dim + 1;
        // Get vertices of the owning element
        Omega_h::Few<Omega_h::LO, 4> verts;
        for (int v = 0; v < nvpe; ++v) {
          verts[v] = elem_verts[owner_elem * nvpe + v];
        }

        // Get vertex coordinates
        if (mesh_dim == 2) {
          Omega_h::Few<Omega_h::Vector<2>, 3> vertex_coords;
          for (int v = 0; v < 3; ++v) {
            vertex_coords[v] = Omega_h::get_vector<2>(mesh_coords, verts[v]);
          }
          Omega_h::Vector<2> point;
          for (int d = 0; d < 2; ++d) {
            point[d] = coords_d(i, d);
          }
          auto local =
            Omega_h::barycentric_from_global<2, 2>(point, vertex_coords);
          for (int d = 0; d < 3; ++d) {
            bary(k, d) = local[d];
          }
        } else if (mesh_dim == 3) {
          Omega_h::Few<Omega_h::Vector<3>, 4> vertex_coords;
          for (int v = 0; v < 4; ++v) {
            vertex_coords[v] = Omega_h::get_vector<3>(mesh_coords, verts[v]);
          }
          Omega_h::Vector<3> point;
          for (int d = 0; d < 3; ++d) {
            point[d] = coords_d(i, d);
          }
          auto local =
            Omega_h::barycentric_from_global<3, 3>(point, vertex_coords);
          for (int d = 0; d < 4; ++d) {
            bary(k, d) = local[d];
          }
        }
      } else {
        LO k = missing_offsets(i);
        missing_indices(k) = i;
      }
    });

  return OmegaHLagrangeLocHint{elem_ids, bary, orig_indices, missing_indices,
                               mode};
}

inline std::variant<GridPointSearch2D, GridPointSearch3D> MakeSearch(
  Omega_h::Mesh& mesh)
{
  if (mesh.dim() == 2)
    return GridPointSearch2D(mesh, 10, 10);
  else if (mesh.dim() == 3)
    return GridPointSearch3D(mesh, 10, 10, 10);
  throw std::invalid_argument(
    "OmegaHLagrangeEvaluatorFactory: only 2D and 3D meshes are supported");
}

} // namespace detail

// ---------------------------------------------------------------------------
// PointEvaluator
// ---------------------------------------------------------------------------

// OmegaHLagrangePointEvaluator<T> implements PointEvaluator<T> for simplex
// meshes backed by Omega_h. Localization results (element IDs, barycentric
// coordinates) are computed once at construction and cached for repeated
// Evaluate calls.
//
// Output shape: [num_query_points][num_components].
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class OmegaHLagrangePointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  OmegaHLagrangePointEvaluator(
    std::shared_ptr<const OmegaHLagrangeLayout> layout,
    OmegaHLagrangeLocHint hint, Real fill_value)
    : layout_(std::move(layout)),
      hint_(std::move(hint)),
      fill_value_(fill_value)
  {
  }

  void Evaluate(
    const Field<T>& field,
    Rank2View<T, DeviceMemorySpace, LayoutPolicy> values) const override
  {
    PCMS_FUNCTION_TIMER;
    auto dof_data = field.GetDOFHolderData();
    LO n_valid = static_cast<LO>(hint_.elem_ids.size());
    int n_comp = layout_->GetNumComponents();

    PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(n_comp));

    if (layout_->GetOrder() == 0) {
      Kokkos::parallel_for(
        "OmegaHLagrangePointEvaluator::EvaluateOrder0",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_valid),
        KOKKOS_CLASS_LAMBDA(LO k) {
          LO orig = hint_.orig_indices(k);
          LO elem = hint_.elem_ids(k);
          for (int c = 0; c < n_comp; ++c) {
            values(orig, c) = dof_data(elem, c);
          }
        });

    } else {
      // Order-1: barycentric interpolation over element vertices
      Omega_h::Mesh& mesh = const_cast<Omega_h::Mesh&>(layout_->GetMesh());
      int mesh_dim = mesh.dim();
      int nvpe = mesh_dim + 1;
      auto elem_verts = mesh.ask_elem_verts();

      Kokkos::parallel_for(
        "OmegaHLagrangePointEvaluator::EvaluateOrder1",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_valid),
        KOKKOS_CLASS_LAMBDA(LO k) {
          LO orig = hint_.orig_indices(k);
          LO elem = hint_.elem_ids(k);
          for (int c = 0; c < n_comp; ++c) {
            T val = T{};
            for (int v = 0; v < nvpe; ++v) {
              LO vert = elem_verts[elem * nvpe + v];
              val += static_cast<T>(hint_.bary(k, v)) * dof_data(vert, c);
            }
            values(orig, c) = val;
          }
        });
    }

    if (hint_.mode == OutOfBoundsMode::FILL) {
      T fill_val = static_cast<T>(fill_value_);
      Kokkos::parallel_for(
        "OmegaHLagrangePointEvaluator::FillOutOfBounds",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
          0, static_cast<LO>(hint_.missing_indices.size())),
        KOKKOS_CLASS_LAMBDA(LO k) {
          LO orig = hint_.missing_indices(k);
          for (int c = 0; c < n_comp; ++c) {
            values(orig, c) = fill_val;
          }
        });
    }
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
  OmegaHLagrangeLocHint hint_;
  Real fill_value_;
};

// ---------------------------------------------------------------------------
// EvaluatorFactory
// ---------------------------------------------------------------------------

// OmegaHLagrangeEvaluatorFactory<T> implements FieldEvaluatorFactory<T> for
// simplex meshes backed by Omega_h. It owns the spatial search structure and
// creates OmegaHLagrangePointEvaluator instances on demand.
template <typename T>
class OmegaHLagrangeEvaluatorFactory : public FieldEvaluatorFactory<T>
{
public:
  explicit OmegaHLagrangeEvaluatorFactory(
    std::shared_ptr<const OmegaHLagrangeLayout> layout)
    : layout_(std::move(layout)),
      search_(detail::MakeSearch(layout_->GetMesh()))
  {
  }

  const FieldLayout& GetLayout() const override { return *layout_; }

  CoordinateSystem GetCoordinateSystem() const override
  {
    return layout_->GetDOFHolderCoordinates().GetCoordinateSystem();
  }

  bool HasDOFHolderCoordinates() const override { return true; }

  bool SupportsNearestBoundary() const override { return false; }

  std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    const EvaluationRequest& request) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto coords = request.coords;
    const auto policy = request.policy;
    if (coords.GetCoordinateSystem() != GetCoordinateSystem()) {
      throw pcms_error(
        "OmegaHLagrangeEvaluatorFactory: coordinate system mismatch");
    }
    if (policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY) {
      throw pcms_error(
        "OmegaHLagrangeEvaluatorFactory: NearestBoundary is not supported");
    }

    auto raw_coords = coords.GetValues();
    LO n_pts = static_cast<LO>(raw_coords.extent(0));
    int mesh_dim = layout_->GetMesh().dim();
    Omega_h::Mesh& mesh = const_cast<Omega_h::Mesh&>(layout_->GetMesh());

    OmegaHLagrangeLocHint hint = std::visit(
      [&](auto& search) {
        using SearchT = std::decay_t<decltype(search)>;
        constexpr int Dim = SearchT::DIM;
        Kokkos::View<Real**, DeviceMemorySpace> coords_d(
          "coords_d", raw_coords.extent(0), raw_coords.extent(1));

        detail::CopyCoordsFunctor<Dim> copy_functor(coords_d, raw_coords);
        Kokkos::parallel_for("copy_coords", n_pts, copy_functor);

        auto results_d = search(coords_d);
        auto owning_ids = search.GetOwningElementIds(results_d);

        return detail::BuildLagrangeLocHint<Dim>(
          mesh, mesh_dim, results_d, coords_d, owning_ids, policy.mode);
      },
      search_);

    return std::make_unique<OmegaHLagrangePointEvaluator<T>>(
      layout_, std::move(hint), policy.fill_value);
  }

  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return layout_->GetDOFHolderCoordinates();
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
  mutable std::variant<GridPointSearch2D, GridPointSearch3D> search_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H
