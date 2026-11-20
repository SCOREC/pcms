// The product of two function spaces: the basis phi_i(a) psi_k(b), DOF (i, k).
// Its evaluation matrix is the Kronecker product of the factors' rows, built
// when points are bound; the product evaluators are then the sparse-row
// evaluators of the 1D spaces.

#include "pcms/field/function_spaces.hpp"
#include "pcms/field/holder_of_gid.hpp"
#include "pcms/geometry/product_geometry.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "product_binding.hpp"
#include "sparse_rows.hpp"
#include <set>

namespace pcms
{

namespace
{

Kokkos::View<Real**, HostMemorySpace> HostCoordinates(const FieldLayout& layout)
{
  const auto values = layout.GetDOFHolderCoordinates().GetValues();
  Kokkos::View<Real**, DeviceMemorySpace> device(
    "holder_coordinates", values.extent(0), values.extent(1));
  pcms::detail::CopyDeviceRank2ViewToDeviceView(device, values);
  const auto mirror =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, device);
  Kokkos::View<Real**, HostMemorySpace> host(
    "holder_coordinates", values.extent(0), values.extent(1));
  for (std::size_t i = 0; i < host.extent(0); ++i) {
    for (std::size_t j = 0; j < host.extent(1); ++j) {
      host(i, j) = mirror(i, j);
    }
  }
  return host;
}

/// The layout of the product of two serial layouts: holder (i, k) is
/// i * (holders of b) + k of the factors' local holders, its global id is
/// gid_a * (global holders of b) + gid_b, it is owned when both are, and its
/// position is the factors' positions side by side.
class ProductLayout final : public FieldLayout
{
public:
  ProductLayout(std::shared_ptr<const Geometry> geometry, const FieldLayout& a,
                const FieldLayout& b)
    : num_global_(a.GetNumGlobalDofHolder() * b.GetNumGlobalDofHolder()),
      dimension_(a.GetDimension() + b.GetDimension())
  {
    if (a.HasGhosts() || b.HasGhosts()) {
      throw pcms_error("CreateProductSpace: the factors' layouts must have "
                       "no ghosts (the product is serial)");
    }
    const LO na = a.GetNumDofHolder();
    const LO nb = b.GetNumDofHolder();
    size_ = na * nb;
    owned_ = Kokkos::View<bool*, HostMemorySpace>("owned", size_);
    gids_ = Kokkos::View<GO*, HostMemorySpace>("gids", size_);
    zeros_ = Kokkos::View<LO*, HostMemorySpace>("zeros", size_);
    coords_ =
      Kokkos::View<Real**, HostMemorySpace>("coords", size_, dimension_);
    const auto gid_a = a.GetGidsHost();
    const auto gid_b = b.GetGidsHost();
    const auto owned_a = a.GetOwnedHost();
    const auto owned_b = b.GetOwnedHost();
    const auto xa = HostCoordinates(a);
    const auto xb = HostCoordinates(b);
    for (LO i = 0; i < na; ++i) {
      for (LO k = 0; k < nb; ++k) {
        const LO holder = i * nb + k;
        owned_(holder) = owned_a(i) && owned_b(k);
        gids_(holder) = gid_a(i) * b.GetNumGlobalDofHolder() + gid_b(k);
        for (int d = 0; d < a.GetDimension(); ++d) {
          coords_(holder, d) = xa(i, d);
        }
        for (int d = 0; d < b.GetDimension(); ++d) {
          coords_(holder, a.GetDimension() + d) = xb(k, d);
        }
      }
    }
    SetCoordinateSystem(geometry->GetCoordinateSystem());
    SetGeometry(std::move(geometry));
  }
  int GetNumComponents() const override { return 1; }
  LO GetNumDofHolder() const override { return size_; }
  GO GetNumGlobalDofHolder() const override { return num_global_; }
  bool IsDistributed() const override { return false; }
  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override
  {
    return MakeRank1View(owned_);
  }
  GlobalIDView<HostMemorySpace> GetGidsHost() const override
  {
    return GlobalIDView<HostMemorySpace>(gids_.data(), size_);
  }
  EntOffsetsArray GetEntOffsets() const override
  {
    const auto n = static_cast<std::size_t>(size_);
    return {0, n, n, n, n};
  }
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    if (device_coords_.extent(0) != static_cast<std::size_t>(size_)) {
      device_coords_ = pcms::detail::ToDeviceView(coords_);
    }
    return CoordinateView<DeviceMemorySpace>(
      GetCoordinateSystem(), MakeConstRank2View(device_coords_));
  }
  int GetDimension() const override { return dimension_; }
  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override
  {
    return MakeRank1View(zeros_);
  }
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override
  {
    return MakeRank1View(zeros_);
  }

private:
  GO num_global_;
  int dimension_;
  LO size_ = 0;
  Kokkos::View<bool*, HostMemorySpace> owned_;
  Kokkos::View<GO*, HostMemorySpace> gids_;
  Kokkos::View<LO*, HostMemorySpace> zeros_;
  Kokkos::View<Real**, HostMemorySpace> coords_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> device_coords_;
};

/// The factors' evaluation matrices at a set of bound points: of the values
/// (v) and, for Gradient, of the partials (g, row p * dimension + j).
struct FactorMatrices
{
  PointMatrix<Real> va, vb, ga, gb;
};

/// The rows of the two factors that output row r combines.
struct RowSelect
{
  int num_out, da, db;
  bool gradient;
  KOKKOS_INLINE_FUNCTION void operator()(const LO r, bool& a_gradient,
                                         bool& b_gradient, LO& row_a,
                                         LO& row_b) const
  {
    const LO p = r / num_out;
    const int c = r % num_out;
    a_gradient = gradient && c < da;
    b_gradient = gradient && c >= da;
    row_a = a_gradient ? p * da + c : p;
    row_b = b_gradient ? p * db + (c - da) : p;
  }
};

/// The Kronecker rows: row (p, c) of the output is the product of a row of a's
/// matrix and a row of b's: the value rows for Value; for Gradient component c
/// of a, its partial row with b's value row, and for b the other way round.
/// (A free function: CUDA does not allow device lambdas in members.)
detail::SparseRows KroneckerRows(
  const FactorMatrices& m, Kokkos::View<const bool*, DeviceMemorySpace> found,
  int da, int db, bool gradient, const detail::HolderOfGid& holder_a,
  const detail::HolderOfGid& holder_b, LO holders_b, GO global_b)
{
  const auto num_points = static_cast<LO>(found.extent(0));
  const int num_out = gradient ? da + db : 1;
  const LO num_rows = num_points * num_out;
  const RowSelect select{num_out, da, db, gradient};
  detail::SparseRows rows;
  rows.num_out = num_out;
  rows.offsets =
    Kokkos::View<LO*, DeviceMemorySpace>("product_offsets", num_rows + 1);
  const auto offsets = rows.offsets;
  Kokkos::parallel_scan(
    "ProductRowOffsets",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_rows + 1),
    KOKKOS_LAMBDA(const LO r, LO& partial, const bool final) {
      if (final) {
        offsets(r) = partial;
      }
      if (r < num_rows) {
        bool ag, bg;
        LO ra, rb;
        select(r, ag, bg, ra, rb);
        const auto& a = ag ? m.ga : m.va;
        const auto& b = bg ? m.gb : m.vb;
        if (found(r / num_out)) {
          partial += (a.row_offsets(ra + 1) - a.row_offsets(ra)) *
                     (b.row_offsets(rb + 1) - b.row_offsets(rb));
        }
      }
    });
  LO nnz = 0;
  Kokkos::deep_copy(nnz, Kokkos::subview(offsets, num_rows));
  rows.columns = Kokkos::View<GO*, DeviceMemorySpace>("product_columns", nnz);
  rows.holders = Kokkos::View<LO*, DeviceMemorySpace>("product_holders", nnz);
  rows.values = Kokkos::View<Real*, DeviceMemorySpace>("product_values", nnz);
  const auto columns = rows.columns;
  const auto holders = rows.holders;
  const auto values = rows.values;
  Kokkos::parallel_for(
    "ProductRows", Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_rows),
    KOKKOS_LAMBDA(const LO r) {
      if (!found(r / num_out)) {
        return;
      }
      bool ag, bg;
      LO ra, rb;
      select(r, ag, bg, ra, rb);
      const auto& a = ag ? m.ga : m.va;
      const auto& b = bg ? m.gb : m.vb;
      LO k = offsets(r);
      for (LO i = a.row_offsets(ra); i < a.row_offsets(ra + 1); ++i) {
        for (LO j = b.row_offsets(rb); j < b.row_offsets(rb + 1); ++j, ++k) {
          columns(k) = a.columns(i) * global_b + b.columns(j);
          holders(k) =
            holder_a(a.columns(i)) * holders_b + holder_b(b.columns(j));
          values(k) = a.values(i) * b.values(j);
        }
      }
    });
  return rows;
}

/// 1 for the points the bound points evaluate, 0 for the missing ones
Kokkos::View<const bool*, DeviceMemorySpace> FoundMask(const BoundPoints& bound)
{
  Kokkos::View<bool*, DeviceMemorySpace> found("product_found",
                                               bound.NumPoints());
  Kokkos::deep_copy(found, true);
  const auto missing = bound.MissingPoints();
  pcms::detail::ForEach(
    "ProductFoundMask", static_cast<LO>(missing.extent(0)),
    KOKKOS_LAMBDA(const LO i) { found(missing(i)) = false; });
  return found;
}

PointMatrix<Real> MatrixOf(const FunctionSpace& space,
                           std::shared_ptr<const BoundPoints> bound,
                           DifferentialOperator op)
{
  const auto evaluator = space.CreatePointEvaluator<Real>(std::move(bound), op);
  if (!evaluator->SupportsMatrix()) {
    throw pcms_error("CreateProductSpace: a factor's evaluator provides no "
                     "matrix (the Lagrange space of Omega_h needs the OmegaH "
                     "backend)");
  }
  return evaluator->GetMatrix();
}

} // namespace

std::shared_ptr<const FunctionSpace> CreateProductSpace(
  std::shared_ptr<const FunctionSpace> a,
  std::shared_ptr<const FunctionSpace> b,
  std::shared_ptr<const Geometry> geometry)
{
  if (a == nullptr || b == nullptr) {
    throw pcms_error("CreateProductSpace: the factors must not be null");
  }
  if (geometry == nullptr) {
    geometry = CreateProductGeometry(a->GetGeometry(), b->GetGeometry());
  }
  const auto product =
    std::dynamic_pointer_cast<const ProductGeometry>(geometry);
  if (product == nullptr || product->FactorA() != a->GetGeometry() ||
      product->FactorB() != b->GetGeometry()) {
    throw pcms_error("CreateProductSpace: the geometry must be the product "
                     "geometry of the factors' geometries");
  }
  if (a->GetLayout()->GetNumComponents() != 1 ||
      b->GetLayout()->GetNumComponents() != 1) {
    throw pcms_error("CreateProductSpace: the factors must be scalar");
  }
  auto layout =
    std::make_shared<ProductLayout>(geometry, *a->GetLayout(), *b->GetLayout());
  const auto holder_a = pcms::detail::MakeHolderOfGid(*a->GetLayout());
  const auto holder_b = pcms::detail::MakeHolderOfGid(*b->GetLayout());
  const LO holders_b = b->GetLayout()->GetNumDofHolder();
  const GO global_b = b->GetLayout()->GetNumGlobalDofHolder();
  const int da = a->GetGeometry()->GetCoordinateSystem()->Dimension();
  const int db = b->GetGeometry()->GetCoordinateSystem()->Dimension();
  std::set<DifferentialOperator> supported = {DifferentialOperator::Value};
  if (a->Supports(DifferentialOperator::Gradient) &&
      b->Supports(DifferentialOperator::Gradient)) {
    supported.insert(DifferentialOperator::Gradient);
  }
  auto evaluator = std::make_shared<detail::SparseRowEvaluatorFactory>(
    geometry, layout, std::move(supported),
    [a, b, da, db, holder_a, holder_b, holders_b,
     global_b](const std::shared_ptr<const BoundPoints>& bound,
               DifferentialOperator op) {
      const auto& binding = static_cast<const ProductBinding&>(*bound);
      const bool gradient = op == DifferentialOperator::Gradient;
      FactorMatrices m;
      m.va = MatrixOf(*a, binding.FactorA(), DifferentialOperator::Value);
      m.vb = MatrixOf(*b, binding.FactorB(), DifferentialOperator::Value);
      if (gradient) {
        m.ga = MatrixOf(*a, binding.FactorA(), DifferentialOperator::Gradient);
        m.gb = MatrixOf(*b, binding.FactorB(), DifferentialOperator::Gradient);
      }
      return KroneckerRows(m, FoundMask(*bound), da, db, gradient, holder_a,
                           holder_b, holders_b, global_b);
    });
  return std::make_shared<detail::QuadratureSpace>(
    layout, std::move(evaluator), geometry, [a, b, product] {
      return product->TensorQuadrature(a->CreateQuadrature(),
                                       b->CreateQuadrature());
    });
}

} // namespace pcms
