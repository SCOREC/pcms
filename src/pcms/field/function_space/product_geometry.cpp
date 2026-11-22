#include "pcms/geometry/product_geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/coordinate_systems/parameter.hpp"
#include "product_binding.hpp"
#include <vector>

namespace pcms
{

namespace
{

/// The columns [first, first + count) of `points` as a device view
/// (a free function: CUDA does not allow device lambdas in members).
Kokkos::View<Real**, DeviceMemorySpace> Columns(
  Kokkos::View<const Real**, DeviceMemorySpace> points, int first, int count)
{
  Kokkos::View<Real**, DeviceMemorySpace> out("product_factor_points",
                                              points.extent(0), count);
  pcms::detail::ForEach(
    "ProductColumns", static_cast<LO>(points.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int d = 0; d < count; ++d) {
        out(p, d) = points(p, first + d);
      }
    });
  return out;
}

CoordinateView<DeviceMemorySpace> FactorView(
  const Geometry& g, const Kokkos::View<Real**, DeviceMemorySpace>& points)
{
  return CoordinateView<DeviceMemorySpace>(g.GetCoordinateSystem(),
                                           MakeConstRank2View(points));
}

/// The parameter system of the product of `a` and `b`, with their
/// identifications.
std::shared_ptr<const CoordinateSystem> ProductSystem(const Geometry& a,
                                                      const Geometry& b)
{
  const int da = a.GetCoordinateSystem()->Dimension();
  std::vector<Identification> identifications;
  for (const auto& id : a.GetCoordinateSystem()->Identifications()) {
    identifications.push_back(id);
  }
  for (auto id : b.GetCoordinateSystem()->Identifications()) {
    if (id.map != nullptr) {
      throw pcms_error("CreateProductGeometry: the second factor has a "
                       "twist-shift identification; pass the product's "
                       "coordinate system");
    }
    id.axis += da;
    identifications.push_back(id);
  }
  return csys::CreateParameterSystem("product",
                                     da + b.GetCoordinateSystem()->Dimension(),
                                     std::move(identifications));
}

} // namespace

ProductBinding::ProductBinding(std::shared_ptr<const Geometry> geometry,
                               const BindRequest& request)
  : BoundPoints(std::move(geometry), request)
{
  Locate(false);
}

void ProductBinding::RebindImpl()
{
  Locate(true);
}

void ProductBinding::Locate(bool rebind)
{
  const auto& product = static_cast<const ProductGeometry&>(*GetGeometry());
  const auto points = LocatedPointValues();
  const int da = product.FactorA()->GetCoordinateSystem()->Dimension();
  const int db = product.FactorB()->GetCoordinateSystem()->Dimension();
  const auto xa = Columns(points, 0, da);
  const auto xb = Columns(points, da, db);
  const auto view_a = FactorView(*product.FactorA(), xa);
  const auto view_b = FactorView(*product.FactorB(), xb);
  if (rebind) {
    a_ = product.FactorA()->Rebind(std::move(a_), view_a);
    b_ = product.FactorB()->Rebind(std::move(b_), view_b);
  } else {
    a_ = product.FactorA()->Bind(view_a, Policy());
    b_ = product.FactorB()->Bind(view_b, Policy());
  }
  if (a_->MissingPoints().extent(0) + b_->MissingPoints().extent(0) > 0) {
    SetFound(pcms::detail::FoundMask(
      NumPoints(), {a_->MissingPoints(), b_->MissingPoints()}));
  }
}

ProductGeometry::ProductGeometry(std::shared_ptr<const Geometry> a,
                                 std::shared_ptr<const Geometry> b,
                                 std::shared_ptr<const CoordinateSystem> system)
  : a_(std::move(a)), b_(std::move(b)), system_(std::move(system))
{
  if (a_ == nullptr || b_ == nullptr) {
    throw pcms_error("CreateProductGeometry: the factors must not be null");
  }
  if (system_ == nullptr) {
    system_ = ProductSystem(*a_, *b_);
  }
  if (system_->HasMetric() ||
      system_->Dimension() != a_->GetCoordinateSystem()->Dimension() +
                                b_->GetCoordinateSystem()->Dimension()) {
    throw pcms_error("CreateProductGeometry: the product's coordinate system "
                     "must have no metric and the summed dimension of the "
                     "factors");
  }
}

std::shared_ptr<BoundPoints> ProductGeometry::BindImpl(
  const BindRequest& request) const
{
  return std::make_shared<ProductBinding>(shared_from_this(), request);
}
Quadrature ProductGeometry::TensorQuadrature(const Quadrature& qa,
                                             const Quadrature& qb) const
{
  const auto pa = qa.points->LocatedPoints().GetValues();
  const auto pb = qb.points->LocatedPoints().GetValues();
  const auto na = static_cast<LO>(pa.extent(0));
  const auto nb = static_cast<LO>(pb.extent(0));
  const int da = static_cast<int>(pa.extent(1));
  const int db = static_cast<int>(pb.extent(1));
  const auto wa = qa.weights;
  const auto wb = qb.weights;
  Kokkos::View<Real**, DeviceMemorySpace> points(
    "tensor_points", static_cast<std::size_t>(na) * nb, da + db);
  Kokkos::View<Real*, DeviceMemorySpace> weights(
    "tensor_weights", static_cast<std::size_t>(na) * nb);
  pcms::detail::ForEach(
    "TensorQuadrature", na * nb, KOKKOS_LAMBDA(const LO i) {
      const LO ia = i / nb;
      const LO ib = i % nb;
      for (int d = 0; d < da; ++d) {
        points(i, d) = pa(ia, d);
      }
      for (int d = 0; d < db; ++d) {
        points(i, da + d) = pb(ib, d);
      }
      weights(i) = wa(ia) * wb(ib);
    });
  return Quadrature{Bind(CoordinateView<DeviceMemorySpace>(
                      system_, MakeConstRank2View(points))),
                    weights,
                    {}};
}

Quadrature ProductGeometry::CreateQuadrature(
  const detail::QuadratureRule& rule) const
{
  return TensorQuadrature(detail::CreateQuadrature(*a_, rule),
                          detail::CreateQuadrature(*b_, rule));
}

std::shared_ptr<const Geometry> CreateProductGeometry(
  std::shared_ptr<const Geometry> a, std::shared_ptr<const Geometry> b,
  std::shared_ptr<const CoordinateSystem> system)
{
  return std::make_shared<ProductGeometry>(std::move(a), std::move(b),
                                           std::move(system));
}

} // namespace pcms
