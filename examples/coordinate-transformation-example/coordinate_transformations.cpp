// Coordinate transformations in PCMS.
//
//   (a) map points and re-express vector components with a CoordinateMap
//   (b) transfer where the source points are cylindrical but both fields store
//       Cartesian components: only the points are mapped
//   (c) transfer where the source points and components are both cylindrical
//       and the target's are Cartesian: points are mapped, values rotated
//
// Each case prints the maximum error against the exact answer.
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/coordinate_systems/cartesian.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/transformed_source.hpp>
#include <pcms/transfer/interpolator.h>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using namespace pcms;
using DeviceView = Kokkos::View<Real**, DeviceMemorySpace>;
using HostView = Kokkos::View<Real**, HostMemorySpace>;

namespace
{

DeviceView ToDevice(const std::vector<Real>& flat, int width)
{
  const int n = static_cast<int>(flat.size()) / width;
  HostView host("host", n, width);
  for (int i = 0; i < n; ++i)
    for (int d = 0; d < width; ++d)
      host(i, d) = flat[i * width + d];
  DeviceView dev("dev", n, width);
  Kokkos::deep_copy(dev, host);
  return dev;
}

HostView ToHost(Rank2View<const Real, DeviceMemorySpace> v)
{
  DeviceView dev("dev", v.extent(0), v.extent(1));
  ConvertMismatchLayoutView2D(dev, v);
  HostView host("host", v.extent(0), v.extent(1));
  DeepCopyMismatchLayouts(host, dev);
  return host;
}

ValueBasis CylindricalVector()
{
  return {csys::CylindricalRThetaZ::Create(), ComponentScaling::Physical,
          values::Vector};
}

/// Sets a 3-component field from f(p0, p1, p2, out) at its DOF holders, where
/// p are the holder coordinates in the field's space.
template <typename F>
void FillField(Field<Real>& field, F f)
{
  auto p = ToHost(field.GetLayout().GetDOFHolderCoordinates().GetValues());
  const int n = static_cast<int>(p.extent(0));
  HostView data("data", n, 3);
  for (int i = 0; i < n; ++i)
    f(p(i, 0), p(i, 1), p(i, 2), &data(i, 0));
  field.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    field.GetData().GetValueBasis(), MakeConstRank2View(data)));
}

/// Max over the field's DOF holders of |b - (x, y, z)|.
Real MaxErrorFromXYZ(const Field<Real>& b)
{
  auto x = ToHost(b.GetLayout().GetDOFHolderCoordinates().GetValues());
  auto v = b.GetDOFHolderDataHost().GetValues();
  Real err = 0.0;
  for (size_t i = 0; i < x.extent(0); ++i)
    for (int d = 0; d < 3; ++d)
      err = std::fmax(err, std::fabs(v(i, d) - x(i, d)));
  return err;
}

/// Source: a box whose coordinates are (r, theta, z) with r in [0, 2] and
/// theta in [0, 1.6]. Target: the Cartesian unit cube, which lies inside it.
void BuildSpaces(Omega_h::Library& lib, Omega_h::Mesh& src_mesh,
                 Omega_h::Mesh& tgt_mesh,
                 std::shared_ptr<LagrangeFunctionSpace>& src,
                 std::shared_ptr<LagrangeFunctionSpace>& tgt)
{
  src_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1, 8,
                                8, 6, false);
  src = LagrangeFunctionSpace::FromMesh(
    src_mesh, 1, 3, csys::CylindricalRThetaZ::Create(), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
  tgt_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 4,
                                4, 4, false);
  tgt = LagrangeFunctionSpace::FromMesh(tgt_mesh, 1, 3,
                                        csys::Cartesian::Create(3), "global",
                                        LagrangeFunctionSpace::Backend::OmegaH);
}

void CaseA()
{
  std::puts("(a) coordinate transformation");
  auto xyz = ToDevice({1.0, 0.0, 0.5, 0.0, 2.0, -1.0, 1.0, 1.0, 0.0}, 3);
  auto bound = CartesianToCylindrical{}.Bind(CoordinateView<DeviceMemorySpace>(
    csys::Cartesian::Create(3), MakeConstRank2View(xyz)));
  auto rtz = ToHost(bound->MappedPoints().GetValues());

  // Values run opposite to points: this map re-expresses cylindrical
  // components as Cartesian ones, in place. Here, e_theta at each point.
  auto v = ToDevice({0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0}, 3);
  auto v_cart = bound->TransformValues(
    ValueView<Real, DeviceMemorySpace>(CylindricalVector(), MakeRank2View(v)));
  auto vh = ToHost(MakeConstRank2View(v));
  Real err = 0.0;
  for (int i = 0; i < 3; ++i) {
    err = std::fmax(err, std::fabs(vh(i, 0) + std::sin(rtz(i, 1))));
    err = std::fmax(err, std::fabs(vh(i, 1) - std::cos(rtz(i, 1))));
    std::printf("    (r,th,z)=(%.3f, %.3f, %.3f)  e_theta=(%.3f, %.3f, %.3f)\n",
                rtz(i, 0), rtz(i, 1), rtz(i, 2), vh(i, 0), vh(i, 1), vh(i, 2));
  }
  std::printf("    output basis '%s', max |e_theta - (-sin, cos, 0)| = %.2e\n",
              std::string(v_cart.GetBasis().system->Kind()).c_str(), err);
}

void CaseB(Omega_h::Library& lib)
{
  std::puts("(b) cylindrical points, Cartesian components on both sides");
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  std::shared_ptr<LagrangeFunctionSpace> src, tgt;
  BuildSpaces(lib, src_mesh, tgt_mesh, src, tgt);

  // A borrowed basis: the space's points are (r, theta, z), but the field's
  // components are Cartesian. B = (x, y, z) = (r cos theta, r sin theta, z).
  auto b_src =
    src->CreateFunction<Real>("b", values::Vector, csys::Cartesian::Create(3));
  FillField(b_src, [](Real r, Real theta, Real z, Real* out) {
    out[0] = r * std::cos(theta);
    out[1] = r * std::sin(theta);
    out[2] = z;
  });
  auto b_tgt = tgt->CreateFunction<Real>("b", values::Vector);

  // The values are already in the target's basis, so the map passes them
  // through and only the points are mapped.
  TransformedSource src_as_cartesian(
    src, std::make_shared<CartesianToCylindrical>());
  Interpolator<Real> op(src_as_cartesian, *tgt);
  op.Apply(b_src, b_tgt);
  // Not exact: B's Cartesian components are not linear in (r, theta, z), so
  // the P1 source field only approximates them.
  std::printf("    via Interpolator:  max |B - (x,y,z)| = %.2e\n",
              MaxErrorFromXYZ(b_tgt));
}

void CaseC(Omega_h::Library& lib)
{
  std::puts("(c) cylindrical points and components -> Cartesian");
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  std::shared_ptr<LagrangeFunctionSpace> src, tgt;
  BuildSpaces(lib, src_mesh, tgt_mesh, src, tgt);

  // B = r e_r + z e_z in physical cylindrical components, which is (x, y, z)
  // in Cartesian ones. Its cylindrical components are linear in
  // (r, theta, z), so the transfer is exact.
  auto b_src =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  FillField(b_src, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });
  auto b_tgt = tgt->CreateFunction<Real>("b", values::Vector);

  // `to_source` must end in the source's system: it carries the target's
  // Cartesian points to (r, theta, z), and the source's cylindrical
  // components back to Cartesian.
  TransformedSource src_as_cartesian(
    src, std::make_shared<CartesianToCylindrical>());
  Interpolator<Real> op(src_as_cartesian, *tgt);
  op.Apply(b_src, b_tgt);
  std::printf("    max |B - (x,y,z)| = %.2e\n", MaxErrorFromXYZ(b_tgt));
}

} // namespace

int main(int argc, char** argv)
{
  Omega_h::Library lib(&argc, &argv);
  CaseA();
  CaseB(lib);
  CaseC(lib);
  return 0;
}
