#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>
#include <cmath>

#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/bounding_box.h"

namespace pcms
{
//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail
{
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform2DGrid[1]> grid,
                              int num_grid_cells);
} // namespace detail

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

/**
 * Compute the Euclidean distance from a point to the closest edge of a triangle
 * using only its barycentric coordinates with respect to that triangle.
 *
 * Given a triangle with vertex coordinates `coords` and a point represented by
 * barycentric coordinates `xi` (with respect to those vertices), the distance
 * from the point to the edge opposite vertex i is `xi[i] * h_i`, where `h_i`
 * is the altitude from vertex i to its opposite edge. The altitude can be
 * computed from the triangle geometry as `h_i = 2A / |e_i|`, where `A` is the
 * triangle area and `|e_i|` is the length of the edge opposite vertex i.
 *
 * This function generalizes to N-dimensional embeddings of a 2-simplex
 * (triangle) by computing the triangle area from two edge vectors using the
 * parallelogram formula valid in any dimension:
 *   area = 0.5 * sqrt(||a||^2 ||b||^2 - (a·b)^2)
 *
 * Template parameter N is the embedding dimension (e.g., 2 for planar, 3 for
 * spatial), while the simplex is always a triangle (3 vertices).
 */
template <int N>
[[nodiscard]] KOKKOS_FUNCTION inline Real
distance_to_closest_edge_from_barycentric(const Omega_h::Matrix<N, 3>& coords,
                                          const Omega_h::Vector<3>& xi)
{
  // Edge vectors originating from vertex 0
  Omega_h::Vector<N> a = coords[1] - coords[0];
  Omega_h::Vector<N> b = coords[2] - coords[0];

  // Parallelogram area magnitude squared = ||a||^2 ||b||^2 - (a·b)^2
  const Real aa = Omega_h::inner_product(a, a);
  const Real bb = Omega_h::inner_product(b, b);
  const Real ab = Omega_h::inner_product(a, b);
  Real parallelogram_area_sq = aa * bb - ab * ab;

  // Guard against degeneracy and tiny negative due to FP roundoff
  if (parallelogram_area_sq <= static_cast<Real>(0))
    return 0;
  const Real area = static_cast<Real>(0.5) * std::sqrt(parallelogram_area_sq);
  const Real two_area = static_cast<Real>(2) * area;

  // Opposite edge lengths as a vector [|e(1,2)|, |e(0,2)|, |e(0,1)|]
  Omega_h::Vector<3> ell;
  for (int i = 0; i < 3; ++i) {
    const int j = (i + 1) % 3;
    const int k = (i + 2) % 3;
    ell[i] = Omega_h::norm(coords[k] - coords[j]);
  }

  // Altitudes h_i = 2A / |e_i| with degeneracy guard, stored as a vector
  const Real eps = static_cast<Real>(1e-15);
  Omega_h::Vector<3> h;
  Omega_h::Vector<3> d;
  for (int i = 0; i < 3; ++i) {
    const Real hi = (ell[i] > eps) ? (two_area / ell[i]) : static_cast<Real>(0);
    h[i] = hi;
    d[i] = xi[i] * hi; // Distances to opposite edges using barycentric weights
  }

  // Return the minimum component (avoid std::min for device-compatibility)
  Real m = d[0];
  for (int i = 1; i < 3; ++i)
    m = (d[i] < m) ? d[i] : m;
  return m;
}

template <int dim>
class PointLocalizationSearch
{
public:
  struct Result
  {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2,
      REGION = 3
    };

    Dimensionality dimensionality;
    LO element_id; // ID of the actual entity (vertex/edge/face/region)
    Omega_h::Vector<dim + 1> parametric_coords;
  };

  static constexpr auto DIM = dim;

  using PointSearchTolerances = Kokkos::View<Real[DIM]>;

  explicit PointLocalizationSearch(const PointSearchTolerances& tolerances)
    : tolerances_(tolerances)
  {
    assert(tolerances_.is_allocated());
  }

  virtual Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [dim]> point) const = 0;
  [[nodiscard]] virtual LO GetOwningElementId(const Result& result) const = 0;
  [[nodiscard]] virtual Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const = 0;
  virtual ~PointLocalizationSearch() = default;

protected:
  PointSearchTolerances tolerances_;
};

using PointLocalizationSearch2D = PointLocalizationSearch<2>;
using PointLocalizationSearch3D = PointLocalizationSearch<3>;

class GridPointSearch2D : public PointLocalizationSearch2D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch2D::Result;

  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                    const PointSearchTolerances& tolerances);

  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) const override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const override;

private:
  [[nodiscard]] LO get_smallest_owner_face_id(
    Result::Dimensionality dimensionality, LO element_id) const;

  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj edges2faces_up_;
  Omega_h::Adj verts2faces_up_;
  Kokkos::View<Uniform2DGrid[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

class GridPointSearch3D : public PointLocalizationSearch3D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch3D::Result;

  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                    const PointSearchTolerances& tolerances);

  /**
   *  Given a point in global coordinates, returns the id of the tetrahedron (3D
   * element) that the point lies within and the parametric coordinate of the
   * point within the tetrahedron. If the point does not lie within any
   * tetrahedron element, then the id will be a negative number and (TODO) will
   * return a negative id of the closest element.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) const override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const override;

private:
  [[nodiscard]] LO get_smallest_owner_region_id(
    Result::Dimensionality dimensionality, LO element_id) const;

  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj verts2regions_up_;
  Omega_h::Adj edges2regions_up_;
  Omega_h::Adj faces2regions_up_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

} // namespace pcms
#endif // PCMS_COUPLING_POINT_SEARCH_H
