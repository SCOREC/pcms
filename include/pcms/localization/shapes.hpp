#ifndef PCMS_LOCALIZATION_SHAPES_HPP
#define PCMS_LOCALIZATION_SHAPES_HPP

#include <Kokkos_Core.hpp>

#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"

namespace pcms
{

/** A point in `Dim`-dimensional space. */
template <int Dim>
struct Point
{
  Kokkos::Array<Real, Dim> coords;
};

/** A line segment in `Dim`-dimensional space, from `start` to `end`. */
template <int Dim>
struct Segment
{
  Kokkos::Array<Real, Dim> start;
  Kokkos::Array<Real, Dim> end;
};

/**
 * A simplex in `Dim`-dimensional space: a triangle in 2D, a tetrahedron in 3D.
 */
template <int Dim>
struct Simplex
{
  // vertices[k][d]: coordinate d of vertex k
  Kokkos::Array<Kokkos::Array<Real, Dim>, Dim + 1> vertices;

  /** Vertex `k`, for `k` in `[0, Dim]`. */
  [[nodiscard]] KOKKOS_INLINE_FUNCTION Point<Dim> Vertex(int k) const
  {
    return {vertices[k]};
  }

  /** The edge from vertex `i` to vertex `j`. */
  [[nodiscard]] KOKKOS_INLINE_FUNCTION Segment<Dim> Edge(int i, int j) const
  {
    return {vertices[i], vertices[j]};
  }
};

/**
 * The spatial dimension of a shape type. A shape type other than the ones
 * here defines a static `dimension`.
 */
template <typename S>
inline constexpr int ShapeDimension = S::dimension;
template <int Dim>
inline constexpr int ShapeDimension<Point<Dim>> = Dim;
template <int Dim>
inline constexpr int ShapeDimension<Segment<Dim>> = Dim;
template <int Dim>
inline constexpr int ShapeDimension<Simplex<Dim>> = Dim;
template <int Dim>
inline constexpr int ShapeDimension<AABBox<Dim>> = Dim;

/**
 * The elements of a simplex mesh as simplices.
 * \param coords vertex coordinates, `Dim` per vertex
 * \param elem_verts the `Dim + 1` vertices of each element
 */
template <int Dim>
[[nodiscard]] Kokkos::View<const Simplex<Dim>*> MeshSimplices(
  Kokkos::View<const Real*> coords, Kokkos::View<const LO*> elem_verts)
{
  const auto nelems = static_cast<LO>(elem_verts.extent(0) / (Dim + 1));
  Kokkos::View<Simplex<Dim>*> elements("mesh_simplices", nelems);
  Kokkos::parallel_for(
    "mesh_simplices", nelems, KOKKOS_LAMBDA(LO e) {
      for (int k = 0; k < Dim + 1; ++k) {
        const LO v = elem_verts(e * (Dim + 1) + k);
        for (int d = 0; d < Dim; ++d) {
          elements(e).vertices[k][d] = coords(v * Dim + d);
        }
      }
    });
  return elements;
}

} // namespace pcms

#include "pcms/localization/detail/shapes.hpp"

#endif // PCMS_LOCALIZATION_SHAPES_HPP
