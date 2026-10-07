#include "arbor_x_localization.hpp"

// Access traits for Omega_h
// constructs bounding boxes for elements
template <int dim>
struct ArborX::AccessTraits<pcms::detail::Omega_h_Mesh_Adapt<dim>>
{
  using memory_space = typename Omega_h::ExecSpace::memory_space;

  static KOKKOS_FUNCTION int size(
    const pcms::detail::Omega_h_Mesh_Adapt<dim>& mesh)
  {
    return mesh.size_;
  }

  static KOKKOS_FUNCTION auto get(
    const pcms::detail::Omega_h_Mesh_Adapt<dim>& mesh, int i)
  {
    ArborX::Point<dim, Omega_h::Real> min = {INFINITY};
    ArborX::Point<dim, Omega_h::Real> max = {-INFINITY};
    for (int j = 0; j < dim + 1; ++j) {
      auto cell_vert_id = mesh.adjacency[(dim + 1) * i + j];
      for (int k = 0; k < dim; ++k) {
        Omega_h::Real curr_coord = mesh.coordinates[cell_vert_id * dim + k];
        if (min[k] > curr_coord)
          min[k] = curr_coord;
        if (max[k] < curr_coord)
          max[k] = curr_coord;
      }
    }
    return ArborX::Box(min, max);
  }
};

template <typename MemorySpace, int dim>
struct ArborX::AccessTraits<
  pcms::detail::Coordinate_View_Adapt<MemorySpace, dim>>
{
  using memory_space = MemorySpace;
  static KOKKOS_FUNCTION int size(
    pcms::detail::Coordinate_View_Adapt<MemorySpace, dim> const& coords)
  {
    return coords.points.extent(0);
  }
  static KOKKOS_FUNCTION auto get(
    pcms::detail::Coordinate_View_Adapt<MemorySpace, dim> const& coords, int i)
  {
    ArborX::Point<dim, double> ax_point;
    for (int j = 0; j < dim; j++)
      ax_point[j] = coords.points(i, j);
    return PredicateWithAttachment(intersects(ax_point), i);
  }
};

namespace pcms
{
namespace detail
{
/**
 * @brief Implementation of the Kronecker delta:
 *		 https://en.wikipedia.org/wiki/Kronecker_delta
 * @param i an integer
 * @param j another integer
 * @returns 1 if i==j, 0 if i i!= j
 */
KOKKOS_INLINE_FUNCTION
LO kronecker(LO i, LO j)
{
  return (LO)(i == j);
}

/**
 * @brief Omega_h::Mesh tagged with the templated dimension required for
 *		 compatability with ArborX::BVH
 */
template <int dim>
struct Omega_h_Mesh_Adapt
{
  const Omega_h::LO size_;
  const Omega_h::LOs adjacency;
  const Omega_h::Reals coordinates;
};

/**
 * @brief pcms CoordinateView tagged with the templated dimension and
 *		 memory space required for compatability with ArborX::BVH
 */
template <typename MemorySpace, int dim>
struct Coordinate_View_Adapt
{
  const Rank2View<const Real, MemorySpace,
                  default_layout_for_memory_space_t<MemorySpace>>
    points;
};

/**
 * @brief Mapping<2> represents a mapping from the global coordinate system
 *		 to the barycentric coordinate system of an element in a 2D
 *Omega_h::Mesh The class Mapping<2> calculates the mapping of any point in
 *global coordinate space to the barycentric coordinate space belonging to a
 *triangle (face) in an Omega_h::Mesh. An instance of Mapping<2> can, from the
 *barycentric coordinates constructed by it, determine whether that point
 *intersects the face, an edge, a vertex, or lies outside the face. It also
 *determines which edge or vertex the point intersects by calculating the offset
 *in the Omega_h::Adj::ab2b
 */
template <>
class Mapping<2>
{
public:
  // the dimension of the space
  static constexpr int DIM = 2;
  // CONSTRUCTORS
  /**
   * @brief Default constructor
   */
  KOKKOS_FUNCTION
  Mapping() = default;
  /**
   * @brief Constructs a barycentric mapping of a triangle in an Omega_h::Mesh
   * @param elem_index the index of the desired triangle or tetrahedron in the
   * Omega_h::Mesh
   * @param mesh the Omega_h::Mesh with the spatial information of the triangle
   */
  Mapping(Omega_h::Matrix<DIM, DIM + 1> const& triangle_,
          Kokkos::View<Real*> const& tolerances)
    : triangle(triangle_)
  {
    auto tol_h = Kokkos::create_mirror_view(tolerances);
    Kokkos::deep_copy(tol_h, tolerances);
    tolerances_[0] = tol_h(0);
    tolerances_[1] = tol_h(1);
    // printf("%d (%lf, %lf)\n", tol_h.size(), tolerances_[0], tolerances_[1]);
    bary_transform = {triangle[0] - triangle[2], triangle[1] - triangle[2]};
    triangle_area = 0.5 * fabs(Omega_h::determinant(bary_transform));
    bary_transform = Omega_h::invert(bary_transform);
  }
  /**
   * @brief Default destructor
   */
  KOKKOS_FUNCTION
  ~Mapping() = default;
  /**
   * @brief Computes the barycentric coordinates of a point in global space
   *
   * This function computes the barycentric coordinate of a point with respect
   * to a triangle in an Omega_h::Mesh
   *
   * @param p the point to compute the barycentric coordinates of
   * @returns an Omega_h::Vector<3> containing the barycentric coordinates
   * 		   of p
   */
  KOKKOS_FUNCTION
  Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const
  {
    Omega_h::Vector<2> coeffs = bary_transform * (p - triangle[2]);
    return {coeffs[0], coeffs[1], 1 - coeffs[0] - coeffs[1]};
  }
  /**
   * @brief Determines which entity of a dimension `ent_dim` a point intersects
   *		 from the barycentric coordinates
   *
   * If the point corresponding to the input barycentric coordinates intersects
   * any entities with dimension `ent_dim` bordering the element this mapping is
   * constructed from, this function determines the offset in the
   * Omega_h::Adj::ab2b structure coresponding to FACE -> `ent_dim` adjacency.
   * If the point does not intersect a border entity, this function returns -1
   *
   * @param ent_dim the dimension of the entities we want to check intersection
   * 				 with
   * @param bary_coords the barycentric coordinates computed by this mapping of
   *					 a point in global space
   * @returns -1 if the point does not intersect any entities, or the offset of
   *		   the intersected entity in the Omega_h::Adj::ab2b structure
   *		   coresponding to FACE -> `ent_dim` adjacency
   */
  KOKKOS_FUNCTION
  int which(int ent_dim, Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    if (dim == Omega_h::VERT) {
      return which_vert(bary_coords);
    }
    if (dim == Omega_h::EDGE) {
      return which_edge(bary_coords);
    }
    if (dim == Omega_h::FACE) {
      return within_elem(bary_coords);
    }
    return -1;
  }

private:
  /**
   * @brief Computes the vertex offset of the point corresponding
   * 		  to the input barycentric coordinates
   * @param bary_coords the input barycentric coordinates
   * @returns the vertex offset if the point corresponding to the given bary.
   * coordinates is within a certain (global) tolerance of a vertex -1 if the
   * point does not lie within the tolerance of any vertex
   */
  KOKKOS_FUNCTION
  int which_vert(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    for (int i = 0; i < 3; i++) {
      // distance of the point to vertex i
      Omega_h::Vector<2> error =
        (bary_coords[0] - kronecker(i, 0)) * triangle[0] +
        (bary_coords[1] - kronecker(i, 1)) * triangle[1] +
        (bary_coords[2] - kronecker(i, 2)) * triangle[2];
      if (Omega_h::norm_squared(error) <= tolerances_[0] * tolerances_[0])
        return i;
    }
    return -1;
  }
  /**
   * @brief Computes the edge offset of the point corresponding
   * 		 to the input barycentric coordinates
   * @param bary_coords the input barycentric coordinates
   * @returns the edge offset if the distnce from the point corresponding to the
   * 			given bary. coordinates to the edge with the respective
   * offset is within a certain (global) tolerance -1 if the point does not lie
   * within the tolerance of any vertex
   */
  KOKKOS_FUNCTION
  int which_edge(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    for (int i = 0; i <= 2; i++) {
      // Distance of the point to the edge opposite the ith vertex,
      // see docs for derivation.
      double dist = (2 * bary_coords[i] * triangle_area) *
                    (2 * bary_coords[i] * triangle_area);
      dist /= opposite_edge_len_sq(i);
      if (dist <= tolerances_[1] * tolerances_[1] &&
          bary_coords[(i + 1) % 3] >= 0 && bary_coords[(i + 2) % 3] >= 0)
        return (i + 1) % 3;
    }
    return -1;
  }
  /**
   * @brief Computes whether a point with the given barycentric coordinates
   * 		  is within a triangle
   * @param bary_coords the input barycentric coordinates
   * @returns 0 if the point is within the highest order element and -1
   * otherwise
   */
  KOKKOS_FUNCTION
  int within_elem(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    return -1 * (int)!(bary_coords[0] > 0 && bary_coords[1] > 0 &&
                       bary_coords[2] > 0);
  }
  /**
   * @brief Calculates the length of the edge opposite vertex i
   * @param i the vertex index opposite the desired edge
   * @returns the length of the ith edge
   */
  KOKKOS_FUNCTION
  double opposite_edge_len_sq(int i) const
  {
    return Omega_h::norm_squared(triangle[(i + 2) % 3] - triangle[(i + 1) % 3]);
  }

  // REPRESENTATION
  Omega_h::Vector<DIM> tolerances_;
  // representation of the barycentric coordinate mapping, source:
  // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Edge_approach
  Omega_h::Matrix<DIM, DIM> bary_transform; // Column-major order
  // Points defining the triangle
  Omega_h::Matrix<DIM, DIM + 1> triangle;
  // area of the triangle
  double triangle_area;
};

template <>
class Mapping<3>
{
public:
  static constexpr int DIM = 3;
  /**
   * @brief Default constructor
   */
  KOKKOS_FUNCTION
  Mapping() = default;

  Mapping(Omega_h::Matrix<DIM, DIM + 1> const& tetrahedron_,
          const Kokkos::View<Real*>& tolerances)
    : tetrahedron(tetrahedron_)
  {
    auto tol_h = Kokkos::create_mirror_view(tolerances);
    Kokkos::deep_copy(tol_h, tolerances);
    tolerances_[0] = tol_h(0);
    tolerances_[1] = tol_h(1);
    tolerances_[2] = tol_h(2);
    set_triangle_areas();
    bary_transform = {tetrahedron[0] - tetrahedron[3],
                      tetrahedron[1] - tetrahedron[3],
                      tetrahedron[2] - tetrahedron[3]};
    tetrahedron_volume = fabs(Omega_h::determinant(bary_transform)) / 6.;
    bary_transform = Omega_h::invert(bary_transform);
  }
  /**
   * @brief Default destructor
   */
  KOKKOS_FUNCTION
  ~Mapping() = default;
  /**
   * @brief Computes the barycentric coordinates of a point in global space
   * @param p the point to compute the barycentric coordinates of
   * @returns an Omega_h::Vector<dim + 1> containing the barycentric coordinates
   * 			of p
   */
  KOKKOS_FUNCTION
  Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const
  {
    Omega_h::Vector<3> coeffs = bary_transform * (p - tetrahedron[3]);
    return {coeffs[0], coeffs[1], coeffs[2],
            1 - coeffs[0] - coeffs[1] - coeffs[2]};
  }
  /**
   * @brief Determines which entity of dimension `ent_dim` a point intersects
   *		 from the barycentric coordinates
   *
   * If the point corresponding to the input barycentric coordinates intersects
   * any entities with dimension `ent_dim` bordering the element this mapping is
   * constructed from, this function determines the offset in the
   * Omega_h::Adj::ab2b structure coresponding to REGION -> `ent_dim` adjacency.
   * If the point does not intersect a border entity, this function returns -1
   *
   * @param ent_dim the dimension of the entities we want to check intersection
   * 				 with
   * @param bary_coords the barycentric coordinates computed by this mapping of
   *					 a point in global space
   * @returns -1 if the point does not intersect any entities, or the offset of
   *		   the intersected entity in the Omega_h::Adj::ab2b structure
   *		   coresponding to REGION -> `ent_dim` adjacency
   */
  KOKKOS_FUNCTION
  int which(int ent_dim, Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    if (dim == Omega_h::VERT) {
      return which_vert(bary_coords);
    }
    if (dim == Omega_h::EDGE) {
      return which_edge(bary_coords);
    }
    if (dim == Omega_h::FACE) {
      return which_face(bary_coords);
    }
    if (dim == Omega_h::REGION) {
      return within_elem(bary_coords);
    }
    return -1;
  }

private:
  /**
   * @brief Computes the vertex offset of the point corresponding
   *		 to the input barycentric coordinates
   * @param bary_coords the input barycentric coordinates
   * @returns the vertex offset if the point corresponding to the given bary.
   *coordinates is within a certain (global) tolerance of a vertex -1 if the
   *point does not lie within the tolerance of any vertex
   */
  KOKKOS_FUNCTION
  int which_vert(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    for (int i = 0; i < 4; i++) {
      // distance to the ith vertex
      Omega_h::Vector<3> error =
        (bary_coords[0] - kronecker(i, 0)) * tetrahedron[0] +
        (bary_coords[1] - kronecker(i, 1)) * tetrahedron[1] +
        (bary_coords[2] - kronecker(i, 2)) * tetrahedron[2] +
        (bary_coords[3] - kronecker(i, 3)) * tetrahedron[3];
      if (Omega_h::norm_squared(error) <= tolerances_[0] * tolerances_[0])
        return i;
    }
    return -1;
  }
  /**
   * @brief Computes the edge offset of the point corresponding
   * 		 to the input barycentric coordinates
   * @param bary_coords the input barycentric coordinates
   * @returns the edge offset if the distnce from the point corresponding to the
   * 		   given bary. coordinates to the edge with the respective
   * offset is within a certain (global) tolerance -1 if the point does not lie
   * within the tolerance of any vertex
   */
  KOKKOS_FUNCTION
  int which_edge(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    int edge_ = 0;
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++) {
        Omega_h::Vector<3> side1 = {0, 0, 0}, side2 = {0, 0, 0},
                           side3 = tetrahedron[i] - tetrahedron[j];
        for (int k = 0; k < 4; k++) {
          side1 += (bary_coords[k] - kronecker(i, k)) * tetrahedron[k];
          side2 += (bary_coords[k] - kronecker(j, k)) * tetrahedron[k];
        }

        // the norm of the cross product of two vectors is twice the area of the
        // triangle those vectors form
        double distnce_sq = Omega_h::norm_squared(Omega_h::cross(side1, side2));
        distnce_sq /= Omega_h::norm_squared(side3);

        if (distnce_sq <= tolerances_[1] * tolerances_[1] &&
            bary_coords[i] > 0 && bary_coords[j] > 0)
          return edge(edge_);
        edge_++;
      }
    }
    return -1;
  }
  /**
   * @brief Computes the face offset of the point corresponding
   * 		 to the input barycentric coordinates
   * @param bary_coords the input barycentric coordinates
   * @returns if the point corresponding to the bary. coords is within a
   * 		   global tolerance of a face, returns the offset of that face
   * 		   -1 if the point does not lie within the tolerance of any face
   * Algorithm source:
   * C.E. Passerello,
   * Interference detection using barycentric coordinates,
   * Mechanics Research Communications,
   * Volume 9, Issue 6, 1982, Pages 373-378,
   * https://doi.org/10.1016/0093-6413(82)90034-9.
   */
  KOKKOS_FUNCTION
  int which_face(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    for (int i = 0; i < 4; i++) {
      if (fabs(3 * tetrahedron_volume * bary_coords[i] / face_areas[i]) <=
            tolerances_[2] &&
          bary_coords[(i + 1) % 4] >= 0 && bary_coords[(i + 2) % 4] >= 0 &&
          bary_coords[(i + 3) % 4] >= 0) {
        return face(i);
      }
    }
    return -1;
  }
  /**
   * @brief Computes whether a point with the given barycentric coordinates
   * 		  is within a triangle
   * @param bary_coords the input barycentric coordinates
   * @returns 0 if the point is within the highet order element and -1 otherwise
   */
  KOKKOS_FUNCTION
  int within_elem(Omega_h::Vector<DIM + 1> const& bary_coords) const
  {
    return -1 * (int)!(bary_coords[0] >= 0 && bary_coords[1] >= 0 &&
                       bary_coords[2] >= 0 && bary_coords[3] >= 0);
  }
  /**
   * @brief Calculates the areas of each face in the tetrahedron the mapping
   * 		 corresponds to
   * @param index the element ID of the tetrahedron
   * @param mesh the Omega_h::Mesh with the spatial information for the
   * 			  tetrahedron corresponding to `index`
   */
  void set_triangle_areas()
  {
    for (int i = 0; i < 4; i++) {
      Omega_h::Vector<3> edge0 =
                           tetrahedron[(i + 1) % 4] - tetrahedron[(i + 3) % 4],
                         edge1 =
                           tetrahedron[(i + 2) % 4] - tetrahedron[(i + 3) % 4];
      Omega_h::Vector<3> cross = Omega_h::cross(edge0, edge1);
      face_areas[i] = 0.5 * Omega_h::norm(cross);
    }
  }
  // returns the actual face/edge offset from the barycentric offset
  KOKKOS_INLINE_FUNCTION int face(int i) const
  {
    return (OFFSETS & 3 << (i * 2)) >> (i * 2);
  }
  KOKKOS_INLINE_FUNCTION int edge(int i) const
  {
    return (i < 4 && i > 0) ? (OFFSETS & 3 << ((i - 1) * 2)) >> ((i - 1) * 2)
                            : i;
  }
  // representation
  Omega_h::Vector<DIM> tolerances_;
  Omega_h::Matrix<DIM, DIM> bary_transform; // Column-major order
  Omega_h::Matrix<DIM, DIM + 1> tetrahedron;
  Omega_h::Vector<DIM + 1> face_areas;
  double tetrahedron_volume;
  // This is used to calculate the face and edge offsets in the cannonical
  // ordering:
  // https://user-images.githubusercontent.com/56453280/74203616-39d27a80-4c3e-11ea-885d-b0260490e184.png
  static const char OFFSETS = 1 << 4 | 3 << 2 | 2;
};

/**
 * Functor to be called when a point intersects a triangle or tetrahedron
 */
class CallOnIntersect3D
{
public:
  using MemorySpace = TreePointSearch::MemorySpace;
  static constexpr int DIM = 3;

  CallOnIntersect3D(const Kokkos::View<Mapping<3>*, MemorySpace>& mappings_,
                    const Omega_h::LOs& adjacencies0,
                    const Omega_h::LOs& adjacencies1,
                    const Omega_h::LOs& adjacencies2,
                    Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace>&
                      dimensionalities_,
                    Kokkos::View<LO*, MemorySpace>& element_ids_,
                    Kokkos::View<Real**, MemorySpace>& parametric_coords_)
    : mappings(MakeRank1View(mappings_)),
      adjacencies{adjacencies0, adjacencies1, adjacencies2},
      dimensionalities(MakeRank1View(dimensionalities_)),
      element_ids(MakeRank1View(element_ids_)),
      parametric_coords(MakeRank2View(parametric_coords_)) {};

  /**
   * Intersection callback,
   */
  template <typename Predicate, typename Value>
  KOKKOS_FUNCTION void operator()(Predicate const& predicate,
                                  Value const& val) const
  {
    ArborX::Point<DIM, Omega_h::Real> const& ax =
      ArborX::getGeometry(predicate);
    Omega_h::Vector<DIM> point{ax[0], ax[1]};
    int point_ind = ArborX::getData(predicate);

    detail::Mapping<2> const& tm = mappings(val.index);

    // calculate the barycentric coefficients of the point
    auto coeffs = tm.get_bary(point);

    for (int i = 0; i < DIM; i++) {
      int elem = tm.which(i, coeffs);
      if (elem >= 0 &&
          dimensionalities(point_ind) > (TreePointSearch::Dimensionality)i) {
        auto elem_ind = adjacencies[i][3 * val.index + elem];
        dimensionalities(point_ind) = (TreePointSearch::Dimensionality)i;
        element_ids(point_ind) = elem_ind;
        for (int j = 0; j < DIM + 1; j++) {
          parametric_coords(point_ind, j) = coeffs(j);
        }
        return;
      }
    }
    if (tm.which(DIM, coeffs) >= 0 &&
        dimensionalities(point_ind) > TreePointSearch::Dimensionality::FACE) {
      dimensionalities(point_ind) = TreePointSearch::Dimensionality::FACE;
      element_ids(point_ind) = (LO)val.index;
      for (int j = 0; j < DIM + 1; j++) {
        parametric_coords(point_ind, j) = coeffs(j);
      }
    }
  }

private:
  Rank1View<const Mapping<3>, MemorySpace> mappings;
  Omega_h::LOs adjacencies[3];
  Rank1View<TreePointSearch::Dimensionality, MemorySpace> dimensionalities;
  Rank1View<LO, MemorySpace> element_ids;
  Rank2View<Real, MemorySpace> parametric_coords;
};

class CallOnIntersect2D
{
public:
  using MemorySpace = TreePointSearch::MemorySpace;
  static constexpr int DIM = 2;

  CallOnIntersect2D(const Kokkos::View<Mapping<2>*, MemorySpace>& mappings_,
                    const Omega_h::LOs& adjacencies0,
                    const Omega_h::LOs& adjacencies1,
                    Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace>&
                      dimensionalities_,
                    Kokkos::View<LO*, MemorySpace>& element_ids_,
                    Kokkos::View<Real**, MemorySpace>& parametric_coords_)
    : mappings(MakeRank1View(mappings_)),
      adjacencies{adjacencies0, adjacencies1},
      dimensionalities(MakeRank1View(dimensionalities_)),
      element_ids(MakeRank1View(element_ids_)),
      parametric_coords(MakeRank2View(parametric_coords_)) {};

  template <typename Predicate, typename Value>
  KOKKOS_FUNCTION void operator()(Predicate const& predicate,
                                  Value const& val) const
  {
    ArborX::Point<DIM, Omega_h::Real> const& ax =
      ArborX::getGeometry(predicate);
    Omega_h::Vector<DIM> point{ax[0], ax[1], ax[2]};
    int point_ind = ArborX::getData(predicate);

    detail::Mapping<3> const& tm = mappings(val.index);

    // calculate the barycentric coefficients of the point
    auto coeffs = tm.get_bary(point);
    int offsets[DIM] = {4, 6, 4};
    for (int i = 0; i < DIM; i++) {
      int elem = tm.which(i, coeffs);
      if (elem >= 0 &&
          dimensionalities(point_ind) > (TreePointSearch::Dimensionality)i) {
        auto elem_ind = adjacencies[i][offsets[i] * val.index + elem];
        dimensionalities(point_ind) = (TreePointSearch::Dimensionality)i;
        element_ids(point_ind) = elem_ind;
        for (int j = 0; j < DIM + 1; j++) {
          parametric_coords(point_ind, j) = coeffs(j);
        }
        return;
      }
    }

    if (tm.which(DIM, coeffs) >= 0 &&
        dimensionalities(point_ind) > TreePointSearch::Dimensionality::REGION) {
      dimensionalities(point_ind) = TreePointSearch::Dimensionality::REGION;
      element_ids(point_ind) = (LO)val.index;
      for (int j = 0; j < DIM + 1; j++) {
        parametric_coords(point_ind, j) = coeffs(j);
      }
    }
  }

private:
  Rank1View<const Mapping<2>, MemorySpace> mappings;
  Omega_h::LOs adjacencies[2];
  Rank1View<TreePointSearch::Dimensionality, MemorySpace> dimensionalities;
  Rank1View<LO, MemorySpace> element_ids;
  Rank2View<Real, MemorySpace> parametric_coords;
};
} // namespace detail

TreePointSearch::Results TreePointSearch::Apply(
  const CoordinateView<TreePointSearch::MemorySpace>& coords) const
{
  if (coords.GetCoordinateSystem() != pcms::CoordinateSystem::Cartesian) {
    throw pcms_error("TreePointSearch::Apply only implemented for"
                     " Cartesian coordinates");
  }
  if (coords.GetValues().extent(1) != mesh_.dim()) {
    throw pcms_error("Input coordinate space dimension " +
                     std::to_string(coords.GetValues().extent(1)) +
                     " does not match query space dimension " +
                     std::to_string(mesh_.dim()));
  }

  if (mesh_.dim() == 2) {
    static constexpr int DIM = 2;
    Omega_h::ExecSpace execution_space;

    Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace> dims(
      "dimensionalities", coords.GetValues().extent(0));
    Kokkos::deep_copy(dims, TreePointSearch::Dimensionality::NO_INTERSECT);

    Kokkos::View<LO*, MemorySpace> elem_ids("element IDs",
                                            coords.GetValues().extent(0));
    Kokkos::deep_copy(elem_ids, -1);

    Kokkos::View<Real**, MemorySpace> parametric_coords(
      "parametric coordinates", coords.GetValues().extent(0),
      coords.GetValues().extent(1) + 1);
    Kokkos::deep_copy(parametric_coords, -1.0);

    tree->get_tree<2>().query(
      execution_space,
      detail::Coordinate_View_Adapt<MemorySpace, DIM>{coords.GetValues()},
      detail::CallOnIntersect2D(tree->get_mappings<2>(),
                                mesh_.get_adj(Omega_h::FACE, 0).ab2b,
                                mesh_.get_adj(Omega_h::FACE, 1).ab2b, dims,
                                elem_ids, parametric_coords));
    return Results{dims, elem_ids, parametric_coords};
  } else {
    static constexpr int DIM = 3;
    Omega_h::ExecSpace execution_space;

    Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace> dims(
      "dimensionalities", coords.GetValues().extent(0));
    Kokkos::deep_copy(dims, TreePointSearch::Dimensionality::NO_INTERSECT);

    Kokkos::View<LO*, MemorySpace> elem_ids("element IDs",
                                            coords.GetValues().extent(0));
    Kokkos::deep_copy(elem_ids, -1);

    Kokkos::View<Real**, MemorySpace> parametric_coords(
      "parametric coordinates", coords.GetValues().extent(0),
      coords.GetValues().extent(1) + 1);
    Kokkos::deep_copy(parametric_coords, -1.0);

    tree->get_tree<3>().query(
      execution_space,
      detail::Coordinate_View_Adapt<MemorySpace, DIM>{coords.GetValues()},
      detail::CallOnIntersect3D(tree->get_mappings<3>(),
                                mesh_.get_adj(Omega_h::REGION, 0).ab2b,
                                mesh_.get_adj(Omega_h::REGION, 1).ab2b,
                                mesh_.get_adj(Omega_h::REGION, 2).ab2b, dims,
                                elem_ids, parametric_coords));
    return Results{dims, elem_ids, parametric_coords};
  }
}

[[nodiscard]] LO TreePointSearch::GetOwningElementId(
  const TreePointSearch::Results& results, int i)
{
  const Kokkos::View<LO[1]> query_id{""};
  const Kokkos::View<Dimensionality[1]> dim{""};
  Kokkos::parallel_for(
    1, KOKKOS_LAMBDA(const int) {
      query_id(0) = (results.element_ids(i) < 0) ? -results.element_ids(i)
                                                 : results.element_ids(i);
      dim(0) = results.dimensionalities(i);
    });
  auto query_id_h = Kokkos::create_mirror_view(query_id);
  auto dim_h = Kokkos::create_mirror_view(dim);
  Kokkos::deep_copy(query_id_h, query_id);
  Kokkos::deep_copy(dim_h, dim);
  return pcms::GetOwningElementId(mesh_, mesh_.dim(),
                                  static_cast<int>(dim_h(0)), query_id_h(0));
}

[[nodiscard]] Kokkos::View<LO*> TreePointSearch::GetOwningElementIds(
  const TreePointSearch::Results& results)
{
  Kokkos::View<LO*> owning_ids("Owning element IDs",
                               results.dimensionalities.size());
  auto vert2elem = mesh_.ask_up(Omega_h::VERT, mesh_.dim());
  auto edge2elem = mesh_.ask_up(Omega_h::EDGE, mesh_.dim());
  auto face2elem = (mesh_.dim() == 3) ? mesh_.ask_up(Omega_h::FACE, mesh_.dim())
                                      : Omega_h::Adj{};
  auto mesh_dim = mesh_.dim();
  Kokkos::parallel_for(
    owning_ids.size(), KOKKOS_LAMBDA(const int i) {
      if (static_cast<int>(results.dimensionalities(i)) == mesh_dim) {
        owning_ids(i) = abs(results.element_ids(i));
      } else if (static_cast<int>(results.dimensionalities(i)) == 2) {
        owning_ids(i) = GetOwningElementIdFromAdj(face2elem, 2, mesh_dim,
                                                  results.element_ids(i));
      } else if (static_cast<int>(results.dimensionalities(i)) == 1) {
        owning_ids(i) = GetOwningElementIdFromAdj(edge2elem, 1, mesh_dim,
                                                  results.element_ids(i));
      } else if (static_cast<int>(results.dimensionalities(i)) == 0) {
        owning_ids(i) = GetOwningElementIdFromAdj(vert2elem, 0, mesh_dim,
                                                  results.element_ids(i));
      } else {
        owning_ids(i) = -1;
      }
    });
  return owning_ids;
}

std::unique_ptr<detail::TreeWrapper> TreePointSearch::make_tree(
  Omega_h::Mesh& mesh) const
{
  ExecSpace execution_space;
  using DeviceType =
    Kokkos::Device<Omega_h::ExecSpace, Omega_h::ExecSpace::memory_space>;

  if (mesh.dim() == 2) {
    detail::Omega_h_Mesh_Adapt<2> tagged_mesh{
      mesh.nelems(), mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b,
      mesh.coords()};

    detail::TreeWrapper::Tree_t<2> tree = detail::TreeWrapper::Tree_t<2>(
      execution_space, ArborX::Experimental::attach_indices(tagged_mesh));

    detail::TreeWrapper::Mappings_t<2> mappings =
      detail::TreeWrapper::Mappings_t<2>("mappings", mesh.nelems());

    auto face2vert =
      Omega_h::HostRead(mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b);
    auto vert_coords = Omega_h::HostRead(mesh.coords());

    auto mappings_h = Kokkos::create_mirror_view(mappings);
    for (int i = 0; i < mesh.nelems(); i++) {
      Omega_h::Matrix<2, 3> triangle = {
        {vert_coords[face2vert[i * 3] * 2],
         vert_coords[face2vert[i * 3] * 2 + 1]},
        {vert_coords[face2vert[i * 3 + 1] * 2],
         vert_coords[face2vert[i * 3 + 1] * 2 + 1]},
        {vert_coords[face2vert[i * 3 + 2] * 2],
         vert_coords[face2vert[i * 3 + 2] * 2 + 1]}};
      mappings_h[i] = detail::Mapping<2>(triangle, tolerances_);
    }
    Kokkos::deep_copy(execution_space, mappings, mappings_h);
    return std::make_unique<detail::TreeWrapper>(mappings, tree);
  }
  if (mesh.dim() == 3) {
    detail::Omega_h_Mesh_Adapt<3> tagged_mesh{
      mesh.nelems(), mesh.ask_down(Omega_h::REGION, Omega_h::VERT).ab2b,
      mesh.coords()};
    detail::TreeWrapper::Tree_t<3> tree = detail::TreeWrapper::Tree_t<3>(
      execution_space, ArborX::Experimental::attach_indices(tagged_mesh));

    detail::TreeWrapper::Mappings_t<3> mappings =
      detail::TreeWrapper::Mappings_t<3>("mappings", mesh.nelems());

    auto region2vert =
      Omega_h::HostRead(mesh.ask_down(Omega_h::REGION, Omega_h::VERT).ab2b);
    auto vert_coords = Omega_h::HostRead(mesh.coords());

    auto mappings_h = Kokkos::create_mirror_view(mappings);
    for (int i = 0; i < mesh.nelems(); i++) {
      Omega_h::Matrix<3, 4> tetrahedron = {
        {vert_coords[region2vert[i * 4] * 3],
         vert_coords[region2vert[i * 4] * 3 + 1],
         vert_coords[region2vert[i * 4] * 3 + 2]},
        {vert_coords[region2vert[i * 4 + 1] * 3],
         vert_coords[region2vert[i * 4 + 1] * 3 + 1],
         vert_coords[region2vert[i * 4 + 1] * 3 + 2]},
        {vert_coords[region2vert[i * 4 + 2] * 3],
         vert_coords[region2vert[i * 4 + 2] * 3 + 1],
         vert_coords[region2vert[i * 4 + 2] * 3 + 2]},
        {vert_coords[region2vert[i * 4 + 3] * 3],
         vert_coords[region2vert[i * 4 + 3] * 3 + 1],
         vert_coords[region2vert[i * 4 + 3] * 3 + 2]}};

      auto tol_h = Kokkos::create_mirror_view(tolerances_);
      Kokkos::deep_copy(tol_h, tolerances_);
      mappings_h[i] = detail::Mapping<3>(tetrahedron, tolerances_);
    }
    Kokkos::deep_copy(execution_space, mappings, mappings_h);
    return std::make_unique<detail::TreeWrapper>(mappings, tree);
  }
  throw pcms_error("Invalid mesh dimension " + std::to_string(mesh.dim()) +
                   ", TreePointSearch only implemented for 2D and 3D");
}

} // namespace pcms
