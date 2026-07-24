#ifndef POINT_LOCALIZATION_H
#define POINT_LOCALIZATION_H

#include <any>

#include <ArborX.hpp>
#include <ArborX_Triangle.hpp>
#include <detail/ArborX_PairValueIndex.hpp>
#include <detail/ArborX_AttachIndices.hpp>

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_simplex.hpp>

#include "pcms/utility/assert.h"
#include "pcms/utility/types.h"
#include "pcms/utility/arrays.h"
#include "pcms/localization/point_search.h"

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
LO kronecker(LO i, LO j) { return (LO)(i==j); }

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
	const Rank2View<const Real, MemorySpace, default_layout_for_memory_space_t<MemorySpace>> points;
};

template <int dim>
class Mapping {};

/**
* @brief Mapping<2> represents a mapping from the global coordinate system
*		 to the barycentric coordinate system of an element in a 2D Omega_h::Mesh
* The class Mapping<2> calculates the mapping of any point in global coordinate
* space to the barycentric coordinate space belonging to a triangle (face) in
* an Omega_h::Mesh.
* An instance of Mapping<2> can, from the barycentric coordinates constructed by it,
* determine whether that point intersects the face, an edge, a vertex, or lies 
* outside the face. It also determines which edge or vertex the point intersects
* by calculating the offset in the Omega_h::Adj::ab2b
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
	* @param elem_index the index of the desired triangle or tetrahedron in the Omega_h::Mesh
	* @param mesh the Omega_h::Mesh with the spatial information of the triangle
	*/

	Mapping(Omega_h::Matrix<DIM,DIM+1> const& triangle_, 
		Kokkos::View<Real*> const& tolerances);
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
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const;
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
	int which(
		int ent_dim, 
		Omega_h::Vector<DIM + 1> const& bary_coords) const;
private:
	/**
	* @brief Computes the vertex offset of the point corresponding 
	* 		  to the input barycentric coordinates
	* @param bary_coords the input barycentric coordinates
	* @returns the vertex offset if the point corresponding to the given bary. coordinates
	* 			is within a certain (global) tolerance of a vertex
	* 			-1 if the point does not lie within the tolerance of any vertex
	*/
	KOKKOS_FUNCTION
	int which_vert(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	/**
	* @brief Computes the edge offset of the point corresponding 
	* 		 to the input barycentric coordinates
	* @param bary_coords the input barycentric coordinates
	* @returns the edge offset if the distnce from the point corresponding to the 
	* 			given bary. coordinates to the edge with the respective offset
	* 			is within a certain (global) tolerance
	* 			-1 if the point does not lie within the tolerance of any vertex
	*/
	KOKKOS_FUNCTION
	int which_edge(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	/**
	* @brief Computes whether a point with the given barycentric coordinates
	* 		  is within a triangle
	* @param bary_coords the input barycentric coordinates
	* @returns 0 if the point is within the highest order element and -1 otherwise
	*/
	KOKKOS_FUNCTION
	int within_elem(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	/**
	* @brief Constructs the triangle the mapping corresponds to
	* @param index the element ID of the triangle
	* @param mesh the Omega_h::Mesh with the spatial information for the triangle
	* 			   corresponding to index
	*/
	void set_mesh_triangle(int index, Omega_h::Mesh const& mesh);
	/**
	* @brief Calculates the length of the edge opposite vertex i
	* @param i the vertex index opposite the desired edge
	* @returns the length of the ith edge
	*/
	KOKKOS_FUNCTION
	double opposite_edge_len_sq(int i) const;
	
	// REPRESENTATION
	Omega_h::Vector<DIM> tolerances_;
	// representation of the barycentric coordinate mapping, source:
	// https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Edge_approach
	Omega_h::Matrix<DIM,DIM> bary_transform; // Column-major order
	// Points defining the triangle
	Omega_h::Matrix<DIM,DIM + 1> triangle;
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

	Mapping(Omega_h::Matrix<DIM,DIM+1> const& tetrahedron_, 
			const Kokkos::View<Real*>& tolerances);
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
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const;

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
	int which(int ent_dim, 
			Omega_h::Vector<DIM+1> const& bary_coords) const;
private:	
	/**
	* @brief Computes the vertex offset of the point corresponding 
	*		 to the input barycentric coordinates
	* @param bary_coords the input barycentric coordinates
	* @returns the vertex offset if the point corresponding to the given bary. coordinates
	* 		   is within a certain (global) tolerance of a vertex
	* 		   -1 if the point does not lie within the tolerance of any vertex
	*/
	KOKKOS_FUNCTION
	int which_vert(Omega_h::Vector<DIM+1> const& bary_coords) const;
	/**
	* @brief Computes the edge offset of the point corresponding 
	* 		 to the input barycentric coordinates
	* @param bary_coords the input barycentric coordinates
	* @returns the edge offset if the distnce from the point corresponding to the
	* 		   given bary. coordinates to the edge with the respective offset
	* 		   is within a certain (global) tolerance
	* 		   -1 if the point does not lie within the tolerance of any vertex
	*/
	KOKKOS_FUNCTION
	int which_edge(Omega_h::Vector<DIM+1> const& bary_coords) const;
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
	int which_face(Omega_h::Vector<DIM+1> const& bary_coords) const;
	/**
	* @brief Computes whether a point with the given barycentric coordinates
	* 		  is within a triangle
	* @param bary_coords the input barycentric coordinates
	* @returns 0 if the point is within the highet order element and -1 otherwise
	*/
	KOKKOS_FUNCTION
	int within_elem(Omega_h::Vector<DIM+1> const& bary_coords) const;
	/**
	* @brief Constructs the tetrahedron the mapping corresponds to
	* @param index the element ID of the tetrahedron
	* @param mesh the Omega_h::Mesh with the spatial information for the 
	* 			  tetrahedron corresponding to `index`
	*/
	void set_mesh_tet(int index, Omega_h::Mesh const& mesh);
	/**
	* @brief Calculates the areas of each face in the tetrahedron the mapping 
	* 		 corresponds to
	* @param index the element ID of the tetrahedron
	* @param mesh the Omega_h::Mesh with the spatial information for the 
	* 			  tetrahedron corresponding to `index`
	*/
	void set_triangle_areas();
	// returns the actual face/edge offset from the barycentric offset
	KOKKOS_INLINE_FUNCTION int face(int i) const 
	{ return (OFFSETS & 3 << (i*2))>>(i*2); }
	KOKKOS_INLINE_FUNCTION int edge(int i) const 
	{ return (i < 4 && i > 0) ? (OFFSETS & 3 << ((i-1)*2))>>((i-1)*2) : i; }
	// representation
	Omega_h::Vector<DIM> tolerances_;
	Omega_h::Matrix<DIM,DIM> bary_transform; // Column-major order
	Omega_h::Matrix<DIM,DIM+1> tetrahedron;
	Omega_h::Vector<DIM+1> face_areas;
	double tetrahedron_volume;
	// This is used to calculate the face and edge offsets in the cannonical ordering: 
	// https://user-images.githubusercontent.com/56453280/74203616-39d27a80-4c3e-11ea-885d-b0260490e184.png
	static const char OFFSETS = 1 << 4 | 3 << 2 | 2;
};

/**
* Wrapper base class for the ArborX::BVH and Mappings classes, this is a *rough*
* workaround for ArborX::BVH being templated on dimension
*/
struct TreeWrapper
{
	template <int dim>
	using Mappings_t = Kokkos::View<Mapping<dim>*, Omega_h::ExecSpace::memory_space>;
	template <int dim>
	using Tree_t = ArborX::BVH<Omega_h::ExecSpace::memory_space,
			ArborX::PairValueIndex<ArborX::Box<dim, double>, unsigned>>;
	
	template <int dim>
	TreeWrapper(const Mappings_t<dim>& mappings, const Tree_t<dim>& tree)
	{
		dim_ = dim;
		mappings_ = std::make_shared<std::any>(mappings);
		tree_ = std::make_shared<std::any>(tree);
	}
	/**
	* @brief returns a pointer to an ArborX tree
	*/
	template <int dim>
	inline Tree_t<dim> get_tree()
	{
		if (dim != dim_)
		{
			throw pcms_error("Requested get_tree return "
				"type does not match internal Tree_t");
		}
		return std::any_cast<Tree_t<dim>>(*tree_);
	}

	/**
	* @brief returns a pointer to a Kokkos::View of Mappings
	*/
	template <int dim>
	inline Mappings_t<dim> get_mappings()
	{
		if (dim != dim_)
		{
			throw pcms_error("Requested get_mappings return "
				"type does not match internal Mappings_t");
		}
		return std::any_cast<Mappings_t<dim>>(*mappings_);
	};
private:
	int dim_;
	std::shared_ptr<std::any> tree_;
	std::shared_ptr<std::any> mappings_;
};
} //namespace detail

class TreePointSearch : public PointSearch
{
public:
	using Results = PointSearch::Results;
	using Dimensionality = Results::Dimensionality;
	using ExecSpace = PointSearch::ExecSpace;
	using MemorySpace = PointSearch::MemorySpace;
	
	TreePointSearch(Omega_h::Mesh& mesh) : 
		PointSearch(PointSearchTolerances("tree point search tolerances", mesh.dim())), 
		mesh_(mesh)
	{
		Kokkos::deep_copy(tolerances_, 1e-12);
		tree = make_tree(mesh);
	}
	TreePointSearch(Omega_h::Mesh& mesh, const PointSearchTolerances& tolerances) : PointSearch(tolerances), mesh_(mesh), tree(make_tree(mesh)) {}
	~TreePointSearch() = default;
	/**
	* Given a set of points in global coordinates give the ids of the entities
	* that the points lie within and the parametric coordinates of each point 
	* within a triangles or tetrahedra adjacent to the intersected entity. 
	* If the point does not lie within any triangle element, the id will 
	* be a negative number.
	*/
	Results apply(const CoordinateView<MemorySpace>& coords) const override;
	[[nodiscard]] LO GetOwningElementId(const Results& results, int i) override;
	[[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
		const Results& results) override;
private:
	std::unique_ptr<detail::TreeWrapper> make_tree(const Omega_h::Mesh& mesh) const;
	// Reference to the input mesh
	Omega_h::Mesh &mesh_;
	std::unique_ptr<detail::TreeWrapper> tree;
};

namespace detail
{

	/**
	* Functor to be called when a point intersects a triangle or tetrahedron
	*/
	class CallOnIntersect3D
	{
	public:
		using MemorySpace = TreePointSearch::MemorySpace;
		static constexpr int DIM = 3;

		CallOnIntersect3D(
			const Kokkos::View<Mapping<3>*, MemorySpace>& mappings_,
			const Omega_h::LOs& adjacencies0,
			const Omega_h::LOs& adjacencies1,
			const Omega_h::LOs& adjacencies2,
			Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace>& dimensionalities_,
			Kokkos::View<LO*, MemorySpace>& element_ids_,
			Kokkos::View<Real**, MemorySpace>& parametric_coords_
		) : mappings(MakeRank1View(mappings_)),  
			adjacencies{adjacencies0, adjacencies1, adjacencies2},
			dimensionalities(MakeRank1View(dimensionalities_)),
			element_ids(MakeRank1View(element_ids_)),
			parametric_coords(MakeRank2View(parametric_coords_))
		{
		};
		
		/**
		* Intersection callback, 
		*/
		template <typename Predicate, typename Value>
		KOKKOS_FUNCTION void operator()(Predicate const &predicate, Value const & val) const;
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

		CallOnIntersect2D(
			const Kokkos::View<Mapping<2>*, MemorySpace>& mappings_,
			const Omega_h::LOs& adjacencies0,
			const Omega_h::LOs& adjacencies1,
			Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace>& dimensionalities_,
			Kokkos::View<LO*, MemorySpace>& element_ids_,
			Kokkos::View<Real**, MemorySpace>& parametric_coords_
		) : mappings(MakeRank1View(mappings_)),  
			adjacencies{adjacencies0, adjacencies1},
			dimensionalities(MakeRank1View(dimensionalities_)),
			element_ids(MakeRank1View(element_ids_)),
			parametric_coords(MakeRank2View(parametric_coords_))
		{
		};
		
		template <typename Predicate, typename Value>
		KOKKOS_FUNCTION void operator()(Predicate const &predicate, Value const & val) const;
	private:
		Rank1View<const Mapping<2>, MemorySpace> mappings;
		Omega_h::LOs adjacencies[2];
		Rank1View<TreePointSearch::Dimensionality, MemorySpace> dimensionalities;
		Rank1View<LO, MemorySpace> element_ids;
		Rank2View<Real, MemorySpace> parametric_coords;
	};

} // namespace detail

} // namespace pcms
#endif // POINT_LOCALIZATION_H
