#ifndef POINT_LOCALIZATION_H
#define POINT_LOCALIZATION_H

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
#include "pcms/field/coordinate_system.h"

namespace pcms
{

class Mapping2D
{
public:
	static constexpr int DIM = 2;
	using MemorySpace = Omega_h::ExecSpace::memory_space;
	Mapping2D() = default;
	Mapping2D(int elem_index, Omega_h::Mesh const& mesh);
	~Mapping2D() = default;
	KOKKOS_FUNCTION
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const;
	KOKKOS_FUNCTION
	int which(int ent_dim, 
					  Omega_h::Vector<DIM + 1> const& bary_coords) const;
private:
	// helpers
	int which_vert(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	int which_edge(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	int within_elem(Omega_h::Vector<DIM + 1> const& bary_coords) const;
	void set_mesh_triangle(int index, Omega_h::Mesh const& mesh);
	double opposite_edge_len_sq(int i) const;
	// representation
	Omega_h::Matrix<DIM,DIM> bary_transform; // Column-major order
	Omega_h::Matrix<DIM,DIM + 1> triangle;
	double triangle_area;
};

class Mapping3D
{
public:
	static constexpr int DIM = 3;
	using MemorySpace = Omega_h::ExecSpace::memory_space;
	Mapping3D() = default;
	Mapping3D(int elem_index, Omega_h::Mesh const& mesh);
	~Mapping3D() = default;
	KOKKOS_FUNCTION
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const;
	KOKKOS_FUNCTION
	int which(int ent_dim, 
					  Omega_h::Vector<DIM+1> const& bary_coords) const;
private:
	// helpers
	int which_vert(Omega_h::Vector<DIM+1> const& bary_coords) const;
	int which_edge(Omega_h::Vector<DIM+1> const& bary_coords) const;
	int which_face(Omega_h::Vector<DIM+1> const& bary_coords) const;
	int within_elem(Omega_h::Vector<DIM+1> const& bary_coords) const;
	void set_mesh_tet(int index, Omega_h::Mesh const& mesh);
	void set_triangle_areas();
	// returns the actual face/edge offset from the barycentric offset
	KOKKOS_INLINE_FUNCTION int face(int i) const 
	{ return (OFFSETS & 3 << (i*2))>>(i*2);}
	KOKKOS_INLINE_FUNCTION int edge(int i) const 
	{ return (i < 4 && i > 0) ? (OFFSETS & 3 << ((i-1)*2))>>((i-1)*2) : i; }
	// representation
	Omega_h::Matrix<DIM,DIM> bary_transform; // Column-major order
	Omega_h::Matrix<DIM,DIM+1> tetrahedron;
	Omega_h::Vector<DIM+1> face_areas;
	double tetrahedron_volume;
	// This is used to calculate the face and edge offsets in the cannonical ordering: 
	// https://user-images.githubusercontent.com/56453280/74203616-39d27a80-4c3e-11ea-885d-b0260490e184.png
	static const char OFFSETS = 1 << 4 | 3 << 2 | 2;
};

struct TreeWrapper
{
	virtual void* get_tree() = 0;
	virtual void* get_mappings() = 0;
};
struct TreeWrapper2D : TreeWrapper
{
	using Mappings_t = Kokkos::View<Mapping2D*, Omega_h::ExecSpace::memory_space>;
	using Tree_t = ArborX::BVH<Omega_h::ExecSpace::memory_space,
			ArborX::PairValueIndex<ArborX::Box<2, double>, unsigned>>;
	
	TreeWrapper2D(const Mappings_t& mappings, const Tree_t& tree) : mappings_(mappings), tree_(tree) {}
	~TreeWrapper2D() = default;
	void* get_tree() override { return &tree_; };
	void* get_mappings() override {return &mappings_; };
private:
	Tree_t tree_;
	Mappings_t mappings_;
};
struct TreeWrapper3D : TreeWrapper
{
	using Mappings_t = Kokkos::View<Mapping3D*, Omega_h::ExecSpace::memory_space>;
	using Tree_t = ArborX::BVH<Omega_h::ExecSpace::memory_space,
			ArborX::PairValueIndex<ArborX::Box<3, double>, unsigned>>;

	TreeWrapper3D(const Mappings_t& mappings, const Tree_t& tree) : mappings_(mappings), tree_(tree) {}
	~TreeWrapper3D() = default;
	void* get_tree() override { return &tree_; };
	void* get_mappings() override {return &mappings_; };
private:
	Tree_t tree_;
	Mappings_t mappings_;
};

/**
 * Point search base class
 */
class PointSearch
{
public:
	using ExecSpace = Omega_h::ExecSpace;
	using MemorySpace = ExecSpace::memory_space;
	struct Result
	{
		enum class Dimensionality
		{
			VERTEX = 0,
			EDGE = 1,
			FACE = 2,
			REGION = 3,
			NO_INTERSECT = 4
		};

		Dimensionality dimensionality;
		LO element_id;
		Omega_h::Vector<4> parametric_coords;
	};

	PointSearch() = default;
	~PointSearch() = default;
	virtual Kokkos::View<Result*> apply(const CoordinateView<MemorySpace>& coords) const = 0;
};

class TreePointSearch : public PointSearch
{
public:
	using Result = PointSearch::Result;
	using Dimensionality = Result::Dimensionality;
	using ExecSpace = PointSearch::ExecSpace;
	using MemorySpace = PointSearch::MemorySpace;
	TreePointSearch(const Omega_h::Mesh& mesh);
	~TreePointSearch() = default;
	Kokkos::View<Result*> apply(
		const CoordinateView<Omega_h::ExecSpace::memory_space>& coords) const override;
private:
	std::unique_ptr<TreeWrapper> make_tree(const Omega_h::Mesh& mesh) const;
	KOKKOS_INLINE_FUNCTION Mapping2D get_mapping2D(LO i) const 
	{ return ((TreeWrapper2D::Mappings_t*)tree->get_mappings())->operator()(i); };
	KOKKOS_INLINE_FUNCTION Mapping3D get_mapping3D(LO i) const 
	{ return ((TreeWrapper3D::Mappings_t*)tree->get_mappings())->operator()(i); };
	// Reference to the input mesh
	Omega_h::Mesh const &mesh_;
	std::unique_ptr<TreeWrapper> tree;
};

} // namespace pcms
#endif // POINT_LOCALIZATION_H
