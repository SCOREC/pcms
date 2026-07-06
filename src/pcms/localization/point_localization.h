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
#include "pcms/field/coordinate_system.h"

namespace pcms
{

template <int dim>
class Mapping
{
public:
	static constexpr int DIM = dim;
	virtual Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const = 0;
	virtual int which(int ent_dim, 
					  Omega_h::Vector<DIM + 1> const& bary_coords) const = 0;
};

// Two dimensional mapping class
// Constructs a mapping of any point in global space to
// the barycentric coordinate system of a given triangle
class Mapping2D : public Mapping<2>
{
public:
	Mapping2D() = default;
	Mapping2D(int elem_index, Omega_h::Mesh const& mesh);
	~Mapping2D() = default;
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const override;
	int which(int ent_dim, 
					  Omega_h::Vector<DIM + 1> const& bary_coords) const override;
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

class Mapping3D : public Mapping<3>
{
public:
	Mapping3D() = default;
	Mapping3D(int elem_index, Omega_h::Mesh const& mesh);
	~Mapping3D() = default;
	Omega_h::Vector<DIM + 1> get_bary(Omega_h::Vector<DIM> const& p) const override;
	int which(int ent_dim, 
					  Omega_h::Vector<DIM+1> const& bary_coords) const override;
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

/**
 * Point search base class
 */
template <int dim>
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
			REGION = 3
		};

		Dimensionality dimensionality;
		LO element_id;
		Omega_h::Vector<dim + 1> parametric_coords;
	};

	static constexpr auto DIM = dim;
	
	PointSearch() = default;
	~PointSearch() = default;
	virtual Kokkos::View<Result*> apply(const CoordinateView<MemorySpace>& coords) const = 0;
};

using PointSearch2D = PointSearch<2>;
using PointSearch3D = PointSearch<3>;

class TreePointSearch2D : public PointSearch2D
{
public:
	using Result_t = PointSearch2D::Result;
	using Dim_t = Result_t::Dimensionality;
	TreePointSearch2D(const Omega_h::Mesh& mesh);
	Kokkos::View<Result_t*> apply(
		const CoordinateView<Omega_h::ExecSpace::memory_space>& coords) const override;
private:
	// Reference to the input mesh
	Omega_h::Mesh const &mesh_;
	// Mapping for each triangle in the mesh
	// (TODO find way to make these the leaf nodes of the tree)
	Kokkos::View<Mapping2D*> mappings;
	// Bounding Volume Hierarchy of input mesh
	ArborX::BVH<Omega_h::ExecSpace::memory_space,
				ArborX::PairValueIndex<ArborX::Box<2, double>, unsigned>> tree;
};

class TreePointSearch3D : public PointSearch3D
{
public:
	using Result_t = PointSearch3D::Result;
	using Dim_t = Result_t::Dimensionality;
	TreePointSearch3D(const Omega_h::Mesh& mesh);
	Kokkos::View<Result_t*> apply(
		const CoordinateView<Omega_h::ExecSpace::memory_space>& coords) const override;
private:
	// Reference to the input mesh
	Omega_h::Mesh const &mesh_;
	// Mapping for each triangle in the mesh
	// (TODO find way to make these the leaf nodes of the tree)
	Kokkos::View<Mapping3D*> mappings;
	// Bounding Volume Hierarchy of input mesh
	ArborX::BVH<Omega_h::ExecSpace::memory_space,
				ArborX::PairValueIndex<ArborX::Box<3, double>, unsigned>> tree;
};

} // namespace pcms
#endif // POINT_LOCALIZATION_H
