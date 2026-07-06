#include "point_localization.h"

#define FACE_TOL 10e-8
#define EDGE_TOL 10e-6
#define VERT_TOL 10e-5

#define kronecker(i,j) (int)(i==j)

// adds template parameters required by ArborX Access traits
template <typename DeviceType, int dim, class coord = float>
struct Omega_h_Mesh_Tagged
{
	const Omega_h::Mesh & m;
};

// Access traits for Omega_h
// constructs bounding boxes for elements
template <typename DeviceType, int dim, class coord>
struct ArborX::AccessTraits<Omega_h_Mesh_Tagged<DeviceType, dim, coord>>
{
	using memory_space = typename DeviceType::memory_space;

	static KOKKOS_FUNCTION int size(const Omega_h_Mesh_Tagged<DeviceType, dim, coord>& mesh)
	{
		return mesh.m.nelems();
	}

	static KOKKOS_FUNCTION auto get(const Omega_h_Mesh_Tagged<DeviceType, dim, coord>& mesh, int i)
	{
		const auto face2verts = mesh.m.get_adj(dim, Omega_h::VERT).ab2b;
		const auto vert_coords = mesh.m.coords();
		ArborX::Point<dim, coord> min = {INFINITY};
		ArborX::Point<dim, coord> max = {-INFINITY};
		for (int j = 0; j < dim + 1; ++j)
		{
			auto cell_vert_id = face2verts[(dim+1)*i+j];
			for (int k = 0; k < dim; ++k)
			{
				coord curr_coord = vert_coords[cell_vert_id*dim+k];
				if (min[k] > curr_coord) min[k] = curr_coord;
				if (max[k] < curr_coord) max[k] = curr_coord;
			}
		}
		return ArborX::Box(min, max);
	}
};

// adds template parameters required by ArborX Access traits
template <typename MemorySpace, int dim>
struct pcms_Coordinate_View_Tagged
{
	const pcms::CoordinateView<MemorySpace>& cv;
};

template <typename MemorySpace, int DIM>
struct ArborX::AccessTraits<pcms_Coordinate_View_Tagged<MemorySpace, DIM>>
{
	using memory_space = MemorySpace;
	static KOKKOS_FUNCTION int size(pcms_Coordinate_View_Tagged<MemorySpace, DIM> const &coords)
	{
		return coords.cv.GetCoordinates().extent(0);
	}
	static KOKKOS_FUNCTION auto get(
		pcms_Coordinate_View_Tagged<MemorySpace, DIM> const &coords, 
		int i)
	{
		const auto points = coords.cv.GetCoordinates();
		constexpr int dim = DIM;
		ArborX::Point<dim, double> ax_point;
		for (int j = 0; j < dim; j++) ax_point[j] = points(i,j);
		return PredicateWithAttachment(intersects(ax_point), i);
	}
};

namespace pcms
{

/**
* @brief Constructs a barycentric mapping of a triangle in an Omega_h::Mesh
* @param elem_index the index of the desired triangle or tetrahedron in the Omega_h::Mesh
* @param mesh the Omega_h::Mesh with the spatial information of the triangle
*/
Mapping2D::Mapping2D(int elem_index, Omega_h::Mesh const& mesh)
{
	set_mesh_triangle(elem_index, mesh);
	bary_transform = { triangle[0] - triangle[2], triangle[1] - triangle[2] };
	triangle_area = 0.5*fabs(Omega_h::determinant(bary_transform));
	bary_transform = Omega_h::invert(bary_transform);
}

/**
	* @brief Computes the barycentric coordinates of a point in global space
	* @param p the point to compute the barycentric coordinates of
	* @returns an Omega_h::Vector<dim + 1> containing the barycentric coordinates
	* 			of p
	*/
KOKKOS_FUNCTION
Omega_h::Vector<Mapping2D::DIM + 1> Mapping2D::get_bary(Omega_h::Vector<Mapping2D::DIM> const& p) const
{
	Omega_h::Vector<2> coeffs = bary_transform*(p - triangle[2]);
	return {coeffs[0], coeffs[1], 1 - coeffs[0] - coeffs[1]};
}

int Mapping2D::which(int dim, 
					 Omega_h::Vector<Mapping2D::DIM + 1> const& bary_coords) const
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

/**
	* @brief Computes the vertex offset of the point corresponding 
	* 		  to the input barycentric coordinates
	* @param bary_coords the input barycentric coordinates
	* @returns the vertex offset if the point corresponding to the given bary. coordinates
	* 			is within a certain (global) tolerance of a vertex
	* 			-1 if the point does not lie within the tolerance of any vertex
	*/
int Mapping2D::which_vert(Omega_h::Vector<Mapping2D::DIM + 1> const& bary_coords) const
{
	for (int i = 0; i < 3; i++)
	{
		// distance of the point to vertex i
		Omega_h::Vector<2> error = (bary_coords[0] - kronecker(i,0)) * triangle[0] 
									+ (bary_coords[1] - kronecker(i,1)) * triangle[1] 
									+ (bary_coords[2] - kronecker(i,2)) * triangle[2];
		if (Omega_h::norm_squared(error) <= VERT_TOL*VERT_TOL) return i;
	}
	return -1;
}

/**
* @brief Computes the edge offset of the point corresponding 
* 		 to the input barycentric coordinates
* @param bary_coords the input barycentric coordinates
* @returns the edge offset if the distnce from the point corresponding to the 
* 			given bary. coordinates to the edge with the respective offset
* 			is within a certain (global) tolerance
* 			-1 if the point does not lie within the tolerance of any vertex
*/
int Mapping2D::which_edge(Omega_h::Vector<Mapping2D::DIM + 1> const& bary_coords) const
{
	for (int i = 0; i <= 2; i++)
	{
		// Distance of the point to the edge opposite the ith vertex,
		// see docs for derivation.
		double dist = (2*bary_coords[i]*triangle_area)*(2*bary_coords[i]*triangle_area);
		dist /= opposite_edge_len_sq(i);
		if (dist <= EDGE_TOL*EDGE_TOL &&
			bary_coords[(i+1)%3] >= 0 && bary_coords[(i+2)%3] >= 0) return (i+1)%3;
	}
	return -1;
}

/**
* @brief Computes whether a point with the given barycentric coordinates
* 		  is within a triangle
* @param bary_coords the input barycentric coordinates
* @returns 0 if the point is within the highet order element and -1 otherwise
*/
int Mapping2D::within_elem(Omega_h::Vector<Mapping2D::DIM + 1> const& bary_coords) const
{
	return -1 * (int)(bary_coords[0] > 0 && bary_coords[1] > 0 && bary_coords[2] > 0);
}

/**
* @brief Constructs the triangle the mapping is for
* @param index the element ID of the triangle
* @param mesh the Omega_h::Mesh with the spatial information for the triangle
* 			   corresponding to index
*/
void Mapping2D::set_mesh_triangle(int index, Omega_h::Mesh const& mesh)
{
	auto face2vert = mesh.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b;
	auto vert_coords = mesh.coords();
	triangle = Omega_h::Matrix<2,3>{{vert_coords[face2vert[index*3]*2], vert_coords[face2vert[index*3]*2 + 1]},
			{vert_coords[face2vert[index*3 + 1]*2], vert_coords[face2vert[index*3 + 1]*2 + 1]},
			{vert_coords[face2vert[index*3 + 2]*2], vert_coords[face2vert[index*3 + 2]*2 + 1]}};
}

/**
* @brief Calculates the length of the edge opposite vertex i
* @param i the vertex index opposite the desired edge
* @returns the length of the ith edge
*/
double Mapping2D::opposite_edge_len_sq(int i) const
{
	return Omega_h::norm_squared(triangle[(i+2)%3] - triangle[(i+1)%3]);
}

/**
* @brief Constructs a barycentric mapping of a triangle in an Omega_h::Mesh
* @param elem_index the index of the desired triangle or tetrahedron in the Omega_h::Mesh
* @param mesh the Omega_h::Mesh with the spatial information of the triangle
*/
Mapping3D::Mapping3D(int elem_index, Omega_h::Mesh const& mesh)
{
	set_mesh_tet(elem_index, mesh);
	set_triangle_areas();
	bary_transform = { tetrahedron[0] - tetrahedron[3], tetrahedron[1] - tetrahedron[3], tetrahedron[2] - tetrahedron[3] };
	tetrahedron_volume = fabs(Omega_h::determinant(bary_transform))/6.;
	bary_transform = Omega_h::invert(bary_transform);
}

/**
* @brief Computes the barycentric coordinates of a point in global space
* @param p the point to compute the barycentric coordinates of
* @returns an Omega_h::Vector<dim + 1> containing the barycentric coordinates
* 			of p
*/
Omega_h::Vector<Mapping3D::DIM + 1> Mapping3D::get_bary(Omega_h::Vector<Mapping3D::DIM> const& p) const
{
	Omega_h::Vector<3> coeffs = bary_transform*(p - tetrahedron[3]);
	return {coeffs[0], coeffs[1], coeffs[2], 1 - coeffs[0] - coeffs[1] - coeffs[2]};
}

int Mapping3D::which(int dim, 
					 Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
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

/**
* @brief Computes the vertex offset of the point corresponding 
		to the input barycentric coordinates
* @param bary_coords the input barycentric coordinates
* @returns the vertex offset if the point corresponding to the given bary. coordinates
* 			is within a certain (global) tolerance of a vertex
* 			-1 if the point does not lie within the tolerance of any vertex
*/
int Mapping3D::which_vert(Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
{
	for (int i = 0; i < 4; i++)
	{
		// distance to the ith vertex
		Omega_h::Vector<3> error = (bary_coords[0] - kronecker(i,0)) * tetrahedron[0] + (bary_coords[1] - kronecker(i,1)) * tetrahedron[1]
									+ (bary_coords[2] - kronecker(i,2)) * tetrahedron[2] + (bary_coords[3] - kronecker(i,3)) * tetrahedron[3];
		if (Omega_h::norm_squared(error) <= VERT_TOL*VERT_TOL) return i;
	}
	return -1;
}

/**
* @brief Computes the edge offset of the point corresponding 
* 		 to the input barycentric coordinates
* @param bary_coords the input barycentric coordinates
* @returns the edge offset if the distnce from the point corresponding to the 
* 			given bary. coordinates to the edge with the respective offset
* 			is within a certain (global) tolerance
* 			-1 if the point does not lie within the tolerance of any vertex
*/
int Mapping3D::which_edge(Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
{
	int edge_ = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = i+1; j < 4; j++)
		{
			Omega_h::Vector<3> side1 = {0,0,0}, side2 = {0,0,0}, side3 = tetrahedron[i] - tetrahedron[j];
			for (int k = 0; k < 4; k++)
			{
				side1 += (bary_coords[k] - kronecker(i,k)) * tetrahedron[k];
				side2 += (bary_coords[k] - kronecker(j,k)) * tetrahedron[k];
			}

			// the norm of the cross product of two vectors is twice the area of the triangle
			// those vectors form
			double distnce_sq = Omega_h::norm_squared(Omega_h::cross(side1, side2));
			distnce_sq /= Omega_h::norm_squared(side3);
			
			if (distnce_sq <= EDGE_TOL*EDGE_TOL && bary_coords[i] > 0 && bary_coords[j] > 0) return edge(edge_);
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
* 			global tolerance of a face, returns the offset of that face
* 			-1 if the point does not lie within the tolerance of any face
* Algorithm source:
* C.E. Passerello,
* Interference detection using barycentric coordinates,
* Mechanics Research Communications,
* Volume 9, Issue 6, 1982, Pages 373-378,
* https://doi.org/10.1016/0093-6413(82)90034-9.
*/
int Mapping3D::which_face(Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
{
	for (int i = 0; i < 4; i++)
	{
		if (fabs(3*tetrahedron_volume*bary_coords[i]/face_areas[i]) <= FACE_TOL
			&& bary_coords[(i+1)%4] >= 0 && bary_coords[(i+2)%4] >= 0 && bary_coords[(i+3)%4] >= 0)
		{
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
int Mapping3D::within_elem(Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
{
	return -1 * (int)(bary_coords[0] >= 0 && bary_coords[1] >= 0 && bary_coords[2] >= 0 && bary_coords[3] >= 0);
}

// constructs the tetrahedron from the mesh's spatial data
void Mapping3D::set_mesh_tet(int index, Omega_h::Mesh const& mesh)
{
	auto region2vert = mesh.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b;
	auto vert_coords = mesh.coords();
	tetrahedron = {{vert_coords[region2vert[index*4]*3], vert_coords[region2vert[index*4]*3 + 1], vert_coords[region2vert[index*4]*3 + 2]},
	{vert_coords[region2vert[index*4 + 1]*3], vert_coords[region2vert[index*4 + 1]*3 + 1], vert_coords[region2vert[index*4 + 1]*3 + 2]},
	{vert_coords[region2vert[index*4 + 2]*3], vert_coords[region2vert[index*4 + 2]*3 + 1], vert_coords[region2vert[index*4 + 2]*3 + 2]},
	{vert_coords[region2vert[index*4 + 3]*3], vert_coords[region2vert[index*4 + 3]*3 + 1], vert_coords[region2vert[index*4 + 3]*3 + 2]}};
}

// calculates the areas of each of the faces of the triangle
// the faces are "named" according to the vertex opposite
// (i.e., face 0 is defined by points 1, 2, 3)
void Mapping3D::set_triangle_areas()
{
	for (int i = 0; i < 4; i++)
	{
		Omega_h::Vector<3> edge0 = tetrahedron[(i+1)%4] - tetrahedron[(i+3)%4], edge1 = tetrahedron[(i+2)%4] - tetrahedron[(i+3)%4];
		Omega_h::Vector<3> cross = Omega_h::cross(edge0, edge1);
		face_areas[i] = 0.5*Omega_h::norm(cross);
	}
}

TreePointSearch2D::TreePointSearch2D(const Omega_h::Mesh& mesh) : mesh_(mesh)
{
	if (mesh.dim() != 2) 
	{
		throw pcms_error("Could not construct TreePointSearch2D, invalid mesh dimension");
	}

	Omega_h::ExecSpace execution_space;
	using DeviceType = Kokkos::Device<Omega_h::ExecSpace, 
									  Omega_h::ExecSpace::memory_space>;
	Omega_h_Mesh_Tagged<DeviceType, 2, double> tagged_mesh{mesh};
	tree = ArborX::BVH(execution_space, 
						ArborX::Experimental::attach_indices(tagged_mesh));

	mappings = Kokkos::View<Mapping2D*, Omega_h::ExecSpace::memory_space>("mappings", mesh.nelems());
	auto mappings_h = Kokkos::create_mirror_view(mappings);
	for (int i = 0; i < mesh.nelems(); i++)
	{
		mappings_h[i] = Mapping2D(i, mesh);
	}
	Kokkos::deep_copy(execution_space, mappings, mappings_h);
}

Kokkos::View<PointSearch2D::Result*> TreePointSearch2D::apply(
	const CoordinateView<Omega_h::ExecSpace::memory_space>& coords) const
{
	if (coords.GetCoordinateSystem() != pcms::CoordinateSystem::Cartesian)
	{
		throw pcms_error("TreePointSearch2D::apply only implemented for"
						 " Cartesian coordinates");
	}
	
	Omega_h::ExecSpace execution_space;
	Kokkos::View<PointSearch2D::Result*> intersection_results("2D intersection results",
																coords.GetCoordinates().extent(0));
	auto results_h = Kokkos::create_mirror_view(intersection_results);
	for (int i = 0; i < intersection_results.size(); i++) 
		results_h[i] = Result_t{
			.dimensionality = Dim_t::REGION,
			.element_id = -1,
			.parametric_coords = {-1, -1, -1}
		};
	Kokkos::deep_copy(execution_space, intersection_results, results_h);

	auto CallOnIntersect = KOKKOS_LAMBDA <typename Predicate, typename Value>
	(Predicate const &predicate, Value const & val)
	{
		ArborX::Point<2, Omega_h::Real> const& ax = ArborX::getGeometry(predicate);
		Omega_h::Vector<2> point{ax[0], ax[1]};
		int point_ind = ArborX::getData(predicate);
		
		Mapping2D const& tm = mappings(val.index);
		
		// calculate the barycentric coefficients of the point
		auto coeffs = tm.get_bary(point);

		for (int i = 0; i < DIM; i++)
		{
			int elem = tm.which(i, coeffs);
			if (elem >= 0 && intersection_results(point_ind).dimensionality > (Dim_t)i)
			{
				auto face2elem = mesh_.get_adj(Omega_h::FACE, i).ab2b;
				auto elem_ind = face2elem[3*val.index + elem];
				intersection_results(point_ind) = PointSearch2D::Result{
					.dimensionality = (Dim_t)i, 
					.element_id = elem_ind, 
					.parametric_coords = coeffs
				};
				return;
			}
		}
		if (tm.which(DIM, coeffs) >= 0 && intersection_results(point_ind).dimensionality > (Dim_t)DIM)
		{
			intersection_results(point_ind) = PointSearch2D::Result{
				.dimensionality = Dim_t::FACE, 
				.element_id = (LO)val.index, 
				.parametric_coords = coeffs
			};
		}
	};

	tree.query(execution_space, 
			   pcms_Coordinate_View_Tagged<Omega_h::ExecSpace::memory_space, 2>{coords},
			   CallOnIntersect);
	return intersection_results;
}

} // namespace pcms


