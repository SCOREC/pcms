#include "point_localization.h"

#define FACE_TOL 10e-8
#define EDGE_TOL 10e-6
#define VERT_TOL 10e-5

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
		const pcms::detail::Omega_h_Mesh_Adapt<dim>& mesh, 
		int i)
	{
		ArborX::Point<dim, Omega_h::Real> min = {INFINITY};
		ArborX::Point<dim, Omega_h::Real> max = {-INFINITY};
		for (int j = 0; j < dim + 1; ++j)
		{
			auto cell_vert_id = mesh.adjacency[(dim+1)*i+j];
			for (int k = 0; k < dim; ++k)
			{
				Omega_h::Real curr_coord = mesh.coordinates[cell_vert_id*dim+k];
				if (min[k] > curr_coord) min[k] = curr_coord;
				if (max[k] < curr_coord) max[k] = curr_coord;
			}
		}
		return ArborX::Box(min, max);
	}
};

template <typename MemorySpace, int dim>
struct ArborX::AccessTraits<pcms::detail::Coordinate_View_Adapt<MemorySpace, dim>>
{
	using memory_space = MemorySpace;
	static KOKKOS_FUNCTION int size(
		pcms::detail::Coordinate_View_Adapt<MemorySpace, dim> const &coords)
	{
		return coords.points.extent(0);
	}
	static KOKKOS_FUNCTION auto get(
		pcms::detail::Coordinate_View_Adapt<MemorySpace, dim> const &coords, 
		int i)
	{
		ArborX::Point<dim, double> ax_point;
		for (int j = 0; j < dim; j++) ax_point[j] = coords.points(i,j);
		return PredicateWithAttachment(intersects(ax_point), i);
	}
};

namespace pcms
{

namespace detail
{

Mapping2D::Mapping2D(int elem_index, Omega_h::Mesh const& mesh)
{
	set_mesh_triangle(elem_index, mesh);
	bary_transform = { triangle[0] - triangle[2], triangle[1] - triangle[2] };
	triangle_area = 0.5*fabs(Omega_h::determinant(bary_transform));
	bary_transform = Omega_h::invert(bary_transform);
}

KOKKOS_FUNCTION
Omega_h::Vector<Mapping2D::DIM + 1> Mapping2D::get_bary(Omega_h::Vector<Mapping2D::DIM> const& p) const
{
	Omega_h::Vector<2> coeffs = bary_transform*(p - triangle[2]);
	return {coeffs[0], coeffs[1], 1 - coeffs[0] - coeffs[1]};
}

KOKKOS_FUNCTION
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

int Mapping2D::within_elem(Omega_h::Vector<Mapping2D::DIM + 1> const& bary_coords) const
{
	return -1 * (int)!(bary_coords[0] > 0 && bary_coords[1] > 0 && bary_coords[2] > 0);
}

void Mapping2D::set_mesh_triangle(int index, Omega_h::Mesh const& mesh)
{
	auto face2vert = Omega_h::HostRead(mesh.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b);
	auto vert_coords = Omega_h::HostRead(mesh.coords());
	triangle = Omega_h::Matrix<2,3>{{vert_coords[face2vert[index*3]*2], vert_coords[face2vert[index*3]*2 + 1]},
			{vert_coords[face2vert[index*3 + 1]*2], vert_coords[face2vert[index*3 + 1]*2 + 1]},
			{vert_coords[face2vert[index*3 + 2]*2], vert_coords[face2vert[index*3 + 2]*2 + 1]}};
}

double Mapping2D::opposite_edge_len_sq(int i) const
{
	return Omega_h::norm_squared(triangle[(i+2)%3] - triangle[(i+1)%3]);
}

Mapping3D::Mapping3D(int elem_index, Omega_h::Mesh const& mesh)
{
	set_mesh_tet(elem_index, mesh);
	set_triangle_areas();
	bary_transform = { tetrahedron[0] - tetrahedron[3], tetrahedron[1] - tetrahedron[3], tetrahedron[2] - tetrahedron[3] };
	tetrahedron_volume = fabs(Omega_h::determinant(bary_transform))/6.;
	bary_transform = Omega_h::invert(bary_transform);
}

KOKKOS_FUNCTION
Omega_h::Vector<Mapping3D::DIM + 1> Mapping3D::get_bary(Omega_h::Vector<Mapping3D::DIM> const& p) const
{
	Omega_h::Vector<3> coeffs = bary_transform*(p - tetrahedron[3]);
	return {coeffs[0], coeffs[1], coeffs[2], 1 - coeffs[0] - coeffs[1] - coeffs[2]};
}

KOKKOS_FUNCTION
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

int Mapping3D::within_elem(Omega_h::Vector<Mapping3D::DIM + 1> const& bary_coords) const
{
	return -1 * (int)!(bary_coords[0] >= 0 && bary_coords[1] >= 0 && bary_coords[2] >= 0 && bary_coords[3] >= 0);
}

void Mapping3D::set_mesh_tet(int index, Omega_h::Mesh const& mesh)
{
	auto region2vert = Omega_h::HostRead(mesh.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b);
	auto vert_coords = Omega_h::HostRead(mesh.coords());
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

template <typename Predicate, typename Value>
KOKKOS_FUNCTION void CallOnIntersect2D::operator()(Predicate const &predicate, Value const & val) const
{
	ArborX::Point<DIM, Omega_h::Real> const& ax = ArborX::getGeometry(predicate);
	Omega_h::Vector<DIM> point{ax[0], ax[1]};
	int point_ind = ArborX::getData(predicate);
	
	detail::Mapping2D const& tm = mappings(val.index);
	
	// calculate the barycentric coefficients of the point
	auto coeffs = tm.get_bary(point);

	for (int i = 0; i < DIM; i++)
	{
		int elem = tm.which(i, coeffs);
		if (elem >= 0 && intersection_results(point_ind).dimensionality > (TreePointSearch::Dimensionality)i)
		{
			auto elem_ind = adjacencies[i][3*val.index + elem];
			intersection_results(point_ind).dimensionality = (TreePointSearch::Dimensionality)i;
			intersection_results(point_ind).element_id = elem_ind;
			for (int j = 0; j < DIM + 1; j++)
				intersection_results(point_ind).parametric_coords(j) = coeffs[j];
			return;
		}
	}
	if (tm.which(DIM, coeffs) >= 0 
		&& intersection_results(point_ind).dimensionality > TreePointSearch::Dimensionality::FACE)
	{
		intersection_results(point_ind).dimensionality = TreePointSearch::Dimensionality::FACE;
		intersection_results(point_ind).element_id = (LO)val.index;
		for (int j = 0; j < DIM + 1; j++)
				intersection_results(point_ind).parametric_coords(j) = coeffs[j];
	}
}

template <typename Predicate, typename Value>
void CallOnIntersect3D::operator()(Predicate const &predicate, Value const & val) const
{
	ArborX::Point<DIM, Omega_h::Real> const& ax = ArborX::getGeometry(predicate);
	Omega_h::Vector<DIM> point{ax[0], ax[1], ax[2]};
	int point_ind = ArborX::getData(predicate);
	
	detail::Mapping3D const& tm = mappings(val.index);
	
	// calculate the barycentric coefficients of the point
	auto coeffs = tm.get_bary(point);
	int offsets[DIM] = {4, 6, 4};
	for (int i = 0; i < DIM; i++)
	{
		int elem = tm.which(i, coeffs);
		if (elem >= 0 && intersection_results(point_ind).dimensionality > (TreePointSearch::Dimensionality)i)
		{
			auto elem_ind = adjacencies[i][offsets[i]*val.index + elem];
			intersection_results(point_ind).dimensionality = (TreePointSearch::Dimensionality)i;
			intersection_results(point_ind).element_id = elem_ind;
			for (int j = 0; j < DIM + 1; j++)
				intersection_results(point_ind).parametric_coords(j) = coeffs[j];
			return;
		}
	}

	if (tm.which(DIM, coeffs) >= 0 
			&& intersection_results(point_ind).dimensionality >= TreePointSearch::Dimensionality::REGION)
	{
		intersection_results(point_ind).dimensionality = TreePointSearch::Dimensionality::REGION;
		intersection_results(point_ind).element_id = (LO)val.index;
		for (int j = 0; j < DIM + 1; j++)
				intersection_results(point_ind).parametric_coords(j) = coeffs[j];
	}
}

} //namespace detail

Kokkos::View<TreePointSearch::Result*> TreePointSearch::apply(
	const CoordinateView<TreePointSearch::MemorySpace>& coords) const
{
	if (coords.GetCoordinateSystem() != pcms::CoordinateSystem::Cartesian)
	{
		throw pcms_error("TreePointSearch::apply only implemented for"
						 " Cartesian coordinates");
	}
	if (coords.GetCoordinates().extent(1) != mesh_.dim())
	{
		throw pcms_error("Input coordinate space dimension " 
			+ std::to_string(coords.GetCoordinates().extent(1))
			+ " does not match query space dimension " 
			+ std::to_string(mesh_.dim()));
	}
	
	if (mesh_.dim() == 2)
	{
		static constexpr int DIM = 2;
		Omega_h::ExecSpace execution_space;
		Kokkos::View<PointSearch::Result*, MemorySpace> 
			intersection_results("2D intersection results",
				coords.GetCoordinates().extent(0));
		auto results_h = Kokkos::create_mirror_view(intersection_results);
		for (int i = 0; i < intersection_results.size(); i++) 
		{
			results_h[i].dimensionality = Dimensionality::REGION,
			results_h[i].element_id = -1,
			results_h[i].parametric_coords = {-1, -1, -1, -1};
		}
		Kokkos::deep_copy(execution_space, intersection_results, results_h);

		((detail::TreeWrapper2D::Tree_t*)tree->get_tree())->query(execution_space, 
					detail::Coordinate_View_Adapt<MemorySpace, DIM>{
					coords.GetCoordinates()},
					detail::CallOnIntersect2D(
						*(detail::TreeWrapper2D::Mappings_t*)tree->get_mappings(),
						mesh_.get_adj(Omega_h::FACE, 0).ab2b,
						mesh_.get_adj(Omega_h::FACE, 1).ab2b,
						intersection_results
					));
		return intersection_results;
	}
	else
	{
		static constexpr int DIM = 3;
		Omega_h::ExecSpace execution_space;
		Kokkos::View<PointSearch::Result*, MemorySpace> 
			intersection_results("3D intersection results",
				coords.GetCoordinates().extent(0));
		auto results_h = Kokkos::create_mirror_view(intersection_results);
		for (int i = 0; i < intersection_results.size(); i++) 
		{
			results_h[i].dimensionality = Dimensionality::REGION,
			results_h[i].element_id = -1,
			results_h[i].parametric_coords = {-1, -1, -1, -1};
		}
		Kokkos::deep_copy(execution_space, intersection_results, results_h);
		
		detail::TreeWrapper3D::Mappings_t mappings = *(detail::TreeWrapper3D::Mappings_t*)tree->get_mappings();
		Omega_h::LOs adjacencies[3] = {
			mesh_.get_adj(Omega_h::REGION, 0).ab2b,
			mesh_.get_adj(Omega_h::REGION, 1).ab2b,
			mesh_.get_adj(Omega_h::REGION, 2).ab2b
		};

		((detail::TreeWrapper3D::Tree_t*)tree->get_tree())->query(execution_space, 
			detail::Coordinate_View_Adapt<MemorySpace, 
			DIM>{
			coords.GetCoordinates()},
			detail::CallOnIntersect3D(
				*(detail::TreeWrapper3D::Mappings_t*)tree->get_mappings(),
				mesh_.get_adj(Omega_h::REGION, 0).ab2b,
				mesh_.get_adj(Omega_h::REGION, 1).ab2b,
				mesh_.get_adj(Omega_h::REGION, 2).ab2b,
				intersection_results	
			));
		return intersection_results;
	}
}

std::unique_ptr<detail::TreeWrapper> TreePointSearch::make_tree(const Omega_h::Mesh& mesh) const
{
	ExecSpace execution_space;
	using DeviceType = Kokkos::Device<Omega_h::ExecSpace, 
									  Omega_h::ExecSpace::memory_space>;
	
	if (mesh.dim() == 2)
	{
		detail::Omega_h_Mesh_Adapt<2> tagged_mesh{
			mesh.nelems(),
			mesh.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b,
			mesh.coords()
		};
		
		detail::TreeWrapper2D::Tree_t tree = detail::TreeWrapper2D::Tree_t(
			execution_space,
			ArborX::Experimental::attach_indices(tagged_mesh));
			
		detail::TreeWrapper2D::Mappings_t mappings = detail::TreeWrapper2D::Mappings_t(
			"mappings", 
			mesh.nelems());
		
		auto mappings_h = Kokkos::create_mirror_view(mappings);
		for (int i = 0; i < mesh.nelems(); i++)
		{
			mappings_h[i] = detail::Mapping2D(i, mesh);
		}
		Kokkos::deep_copy(execution_space, mappings, mappings_h);
		return std::make_unique<detail::TreeWrapper2D>(mappings, tree);
	}
	if (mesh.dim() == 3)
	{
		detail::Omega_h_Mesh_Adapt<3> tagged_mesh{
			mesh.nelems(),
			mesh.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b,
			mesh.coords()
		};
		detail::TreeWrapper3D::Tree_t tree = detail::TreeWrapper3D::Tree_t(
			execution_space, 
			ArborX::Experimental::attach_indices(tagged_mesh));
	
		detail::TreeWrapper3D::Mappings_t mappings = detail::TreeWrapper3D::Mappings_t(
			"mappings", 
			mesh.nelems());
		
		auto mappings_h = Kokkos::create_mirror_view(mappings);
		for (int i = 0; i < mesh.nelems(); i++)
		{
			mappings_h[i] = detail::Mapping3D(i, mesh);
		}
		Kokkos::deep_copy(execution_space, mappings, mappings_h);
		return std::make_unique<detail::TreeWrapper3D>(mappings, tree);
	}
	throw pcms_error("Invalid mesh dimension " + std::to_string(mesh.dim()) + ", TreePointSearch only implemented for 2D and 3D");
}

} // namespace pcms
