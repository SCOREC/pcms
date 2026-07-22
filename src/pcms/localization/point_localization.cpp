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

Mapping<2>::Mapping(Omega_h::Matrix<Mapping<2>::DIM,Mapping<2>::DIM+1> const& triangle_) : triangle(triangle_)
{
	bary_transform = { triangle[0] - triangle[2], triangle[1] - triangle[2] };
	triangle_area = 0.5*fabs(Omega_h::determinant(bary_transform));
	bary_transform = Omega_h::invert(bary_transform);
}

KOKKOS_FUNCTION
Omega_h::Vector<Mapping<2>::DIM + 1> Mapping<2>::get_bary(Omega_h::Vector<Mapping<2>::DIM> const& p) const
{
	Omega_h::Vector<2> coeffs = bary_transform*(p - triangle[2]);
	return {coeffs[0], coeffs[1], 1 - coeffs[0] - coeffs[1]};
}

KOKKOS_FUNCTION
int Mapping<2>::which(int dim, 
					 Omega_h::Vector<Mapping<2>::DIM + 1> const& bary_coords) const
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

int Mapping<2>::which_vert(Omega_h::Vector<Mapping<2>::DIM + 1> const& bary_coords) const
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

int Mapping<2>::which_edge(Omega_h::Vector<Mapping<2>::DIM + 1> const& bary_coords) const
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

int Mapping<2>::within_elem(Omega_h::Vector<Mapping<2>::DIM + 1> const& bary_coords) const
{
	return -1 * (int)!(bary_coords[0] > 0 && bary_coords[1] > 0 && bary_coords[2] > 0);
}

double Mapping<2>::opposite_edge_len_sq(int i) const
{
	return Omega_h::norm_squared(triangle[(i+2)%3] - triangle[(i+1)%3]);
}

Mapping<3>::Mapping(Omega_h::Matrix<Mapping<3>::DIM,Mapping<3>::DIM+1> const& tetrahedron_)
	: tetrahedron(tetrahedron_)
{
	set_triangle_areas();
	bary_transform = { tetrahedron[0] - tetrahedron[3], tetrahedron[1] - tetrahedron[3], tetrahedron[2] - tetrahedron[3] };
	tetrahedron_volume = fabs(Omega_h::determinant(bary_transform))/6.;
	bary_transform = Omega_h::invert(bary_transform);
}

KOKKOS_FUNCTION
Omega_h::Vector<Mapping<3>::DIM + 1> Mapping<3>::get_bary(Omega_h::Vector<Mapping<3>::DIM> const& p) const
{
	Omega_h::Vector<3> coeffs = bary_transform*(p - tetrahedron[3]);
	return {coeffs[0], coeffs[1], coeffs[2], 1 - coeffs[0] - coeffs[1] - coeffs[2]};
}

KOKKOS_FUNCTION
int Mapping<3>::which(int dim, 
					 Omega_h::Vector<Mapping<3>::DIM + 1> const& bary_coords) const
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

int Mapping<3>::which_vert(Omega_h::Vector<Mapping<3>::DIM + 1> const& bary_coords) const
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


int Mapping<3>::which_edge(Omega_h::Vector<Mapping<3>::DIM + 1> const& bary_coords) const
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

int Mapping<3>::which_face(Omega_h::Vector<Mapping<3>::DIM + 1> const& bary_coords) const
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

int Mapping<3>::within_elem(Omega_h::Vector<Mapping<3>::DIM + 1> const& bary_coords) const
{
	return -1 * (int)!(bary_coords[0] >= 0 && bary_coords[1] >= 0 && bary_coords[2] >= 0 && bary_coords[3] >= 0);
}

// calculates the areas of each of the faces of the triangle
// the faces are "named" according to the vertex opposite
// (i.e., face 0 is defined by points 1, 2, 3)
void Mapping<3>::set_triangle_areas()
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
	
	detail::Mapping<2> const& tm = mappings(val.index);
	
	// calculate the barycentric coefficients of the point
	auto coeffs = tm.get_bary(point);

	for (int i = 0; i < DIM; i++)
	{
		int elem = tm.which(i, coeffs);
		if (elem >= 0 && dimensionalities(point_ind) > (TreePointSearch::Dimensionality)i)
		{
			auto elem_ind = adjacencies[i][3*val.index + elem];
			dimensionalities(point_ind) = (TreePointSearch::Dimensionality)i;
			element_ids(point_ind) = elem_ind;
			for (int j = 0; j < DIM + 1; j++)
			{
				parametric_coords(point_ind, j) = coeffs(j);
			}
			return;
		}
	}
	if (tm.which(DIM, coeffs) >= 0 
		&& dimensionalities(point_ind) > TreePointSearch::Dimensionality::FACE)
	{
		dimensionalities(point_ind) = TreePointSearch::Dimensionality::FACE;
		element_ids(point_ind) = (LO)val.index;
		for (int j = 0; j < DIM + 1; j++)
		{
			parametric_coords(point_ind, j) = coeffs(j);
		}
	}
}

template <typename Predicate, typename Value>
void CallOnIntersect3D::operator()(Predicate const &predicate, Value const & val) const
{
	ArborX::Point<DIM, Omega_h::Real> const& ax = ArborX::getGeometry(predicate);
	Omega_h::Vector<DIM> point{ax[0], ax[1], ax[2]};
	int point_ind = ArborX::getData(predicate);
	
	detail::Mapping<3> const& tm = mappings(val.index);
	
	// calculate the barycentric coefficients of the point
	auto coeffs = tm.get_bary(point);
	int offsets[DIM] = {4, 6, 4};
	for (int i = 0; i < DIM; i++)
	{
		int elem = tm.which(i, coeffs);
		if (elem >= 0 && dimensionalities(point_ind) > (TreePointSearch::Dimensionality)i)
		{
			auto elem_ind = adjacencies[i][offsets[i]*val.index + elem];
			dimensionalities(point_ind) = (TreePointSearch::Dimensionality)i;
			element_ids(point_ind) = elem_ind;
			for (int j = 0; j < DIM + 1; j++)
			{
				parametric_coords(point_ind, j) = coeffs(j);
			}
			return;
		}
	}

	if (tm.which(DIM, coeffs) >= 0 
			&& dimensionalities(point_ind) > TreePointSearch::Dimensionality::REGION)
	{
		dimensionalities(point_ind) = TreePointSearch::Dimensionality::REGION;
		element_ids(point_ind) = (LO)val.index;
		for (int j = 0; j < DIM + 1; j++)
		{
			parametric_coords(point_ind, j) = coeffs(j);
		}
	}
}

} //namespace detail

TreePointSearch::Results TreePointSearch::apply(
	const CoordinateView<TreePointSearch::MemorySpace>& coords) const
{
	if (coords.GetCoordinateSystem() != pcms::CoordinateSystem::Cartesian)
	{
		throw pcms_error("TreePointSearch::apply only implemented for"
						 " Cartesian coordinates");
	}
	if (coords.GetValues().extent(1) != mesh_.dim())
	{
		throw pcms_error("Input coordinate space dimension " 
			+ std::to_string(coords.GetValues().extent(1))
			+ " does not match query space dimension " 
			+ std::to_string(mesh_.dim()));
	}
	
	if (mesh_.dim() == 2)
	{
		static constexpr int DIM = 2;
		Omega_h::ExecSpace execution_space;

		Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace> dims("dimensionalities", coords.GetValues().extent(0));
		Kokkos::deep_copy(dims, TreePointSearch::Dimensionality::NO_INTERSECT);

		
		Kokkos::View<LO*, MemorySpace> elem_ids("element IDs", coords.GetValues().extent(0));
		Kokkos::deep_copy(elem_ids, -1);

		Kokkos::View<Real**, MemorySpace> parametric_coords("parametric coordinates", coords.GetValues().extent(0), coords.GetValues().extent(1) + 1);
		Kokkos::deep_copy(parametric_coords, -1.0);


		tree->get_tree<2>()->query(
			execution_space, 
			detail::Coordinate_View_Adapt<MemorySpace, DIM>{coords.GetValues()},
			detail::CallOnIntersect2D(
				*tree->get_mappings<2>(),
				mesh_.get_adj(Omega_h::FACE, 0).ab2b,
				mesh_.get_adj(Omega_h::FACE, 1).ab2b,
				dims,
				elem_ids,
				parametric_coords
			));
		return Results{dims, elem_ids, parametric_coords};
	}
	else
	{
		static constexpr int DIM = 3;
		Omega_h::ExecSpace execution_space;

		Kokkos::View<TreePointSearch::Dimensionality*, MemorySpace> dims("dimensionalities", coords.GetValues().extent(0));
		Kokkos::deep_copy(dims, TreePointSearch::Dimensionality::NO_INTERSECT);


		Kokkos::View<LO*, MemorySpace> elem_ids("element IDs", coords.GetValues().extent(0));
		Kokkos::deep_copy(elem_ids, -1);

		Kokkos::View<Real**, MemorySpace> parametric_coords("parametric coordinates", coords.GetValues().extent(0), coords.GetValues().extent(1) + 1);
		Kokkos::deep_copy(parametric_coords, -1.0);

		tree->get_tree<3>()->query(
			execution_space, 
			detail::Coordinate_View_Adapt<MemorySpace, DIM>{coords.GetValues()},
			detail::CallOnIntersect3D(
				*tree->get_mappings<3>(),
				mesh_.get_adj(Omega_h::REGION, 0).ab2b,
				mesh_.get_adj(Omega_h::REGION, 1).ab2b,
				mesh_.get_adj(Omega_h::REGION, 2).ab2b,
				dims,
				elem_ids,
				parametric_coords
			));
		return Results{dims, elem_ids, parametric_coords};
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
		
		detail::TreeWrapper::Tree_t<2> tree = detail::TreeWrapper::Tree_t<2>(
			execution_space,
			ArborX::Experimental::attach_indices(tagged_mesh));
			
		detail::TreeWrapper::Mappings_t<2> mappings = detail::TreeWrapper::Mappings_t<2>(
			"mappings", 
			mesh.nelems());
		
		auto face2vert = Omega_h::HostRead(mesh.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b);
		auto vert_coords = Omega_h::HostRead(mesh.coords());
		
		auto mappings_h = Kokkos::create_mirror_view(mappings);
		for (int i = 0; i < mesh.nelems(); i++)
		{
			Omega_h::Matrix<2,3> triangle = {
				{vert_coords[face2vert[i*3]*2], vert_coords[face2vert[i*3]*2 + 1]},
				{vert_coords[face2vert[i*3 + 1]*2], vert_coords[face2vert[i*3 + 1]*2 + 1]},
				{vert_coords[face2vert[i*3 + 2]*2], vert_coords[face2vert[i*3 + 2]*2 + 1]}};
			mappings_h[i] = detail::Mapping<2>(triangle);
		}
		Kokkos::deep_copy(execution_space, mappings, mappings_h);
		return std::make_unique<detail::TreeWrapper>(mappings, tree);
	}
	if (mesh.dim() == 3)
	{
		detail::Omega_h_Mesh_Adapt<3> tagged_mesh{
			mesh.nelems(),
			mesh.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b,
			mesh.coords()
		};
		detail::TreeWrapper::Tree_t<3> tree = detail::TreeWrapper::Tree_t<3>(
			execution_space, 
			ArborX::Experimental::attach_indices(tagged_mesh));
	
		detail::TreeWrapper::Mappings_t<3> mappings = detail::TreeWrapper::Mappings_t<3>(
			"mappings", 
			mesh.nelems());
		
		auto region2vert = Omega_h::HostRead(mesh.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b);
		auto vert_coords = Omega_h::HostRead(mesh.coords());

		auto mappings_h = Kokkos::create_mirror_view(mappings);
		for (int i = 0; i < mesh.nelems(); i++)
		{
			Omega_h::Matrix<3,4> tetrahedron = {
				{vert_coords[region2vert[i*4]*3], vert_coords[region2vert[i*4]*3 + 1], vert_coords[region2vert[i*4]*3 + 2]},
				{vert_coords[region2vert[i*4 + 1]*3], vert_coords[region2vert[i*4 + 1]*3 + 1], vert_coords[region2vert[i*4 + 1]*3 + 2]},
				{vert_coords[region2vert[i*4 + 2]*3], vert_coords[region2vert[i*4 + 2]*3 + 1], vert_coords[region2vert[i*4 + 2]*3 + 2]},
				{vert_coords[region2vert[i*4 + 3]*3], vert_coords[region2vert[i*4 + 3]*3 + 1], vert_coords[region2vert[i*4 + 3]*3 + 2]}};
			mappings_h[i] = detail::Mapping<3>(tetrahedron);
		}
		Kokkos::deep_copy(execution_space, mappings, mappings_h);
		return std::make_unique<detail::TreeWrapper>(mappings, tree);
	}
	throw pcms_error("Invalid mesh dimension " + std::to_string(mesh.dim()) + ", TreePointSearch only implemented for 2D and 3D");
}

} // namespace pcms
