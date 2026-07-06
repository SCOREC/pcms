#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

#include <pcms/localization/point_localization.h>

TEST_CASE ("Test 1x1 2D grid point classification") {
	// Setup for 1x1 grid test cases
	auto lib = Omega_h::Library{};
	auto world = lib.world();
	Omega_h::Mesh mesh_1x1 = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 1, 1, 0, false);
	REQUIRE(mesh_1x1.dim() == 2);
	Omega_h::ExecSpace execs;
	
	pcms::TreePointSearch2D tree_search(mesh_1x1);
	// End setup for 1x1 grid test cases

	auto check_verts = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nverts(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch2D::Dim_t::VERTEX);
			CHECK(results(i).element_id == i);
			CHECK_THAT(results(i).parametric_coords,
				Catch::Matchers::Contains(Catch::Matchers::WithinAbs(1.0, 10e-7)));
		}
	};

	auto check_edges = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch2D::Dim_t::EDGE);
			CHECK(results(i).element_id == i);
		}
	};

	SECTION ("Vertex intersection") {
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> vert_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(mesh_1x1.coords(), 2));
		Kokkos::View<pcms::TreePointSearch2D::Result_t*> results 
			= tree_search.apply(vert_cv);
		check_verts(results);
	}
	SECTION ("Edge intersection") {
		Omega_h::Write<pcms::Real> coords(mesh_1x1.nedges()*2);
		auto edge2vert = mesh_1x1.get_adj(Omega_h::EDGE, Omega_h::VERT).ab2b;
		auto meshcoords = mesh_1x1.coords();
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			Omega_h::Vector<2> v0{
				meshcoords[edge2vert[2*i]*2], meshcoords[edge2vert[2*i]*2 + 1]
			};
			Omega_h::Vector<2> v1{
				meshcoords[edge2vert[2*i + 1]*2], meshcoords[edge2vert[2*i + 1]*2 + 1]
			};
			Omega_h::Vector<2> midpoint = (v0 + v1)/2.;
			coords[i*2] = midpoint[0];
			coords[i*2 + 1] = midpoint[1];
		}

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> edge_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 2));
		Kokkos::View<pcms::TreePointSearch2D::Result_t*> results 
			= tree_search.apply(edge_cv);
		check_edges(results);
	}
	SECTION ("Face intersection") {
		
	}
}
