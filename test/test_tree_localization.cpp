#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

#include <pcms/localization/point_localization.h>

TEST_CASE ("Test 10x10 2D tree point classification") {
	// Setup for 10x10 tree test cases
	auto lib = Omega_h::Library{};
	auto world = lib.world();
	Omega_h::Mesh mesh_1x1 = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
	REQUIRE(mesh_1x1.dim() == 2);
	Omega_h::ExecSpace execs;
	
	pcms::TreePointSearch tree_search(mesh_1x1);
	// End setup for 10x10 tree test cases

	auto check_res = [] (auto const& results, int dim)
	{
		auto dims = Kokkos::create_mirror_view(results.dimensionalities);
		auto ids = Kokkos::create_mirror_view(results.element_ids);
		Kokkos::deep_copy(dims, results.dimensionalities);
		Kokkos::deep_copy(ids, results.element_ids);
		for (int i = 0; i < dims.size(); i++)
		{
			CHECK(dims(i) == (pcms::TreePointSearch::Dimensionality)dim);
			CHECK(ids(i) == i);
			// CHECK_THAT(results(i).parametric_coords,
			// 	Catch::Matchers::Contains(Catch::Matchers::WithinAbs(1.0, 10e-7)));
		}
	};

	SECTION ("Vertex intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("vertex coordinates", mesh_1x1.nverts(), 2);
		auto coords_h = Kokkos::create_mirror_view(coords);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		for (int i = 0; i < mesh_1x1.nverts(); i++)
		{
			coords_h(i,0) = meshcoords[i*2];
			coords_h(i,1) = meshcoords[i*2 + 1];
		}
		Kokkos::deep_copy(execs, coords, coords_h);
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> vert_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(vert_cv);
		check_res(results, 0);
	}
	SECTION ("Edge intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("edge coordinates", mesh_1x1.nedges(), 2);
		auto coords_h = Kokkos::create_mirror_view(coords);
		auto edge2vert = Omega_h::HostRead(mesh_1x1.get_adj(Omega_h::EDGE, Omega_h::VERT).ab2b);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			Omega_h::Vector<2> v0{
				meshcoords[edge2vert[2*i]*2], meshcoords[edge2vert[2*i]*2 + 1]
			};
			Omega_h::Vector<2> v1{
				meshcoords[edge2vert[2*i + 1]*2], meshcoords[edge2vert[2*i + 1]*2 + 1]
			};
			Omega_h::Vector<2> midpoint = (v0 + v1)/2.;
			coords_h(i,0) = midpoint[0];
			coords_h(i,1) = midpoint[1];
		}
		Kokkos::deep_copy(execs, coords, coords_h);

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> edge_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(edge_cv);
		check_res(results, 1);
	}
	SECTION ("Face intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("face coordinates", mesh_1x1.nfaces(), 2);
		auto face2vert = Omega_h::HostRead(mesh_1x1.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		auto coords_h = Kokkos::create_mirror_view(coords);
		for (int i = 0; i < mesh_1x1.nfaces(); i++)
		{
			Omega_h::Vector<2> v0{
				meshcoords[face2vert[3*i]*2], meshcoords[face2vert[3*i]*2 + 1]
			};
			Omega_h::Vector<2> v1{
				meshcoords[face2vert[3*i + 1]*2], meshcoords[face2vert[3*i + 1]*2 + 1]
			};
			Omega_h::Vector<2> v2{
				meshcoords[face2vert[3*i + 2]*2], meshcoords[face2vert[3*i + 2]*2 + 1]
			};
			Omega_h::Vector<2> midpoint = (v0 + v1 + v2)/3.;
			coords_h(i,0) = midpoint[0];
			coords_h(i,1) = midpoint[1];
		}
		Kokkos::deep_copy(execs, coords, coords_h);

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> face_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(face_cv);
		check_res(results, 2);
	}
	// SECTION("Outside Mesh") {
	// 	Omega_h::Write<pcms::Real> coords
	// 	{
	// 		-0.5, -0.5,
	// 		-0.5, 0.0,
	// 		-0.5, 0.5,
	// 		-0.5, 1.0,
	// 		-0.5, 1.5,
	// 		0.0, -0.5,
	// 		0.0, 1.5,
	// 		0.5, -0.5,
	// 		0.5, 1.5,
	// 		1.0, -0.5,
	// 		1.0, 1.5,
	// 		1.5, -0.5,
	// 		1.5, 0.0,
	// 		1.5, 0.5,
	// 		1.5, 1.0,
	// 		1.5, 1.5
	// 	};
	// 	pcms::CoordinateView<pcms::PointSearch::MemorySpace> outside_cv(
	// 		pcms::CoordinateSystem::Cartesian,
	// 		pcms::MakeConstRank2View<pcms::Real>(Omega_h::Read<pcms::Real>(coords), 2));
	// 	pcms::TreePointSearch::Results results = tree_search.apply(outside_cv);
	// 	check_res(results, 4);
	// }
}

TEST_CASE ("Test 5x5x5 3D tree point classification") {
	// Setup for 5x5x5 tree test cases
	auto lib = Omega_h::Library{};
	auto world = lib.world();
	Omega_h::Mesh mesh_1x1 = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 5, 5, false);
	REQUIRE(mesh_1x1.dim() == 3);
	Omega_h::ExecSpace execs;
	
	pcms::TreePointSearch tree_search(mesh_1x1);

	auto check_res = [] (auto const& results, int dim)
	{
		auto dims = Kokkos::create_mirror_view(results.dimensionalities);
		auto ids = Kokkos::create_mirror_view(results.element_ids);
		Kokkos::deep_copy(dims, results.dimensionalities);
		Kokkos::deep_copy(ids, results.element_ids);
		for (int i = 0; i < dims.size(); i++)
		{
			CHECK(dims(i) == (pcms::TreePointSearch::Dimensionality)dim);
			CHECK(ids(i) == i);
			// CHECK_THAT(results(i).parametric_coords,
			// 	Catch::Matchers::Contains(Catch::Matchers::WithinAbs(1.0, 10e-7)));
		}
	};
	
	SECTION ("Vertex intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("vertex coordinates", mesh_1x1.nverts(), 3);
		auto coords_h = Kokkos::create_mirror_view(coords);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		for (int i = 0; i < mesh_1x1.nverts(); i++)
		{
			coords_h(i, 0) = meshcoords[i*3];
			coords_h(i, 1) = meshcoords[i*3 + 1];
			coords_h(i, 2) = meshcoords[i*3 + 2];
		}
		Kokkos::deep_copy(execs, coords, coords_h);
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> vert_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(vert_cv);
		
		check_res(results, 0);
	}
	SECTION ("Edge intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("intersection coords", mesh_1x1.nedges(), 3);
		auto edge2vert = Omega_h::HostRead(mesh_1x1.get_adj(Omega_h::EDGE, Omega_h::VERT).ab2b);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		auto coords_h = Kokkos::create_mirror_view(coords);

		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			Omega_h::Vector<3> v0{
				meshcoords[edge2vert[2*i]*3], meshcoords[edge2vert[2*i]*3 + 1], meshcoords[edge2vert[2*i]*3 + 2]
			};
			Omega_h::Vector<3> v1{
				meshcoords[edge2vert[2*i + 1]*3], meshcoords[edge2vert[2*i + 1]*3 + 1], meshcoords[edge2vert[2*i + 1]*3 + 2]
			};
			Omega_h::Vector<3> midpoint = (v0 + v1)/2.;
			coords_h(i, 0) = midpoint[0];
			coords_h(i, 1) = midpoint[1];
			coords_h(i, 2) = midpoint[2];
		}

		Kokkos::deep_copy(execs, coords, coords_h);

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> edge_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(edge_cv);
		check_res(results, 1);
	}
	SECTION ("Face intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("intersection coords", mesh_1x1.nfaces(), 3);
		auto face2vert = Omega_h::HostRead(mesh_1x1.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		auto coords_h = Kokkos::create_mirror_view(coords);

		for (int i = 0; i < mesh_1x1.nfaces(); i++)
		{
			Omega_h::Vector<3> v0{
				meshcoords[face2vert[3*i]*3], meshcoords[face2vert[3*i]*3 + 1], meshcoords[face2vert[3*i]*3 + 2]
			};
			Omega_h::Vector<3> v1{
				meshcoords[face2vert[3*i + 1]*3], meshcoords[face2vert[3*i + 1]*3 + 1], meshcoords[face2vert[3*i + 1]*3 + 2]
			};
			Omega_h::Vector<3> v2{
				meshcoords[face2vert[3*i + 2]*3], meshcoords[face2vert[3*i + 2]*3 + 1], meshcoords[face2vert[3*i + 2]*3 + 2]
			};
			Omega_h::Vector<3> midpoint = (v0 + v1 + v2)/3.;
			coords_h(i, 0) = midpoint[0];
			coords_h(i, 1) = midpoint[1];
			coords_h(i, 2) = midpoint[2];
		}
		Kokkos::deep_copy(execs, coords, coords_h);

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> face_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(face_cv);
		check_res(results, 2);
	}
	SECTION ("Region intersection") {
		Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords("intersection coords", mesh_1x1.nelems(), 3);
		auto region2vert = Omega_h::HostRead(mesh_1x1.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b);
		auto meshcoords = Omega_h::HostRead(mesh_1x1.coords());
		auto coords_h = Kokkos::create_mirror_view(coords);

		for (int i = 0; i < mesh_1x1.nelems(); i++)
		{
			Omega_h::Vector<3> v0{
				meshcoords[region2vert[4*i]*3], meshcoords[region2vert[4*i]*3 + 1], meshcoords[region2vert[4*i]*3 + 2]
			};
			Omega_h::Vector<3> v1{
				meshcoords[region2vert[4*i + 1]*3], meshcoords[region2vert[4*i + 1]*3 + 1], meshcoords[region2vert[4*i + 1]*3 + 2]
			};
			Omega_h::Vector<3> v2{
				meshcoords[region2vert[4*i + 2]*3], meshcoords[region2vert[4*i + 2]*3 + 1], meshcoords[region2vert[4*i + 2]*3 + 2]
			};
			Omega_h::Vector<3> v3{
				meshcoords[region2vert[4*i + 3]*3], meshcoords[region2vert[4*i + 3]*3 + 1], meshcoords[region2vert[4*i + 3]*3 + 2]
			};
			Omega_h::Vector<3> midpoint = (v0 + v1 + v2 + v3)/4.;
			coords_h(i, 0) = midpoint[0];
			coords_h(i, 1) = midpoint[1];
			coords_h(i, 2) = midpoint[2];
		}
		Kokkos::deep_copy(execs, coords, coords_h);

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> region_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(coords));
		pcms::TreePointSearch::Results results = tree_search.apply(region_cv);
		check_res(results, 3);
	}

	// SECTION ("Region intersection") {
	// 	Kokkos::View<pcms::Real**, pcms::TreePointSearch::MemorySpace> coords
	// 	{
	// 		-0.5, -0.5, -0.5,
	// 		-0.5, -0.5, 1.5,
	// 		-0.5, 1.5, -0.5,
	// 		-0.5, 1.5, 1.5,
	// 		1.5, -0.5, -0.5,
	// 		1.5, -0.5, 1.5,
	// 		1.5, 1.5, -0.5,
	// 		1.5, 1.5, 1.5
	// 	};

	// 	pcms::CoordinateView<Omega_h::ExecSpace::memory_space, pcms::detail::default_layout_for_memory_space_t<Omega_h::ExecSpace::memory_space>> outside_cv(
	// 		pcms::CoordinateSystem::Cartesian,
	// 		pcms::MakeConstRank2View(coords));
	// 	Kokkos::View<pcms::TreePointSearch::Result*> results 
	// 		= tree_search.apply(outside_cv);
// auto results_h = Kokkos::create_mirror_view(results);
// Kokkos::deep_copy(execs, results_h, results);
	// 	check_outside(results_h);
	// }
}
