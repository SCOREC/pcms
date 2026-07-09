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
	
	pcms::TreePointSearch tree_search(mesh_1x1);
	// End setup for 1x1 grid test cases

	auto check_verts = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nverts(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::VERTEX);
			CHECK(results(i).element_id == i);
			// CHECK_THAT(results(i).parametric_coords,
			// 	Catch::Matchers::Contains(Catch::Matchers::WithinAbs(1.0, 10e-7)));
		}
	};

	auto check_edges = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::EDGE);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_faces = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < results.size(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::FACE);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_outside = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < results.size(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::REGION);
			CHECK(results(i).element_id == -1);
		}
	};

	SECTION ("Vertex intersection") {
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> vert_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(mesh_1x1.coords(), 2));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
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
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(edge_cv);
		check_edges(results);
	}
	SECTION ("Face intersection") {
		Omega_h::Write<pcms::Real> coords(mesh_1x1.nfaces()*2);
		auto face2vert = mesh_1x1.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b;
		auto meshcoords = mesh_1x1.coords();

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
			coords[i*2] = midpoint[0];
			coords[i*2 + 1] = midpoint[1];
		}

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> face_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 2));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(face_cv);
		check_faces(results);
	}
	SECTION("Outside Mesh") {
		Omega_h::Write<pcms::Real> coords
		{
			-0.5, -0.5,
			-0.5, 0.0,
			-0.5, 0.5,
			-0.5, 1.0,
			-0.5, 1.5,
			0.0, -0.5,
			0.0, 1.5,
			0.5, -0.5,
			0.5, 1.5,
			1.0, -0.5,
			1.0, 1.5,
			1.5, -0.5,
			1.5, 0.0,
			1.5, 0.5,
			1.5, 1.0,
			1.5, 1.5
		};
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> outside_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 2));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(outside_cv);
		check_outside(results);
	}
}

TEST_CASE ("Test 1x1x1 3D grid point classification") {
	// Setup for 1x1x1 grid test cases
	auto lib = Omega_h::Library{};
	auto world = lib.world();
	Omega_h::Mesh mesh_1x1 = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 1, 1, 1, false);
	REQUIRE(mesh_1x1.dim() == 3);
	Omega_h::ExecSpace execs;
	
	pcms::TreePointSearch tree_search(mesh_1x1);
	// End setup for 1x1 grid test cases

	auto check_verts = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nverts(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::VERTEX);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_edges = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::EDGE);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_faces = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < results.size(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::FACE);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_regions = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < results.size(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::REGION);
			CHECK(results(i).element_id == i);
		}
	};

	auto check_outside = KOKKOS_LAMBDA (auto const& results)
	{
		for (int i = 0; i < results.size(); i++)
		{
			CHECK(results(i).dimensionality == pcms::TreePointSearch::Dimensionality::REGION);
			CHECK(results(i).element_id == -1);
		}
	};

	SECTION ("Vertex intersection") {
		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> vert_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(mesh_1x1.coords(), 3));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(vert_cv);
		check_verts(results);
	}
	SECTION ("Edge intersection") {
		Omega_h::Write<pcms::Real> coords(mesh_1x1.nedges()*3);
		auto edge2vert = mesh_1x1.get_adj(Omega_h::EDGE, Omega_h::VERT).ab2b;
		auto meshcoords = mesh_1x1.coords();
		for (int i = 0; i < mesh_1x1.nedges(); i++)
		{
			Omega_h::Vector<3> v0{
				meshcoords[edge2vert[2*i]*3], meshcoords[edge2vert[2*i]*3 + 1], meshcoords[edge2vert[2*i]*3 + 2]
			};
			Omega_h::Vector<3> v1{
				meshcoords[edge2vert[2*i + 1]*3], meshcoords[edge2vert[2*i + 1]*3 + 1], meshcoords[edge2vert[2*i + 1]*3 + 2]
			};
			Omega_h::Vector<3> midpoint = (v0 + v1)/2.;
			coords[i*3] = midpoint[0];
			coords[i*3 + 1] = midpoint[1];
			coords[i*3 + 2] = midpoint[2];
		}

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> edge_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 3));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(edge_cv);
		check_edges(results);
	}
	SECTION ("Face intersection") {
		Omega_h::Write<pcms::Real> coords(mesh_1x1.nfaces()*3);
		auto face2vert = mesh_1x1.get_adj(Omega_h::FACE, Omega_h::VERT).ab2b;
		auto meshcoords = mesh_1x1.coords();

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
			coords[i*3] = midpoint[0];
			coords[i*3 + 1] = midpoint[1];
			coords[i*3 + 2] = midpoint[2];
		}

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> face_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 3));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(face_cv);
		check_faces(results);
	}
	SECTION ("Region intersection") {
		Omega_h::Write<pcms::Real> coords(mesh_1x1.nelems()*3);
		auto region2vert = mesh_1x1.get_adj(Omega_h::REGION, Omega_h::VERT).ab2b;
		auto meshcoords = mesh_1x1.coords();

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
			coords[i*3] = midpoint[0];
			coords[i*3 + 1] = midpoint[1];
			coords[i*3 + 2] = midpoint[2];
		}

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> region_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 3));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(region_cv);
		check_regions(results);
	}

	SECTION ("Region intersection") {
		Omega_h::Write<pcms::Real> coords
		{
			-0.5, -0.5, -0.5,
			-0.5, -0.5, 1.5,
			-0.5, 1.5, -0.5,
			-0.5, 1.5, 1.5,
			1.5, -0.5, -0.5,
			1.5, -0.5, 1.5,
			1.5, 1.5, -0.5,
			1.5, 1.5, 1.5
		};

		pcms::CoordinateView<Omega_h::ExecSpace::memory_space> outside_cv(
			pcms::CoordinateSystem::Cartesian,
			pcms::MakeConstRank2View(Omega_h::Read<pcms::Real>(coords), 3));
		Kokkos::View<pcms::TreePointSearch::Result*> results 
			= tree_search.apply(outside_cv);
		check_outside(results);
	}
}
