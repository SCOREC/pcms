#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/MLSCoefficients.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

TEST_CASE("basis test")
{
    SECTION("check basis slice lengths, size of the basis vector and basis vector for degree = 3 and dim = 3")
    {

	const int dim = 3;
	const int degree = 3;
	Kokkos::View<int**,Kokkos::HostSpace> array("array", degree, dim);
	Kokkos::deep_copy(array,0);
	basisSliceLengths(array);
	auto size = basisSize(array);  

	int expected[degree][dim] = {
	    {1, 1, 1},
	    {1, 2, 3},
	    {1, 3, 6}
	    };

	for (int i = 0; i < degree; ++i) {
	    for (int j = 0; j < dim; ++j) {
	    REQUIRE(array(i, j) == expected[i][j]);
	    }
	}

	REQUIRE(size == 20);


	MatViewType d_array("array in device", degree, dim);
	auto array_hd = Kokkos::create_mirror_view(d_array);
	Kokkos::deep_copy(array_hd,array);
	Kokkos::deep_copy(d_array, array_hd);

	int nvertices_target = 2;
	Omega_h::Matrix<3, 2> coords{{2.0,3.0,4.0}, {0.4, 0.5, 0.2}};
	Kokkos::View<double**> results("stores results", nvertices_target, size);

	auto host_results = Kokkos::create_mirror_view(results);

	team_policy tp(nvertices_target, Kokkos::AUTO);
	Kokkos::parallel_for("inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
	KOKKOS_LAMBDA(const member_type& team){
	    int i = team.league_rank();
	    ScratchVecView basis_vector(team.team_scratch(0), size);
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,size),[=](int j){
	    basis_vector(j) = 0;
	    });
	
	    Coord target_point;
	    target_point.x = coords[i][0];
	    target_point.y = coords[i][1];
	    target_point.z = coords[i][2];
	    
	    BasisPoly(basis_vector, d_array, target_point);
	    
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,size),[=](int j){
	    results(i,j) = basis_vector(j);
	    });


	});

	Kokkos::deep_copy(host_results, results);
	Kokkos::View<double*, Kokkos::HostSpace> expected_basis("expected vector", size);
	std::vector<std::string> names = {"1", "x", "y", "z", "x^2", "xy", "y^2", "xz", "yz", "z^2", "x^3", "x^2y", "xy^2", "Y^3", "x^2z", "xyz", "y^2z", "xz^2", "yz^2", "z^3"};

	for (int i = 0; i < nvertices_target; ++i){
	    double x = coords[i][0];
	    double y = coords[i][1];
	    double z = coords[i][2];
	    expected_basis(0)  = 1;
	    expected_basis(1)  = x;         // (x)
	    expected_basis(2)  = y;         // (y)
	    expected_basis(3)  = z;         // (z)
	    expected_basis(4)  = x * x;   // (x^2)
	    expected_basis(5)  = x * y;   // (xy)
	    expected_basis(6)  = y * y;   // (y^2)
	    expected_basis(7)  = x * z;   // (xz)
	    expected_basis(8)  = y * z;   // (yz)
	    expected_basis(9)  = z * z;   // (z^2)
	    expected_basis(10) = x * x * x;   // (x^3)
	    expected_basis(11) = x * x * y;   // (x^2y)
	    expected_basis(12) = x * y * y;   // (xy^2)
	    expected_basis(13) = y * y * y;   // (y^3)
	    expected_basis(14) = x * x * z;   // (x^2z)
	    expected_basis(15) = x * y * z;   // (xyz)
	    expected_basis(16) = y * y * z;   // (y^2z)
	    expected_basis(17) = x * z * z;   // (xz^2)
	    expected_basis(18) = y * z * z;   // (yz^2)
	    expected_basis(19) = z * z * z;   // (z^3)
	    for (int j = 0; j < size; ++j){
	    std::cout<<names[j]<<" "<<expected_basis(j)<<" "<<host_results(i,j)<<"\n";
	    REQUIRE(expected_basis(j) == host_results(i,j));
	    }
	}

    
}


    SECTION("check basis slice lengths, size of the basis vector and basis vector for degree = 1 and dim = 3")
    {

	const int dim = 3;
	const int degree = 1;
	Kokkos::View<int**,Kokkos::HostSpace> array("array", degree, dim);
	Kokkos::deep_copy(array,0);
	basisSliceLengths(array);

	auto size = basisSize(array);  

	int expected[degree][dim] = {
		{1, 1, 1},
	    };

	for (int i = 0; i < degree; ++i) {
	    for (int j = 0; j < dim; ++j) {
		REQUIRE(array(i, j) == expected[i][j]);
	    }
	}

	REQUIRE(size == 4);


	MatViewType d_array("array in device", degree, dim);
	auto array_hd = Kokkos::create_mirror_view(d_array);
	Kokkos::deep_copy(array_hd,array);
	Kokkos::deep_copy(d_array, array_hd);
    
	int nvertices_target = 2;
	Omega_h::Matrix<3, 2> coords{{2.0,3.0,4.0}, {0.4, 0.5, 0.2}};
	Kokkos::View<double**> results("stores results", nvertices_target, size);

	auto host_results = Kokkos::create_mirror_view(results);

	team_policy tp(nvertices_target, Kokkos::AUTO);
	Kokkos::parallel_for("inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
	KOKKOS_LAMBDA(const member_type& team){
	    int i = team.league_rank();
	    ScratchVecView basis_vector(team.team_scratch(0), size);
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,size),[=](int j){
		basis_vector(j) = 0;
	    });
	
	    Coord target_point;
	    target_point.x = coords[i][0];
	    target_point.y = coords[i][1];
	    target_point.z = coords[i][2];
	    
	    BasisPoly(basis_vector, d_array, target_point);
	    
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,size),[=](int j){
		results(i,j) = basis_vector(j);
	    });


	});

	Kokkos::deep_copy(host_results, results);
	Kokkos::View<double*, Kokkos::HostSpace> expected_basis("expected vector", size);
	std::vector<std::string> names = {"1", "x", "y", "z"};

	for (int i = 0; i < nvertices_target; ++i){
	    double x = coords[i][0];
	    double y = coords[i][1];
	    double z = coords[i][2];
	    expected_basis(0)  = 1;
	    expected_basis(1)  = x;         // (x)
	    expected_basis(2)  = y;         // (y)
	    expected_basis(3)  = z;         // (z)
	    for (int j = 0; j < size; ++j){
		std::cout<<names[j]<<" "<<expected_basis(j)<<" "<<host_results(i,j)<<"\n";
 		REQUIRE(expected_basis(j) == host_results(i,j));
	    }
	}

    }
  

}

