#include <catch2/catch_test_macros.hpp>
#include <pcms/gnet/gnet_grid.hpp>
#include <pcms/gnet/gnet_io.hpp>
#include <pcms/field/evaluator/mls_interpolation_impl.hpp>
#include <cstdlib>
#include <string>

using namespace pcms::gnet;

namespace
{

// get the grid ids for the given radial bins, poloidal bins and toroidal bins
Omega_h::LOs get_grid_ids(const Omega_h::LOs& radial_bins, const Omega_h::LOs& poloidal_bins, const Omega_h::LOs& toroidal_bins, const size_t& mrho_max, const size_t& mtor_max){
			
		Omega_h::Write<Omega_h::LO> d_results_ids(4, 1, "Result Ids");

		auto pol_max_scan_array = get_pol_max_scan_array(mrho_max);	
		auto mpol_max_array = get_mpol_max_array(mrho_max);
		
		
		Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {
			int mpol_sum_at_current_radial = pol_max_scan_array[radial_bins[i]];
			auto grid_id = eval_grid_id(mpol_sum_at_current_radial, poloidal_bins[i], toroidal_bins[i], mtor_max); 	
			d_results_ids[i] = grid_id;	

		  });

		return Omega_h::read(d_results_ids);
}

void  test_radial_bin_neighbors(const size_t num_rho_neighbors, 
								const Omega_h::Reals& test_rho_coordinates, const int mrho_max, 
								Omega_h::Write<Omega_h::LO>& neighbors_ids) {

			size_t total_neighbors = 4 * num_rho_neighbors;

			
			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {	
				//printf("i = %d\n", i);
				size_t start_counter = num_rho_neighbors * i;
            	Omega_h::LO radial_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};
				get_radial_bin_neighbors(num_rho_neighbors, mrho_max, test_rho_coordinates[i], radial_bin_neighbors);

				for (int j = 0; j < num_rho_neighbors; ++j){
					neighbors_ids[start_counter + j] = radial_bin_neighbors[j];		
				}
			});


}

void  test_poloidal_bin_neighbors(const size_t num_theta_neighbors, 
								const Omega_h::LOs& radial_bins, 
								const Omega_h::Reals& test_theta_coordinates, 
								Omega_h::Write<Omega_h::LO>& neighbors_ids) {


			
			
			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {	
				//printf("i = %d\n", i);
				size_t start_counter = num_theta_neighbors * i;
            	Omega_h::LO poloidal_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

				get_poloidal_bin_neighbors(num_theta_neighbors, radial_bins[i],
					test_theta_coordinates[i], poloidal_bin_neighbors); 
			

				for (int j = 0; j < num_theta_neighbors; ++j){
					neighbors_ids[start_counter + j] = poloidal_bin_neighbors[j];		
				}
			});

			//Omega_h::read(neighbors_ids);

} 

void test_toroidal_bin_neighbors(const size_t num_zeta_neighbors, 
								 const size_t mtor_max,
								 const double zeta_max,
								 const Omega_h::Reals& test_zeta_coordinates, 
								 Omega_h::Write<Omega_h::LO>& neighbors_ids){
	
			
			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {	
				//printf("i = %d\n", i);
				size_t start_counter = num_zeta_neighbors * i;
            	Omega_h::LO toroidal_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

				get_toroidal_bin_neighbors(num_zeta_neighbors, mtor_max, zeta_max,
						test_zeta_coordinates[i], toroidal_bin_neighbors); 
				for (int j = 0; j < num_zeta_neighbors; ++j){
					neighbors_ids[start_counter + j] = toroidal_bin_neighbors[j];		
				}
			});



	
}

void  test_pitch_bin_neighbors(const size_t num_pitch_neighbors, 
								const Omega_h::Reals& test_pitch_coordinates, const size_t mlmd_max, 
								Omega_h::Write<Omega_h::LO>& neighbors_ids) {

			size_t total_neighbors = 4 * num_pitch_neighbors;

			
			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {	
				size_t start_counter = num_pitch_neighbors * i;
            	Omega_h::LO pitch_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

				get_pitch_bin_neighbors(num_pitch_neighbors, mlmd_max,test_pitch_coordinates[i], pitch_bin_neighbors); 
			

				for (int j = 0; j < num_pitch_neighbors; ++j){
					neighbors_ids[start_counter + j] = pitch_bin_neighbors[j];		
				}
			});
			
	
}

void  test_energy_bin_neighbors(const size_t num_energy_neighbors, 
								const Omega_h::Reals& test_energy_coordinates, const size_t meng_max, 
								double emax, Omega_h::Write<Omega_h::LO>& neighbors_ids) {


			
			
			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {	
				//printf("i = %d\n", i);
				size_t start_counter = num_energy_neighbors * i;
            	Omega_h::LO energy_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

				get_energy_bin_neighbors(num_energy_neighbors, meng_max, emax, test_energy_coordinates[i], energy_bin_neighbors); 
			

				for (int j = 0; j < num_energy_neighbors; ++j){
					neighbors_ids[start_counter + j] = energy_bin_neighbors[j];		
				}
			});	


}

} // namespace

TEST_CASE("gnet bin id evaluation", "[gnet]")
{


		int mrho_max = 14;
		int mtor_max = 90;
		double zeta_max = 1.257;
		int mlmd_max = 33;
		int meng_max = 50;
		double emax = 200000;
		int dim = 5;
		int degree = 2;
		pcms::IntHostMatView slice_lengths_array("array", degree, dim);
		Kokkos::deep_copy(slice_lengths_array, 0);
		pcms::detail::calculate_basis_slice_lengths(slice_lengths_array);
	   auto min_num_supports =  pcms::detail::calculate_basis_vector_size(slice_lengths_array);

		BinNeighbors binNeighbors{};
		binNeighbors.eval(min_num_supports);

		auto num_rho_neighbor_bins = binNeighbors.num_rho_neighbor_bins;
		auto num_theta_neighbor_bins = binNeighbors.num_theta_neighbor_bins;
		auto num_zeta_neighbor_bins = binNeighbors.num_zeta_neighbor_bins;
		auto num_pitch_neighbor_bins = binNeighbors.num_pitch_neighbor_bins;
		auto num_energy_neighbor_bins = binNeighbors.num_energy_neighbor_bins;


		Omega_h::Write<Omega_h::Real> test_rho_coordinates({0.56, 0.8234, 0.92, 0.221}, "TEST RHO COORDINATES");

		Omega_h::Write<Omega_h::Real> test_theta_coordinates({1.4313, 3.0762, 5.224, 0.221}, "TEST THETA COORDINATES");
		
		Omega_h::Write<Omega_h::Real> test_zeta_coordinates({0.0456, 0.123, 0.563, 1.134}, "TEST ZETA COORDINATES");

		Omega_h::Write<Omega_h::Real> test_pitch_coordinates({-0.845, 0.123, -0.38, 0.45}, "TEST PITCH COORDINATES");
		
		Omega_h::Write<Omega_h::Real> test_energy_coordinates({2983.4, 105000.24, 33428, 62346}, "TEST ENERGY COORDINATES");

		
		Omega_h::LOs radial_bins({7, 11, 12, 2}, "TEST RADIAL BINS");
		Omega_h::LOs poloidal_bins({6, 23, 42, 0}, "TEST POLOIDAL BINS");
		Omega_h::LOs toroidal_bins({2, 8, 39, 80}, "TEST TOROIDAL BINS");
		Omega_h::LOs pitch_bins({2, 18, 9, 23}, "TEST PITCH BINS"); 
		Omega_h::LOs energy_bins({0, 25, 7, 15}, "TEST ENERGY BINS"); 
		
		Omega_h::HostRead<Omega_h::LO> host_radial_bins(radial_bins);
		Omega_h::HostRead<Omega_h::LO> host_poloidal_bins(poloidal_bins);
		Omega_h::HostRead<Omega_h::LO> host_toroidal_bins(toroidal_bins);
		Omega_h::HostRead<Omega_h::LO> host_pitch_bins(pitch_bins);
		Omega_h::HostRead<Omega_h::LO> host_energy_bins(energy_bins);
		
		auto pol_max_scan_array = get_pol_max_scan_array(mrho_max);	
		
		Omega_h::HostRead<Omega_h::LO> host_pol_max_scan_array(pol_max_scan_array);
	
		auto mpol_max_array = get_mpol_max_array(mrho_max);

		Omega_h::HostRead<Omega_h::LO> host_mpol_max_array(mpol_max_array);
		


		SECTION("check radial bin id evaluation")
		{
		
			Omega_h::read(test_rho_coordinates);		
	    	Omega_h::Write<Omega_h::LO> d_result_ids(4, 1, "Result Ids");	
			Omega_h::parallel_for(test_rho_coordinates.size(), OMEGA_H_LAMBDA(int i){
				size_t id = search_radial_bin(mrho_max, test_rho_coordinates[i]);
				d_result_ids[i] = id;
		});
	
			Omega_h::HostRead<Omega_h::LO> host_result_ids(d_result_ids);

			auto d_grid_ids = get_grid_ids(radial_bins, poloidal_bins, toroidal_bins, mrho_max, mtor_max);
		
			for (int i = 0; i < host_result_ids.size(); ++i){
				REQUIRE(host_result_ids[i] == host_radial_bins[i]);
			}
  		}

		SECTION("check poloidal bin id evaluation")
		{
			
//			Omega_h::read(test_theta_coordinates);		
			

			Omega_h::Write<Omega_h::LO> d_result_ids(4, 1, "Result Ids");	
			Omega_h::parallel_for(test_theta_coordinates.size(), OMEGA_H_LAMBDA(int i){
					size_t id = search_poloidal_bin(radial_bins[i], test_theta_coordinates[i]);
					d_result_ids[i] = id;
			});
			
			Omega_h::HostRead<Omega_h::LO> host_result_ids(d_result_ids);

				
			for (int i = 0; i < host_result_ids.size(); ++i){
				REQUIRE(host_result_ids[i] == host_poloidal_bins[i]);
			}
		}


	
		SECTION("check toroidal bin id evaluation")
		{	
			Omega_h::read(test_zeta_coordinates);		
			Omega_h::Write<Omega_h::LO> d_result_ids(4, 1, "Result Ids");	
			Omega_h::parallel_for(test_zeta_coordinates.size(), OMEGA_H_LAMBDA(int i){
					size_t id = search_toroidal_bin(mtor_max, zeta_max, test_zeta_coordinates[i]);
					d_result_ids[i] = id;
			});
		
			Omega_h::HostRead<Omega_h::LO> host_result_ids(d_result_ids);

			
			for (int i = 0; i < host_result_ids.size(); ++i){
				REQUIRE(host_result_ids[i] == host_toroidal_bins[i]);
			}
		}
		
		
		SECTION("check pitch bin id evaluation")
		{	
			Omega_h::read(test_pitch_coordinates);		
			Omega_h::Write<Omega_h::LO> d_result_ids(4, 1, "Result Ids");	
			Omega_h::parallel_for(test_pitch_coordinates.size(), OMEGA_H_LAMBDA(int i){
					size_t id = search_pitch_bin(mlmd_max, test_pitch_coordinates[i]);
					d_result_ids[i] = id;
			});
		
			Omega_h::HostRead<Omega_h::LO> host_result_ids(d_result_ids);
			
			for (int i = 0; i < host_result_ids.size(); ++i){
				REQUIRE(host_result_ids[i] == host_pitch_bins[i]);
			}
		}



		SECTION("check energy bin id evaluation")
		{	
			Omega_h::read(test_energy_coordinates);		
			Omega_h::Write<Omega_h::LO> d_result_ids(4, 1, "Result Ids");	
			Omega_h::parallel_for(test_energy_coordinates.size(), OMEGA_H_LAMBDA(int i){
					size_t id = search_energy_bin(meng_max, emax, test_energy_coordinates[i]);
					d_result_ids[i] = id;
			});
		
			Omega_h::HostRead<Omega_h::LO> host_result_ids(d_result_ids);

			for (int i = 0; i < host_result_ids.size(); ++i){
				REQUIRE(host_result_ids[i] == host_energy_bins[i]);
			}
		}

		SECTION(" check grid id evaluation")
		{
			
			int expected_ids[4] = {0, 0, 0, 0};

			auto d_grid_ids = get_grid_ids(radial_bins, poloidal_bins, toroidal_bins, mrho_max, mtor_max);
	
			Omega_h::HostRead<Omega_h::LO> host_grid_ids(d_grid_ids);

			for (int i = 0; i < 4; ++i){
				int rho_bin_id = host_radial_bins[i];
				auto mpol_sum_at_radial_bin = host_pol_max_scan_array[rho_bin_id];
				auto theta_bin_id = host_poloidal_bins[i];
				auto zeta_bin_id = host_toroidal_bins[i];

				expected_ids[i] = (mpol_sum_at_radial_bin + theta_bin_id) * mtor_max + zeta_bin_id;
				REQUIRE(host_grid_ids[i] == expected_ids[i]); 
			}
			
			
		}

		SECTION(" check marker id evaluation")
		{
			
			
			auto d_grid_ids = get_grid_ids(radial_bins, poloidal_bins, toroidal_bins, mrho_max, mtor_max);
			
			Omega_h::Write<Omega_h::LO> d_results_ids(4, 1, "Result Ids");
			
			int expected_ids[4] = {0, 0, 0, 0};

	
			Omega_h::HostRead<Omega_h::LO> host_grid_ids(d_grid_ids);

			Omega_h::parallel_for(4, OMEGA_H_LAMBDA(int i) {
				int mpol_sum_at_current_radial = pol_max_scan_array[radial_bins[i]];
				auto marker_id = eval_marker_id(mlmd_max, meng_max, d_grid_ids[i], pitch_bins[i], energy_bins[i] ); 	
				d_results_ids[i] = marker_id;	

			  });
			
			Omega_h::HostRead<Omega_h::LO> host_results_ids(d_results_ids);

			for (int i = 0; i < 4; ++i){
				auto grid_id = host_grid_ids[i];
				auto pitch_bin_id = host_pitch_bins[i];
				auto energy_bin_id = host_energy_bins[i];
				
				expected_ids[i] = grid_id * mlmd_max * meng_max +
								energy_bin_id * mlmd_max + pitch_bin_id;
									
				REQUIRE(host_results_ids[i] == expected_ids[i]); 
			}
				
		}




		SECTION("check radial bin neighbors")
		{	
			Omega_h::Write<Omega_h::LO> neighbors_ids(32, 0, "poloidal bin neighbors");
			test_radial_bin_neighbors(8, Omega_h::read(test_rho_coordinates), mrho_max, neighbors_ids);

			int expected_ids[32] = {7, 8, 6, 9, 5, 10, 4, 11,
									11, 12, 10, 13, 9, 8, 7, 6,
									12, 13, 11, 10, 9, 8, 7, 6,
									2, 3, 1, 4, 0, 5, 6, 7};
			
			Omega_h::HostRead<Omega_h::LO> host_neighbors_ids(neighbors_ids);
			int total_neighbors = 32;
			for (int i = 0; i < total_neighbors; ++i){
				//printf("i = %d, neighbors_ids[%d] = %d\n", i, i, host_neighbors_ids[i]); 
				CAPTURE(i,host_neighbors_ids[i], expected_ids[i]);
				REQUIRE(host_neighbors_ids[i] == expected_ids[i]);
			}

		}


		SECTION("check poloidal bin neighbors")
		{
			Omega_h::Write<Omega_h::LO> neighbors_ids(32, 0, "poloidal bin neighbors");
			test_poloidal_bin_neighbors(8, radial_bins, Omega_h::read(test_theta_coordinates), neighbors_ids);
			
			int expected_ids[32] = {6, 7, 5, 8, 4, 9, 3, 10,
									23, 24, 22, 25, 21, 26, 20, 27,
									42, 43, 41, 44, 40, 45, 39, 46,
									0, 1, 11, 2, 10, 3, 9, 4};
			
			Omega_h::HostRead<Omega_h::LO> host_neighbors_ids(neighbors_ids);
			//std::cout<<"HOST NEIGHBORS\n";
		//	for(int i=0; i<host_neighbors_ids.size(); ++i) {
		//		std::cout<<host_neighbors_ids[i]<<" ";
		//	}
		//	std::cout<<std::endl;
			int total_neighbors = 32;
			for (int i = 0; i < total_neighbors; ++i){
				//printf("i = %d, neighbors_ids[%d] = %d\n", i, i, host_neighbors_ids[i]); 
				CAPTURE(i,host_neighbors_ids[i], expected_ids[i]);
				REQUIRE(host_neighbors_ids[i] == expected_ids[i]);
			}

		}
		
		
	SECTION("check toroidal bin neighbors")
		{
			Omega_h::Write<Omega_h::LO> neighbors_ids(32, 0, "poloidal bin neighbors");
			
			test_toroidal_bin_neighbors(8, mtor_max, zeta_max, Omega_h::read(test_zeta_coordinates), neighbors_ids);
			
			int expected_ids[32] = {2, 3, 1, 4, 0, 5, 89, 6,
									8, 9, 7, 10, 6, 11, 5, 12,
									39, 40, 38, 41, 37, 42, 36, 43, 
                                    80, 81, 79, 82, 78, 83, 77, 84};
			
			Omega_h::HostRead<Omega_h::LO> host_neighbors_ids(neighbors_ids);

			size_t total_neighbors =  32;
			for (int i = 0; i < total_neighbors; ++i){
				//printf("i = %d, neighbors_ids[%d] = %d\n", i, i, host_neighbors_ids[i]); 
				CAPTURE(i,host_neighbors_ids[i], expected_ids[i]);
				REQUIRE(host_neighbors_ids[i] == expected_ids[i]);
			}
	
		}

		SECTION("check pitch bin neighbors")
		{
			Omega_h::Write<Omega_h::LO> neighbors_ids(32, 0, "poloidal bin neighbors");
			test_pitch_bin_neighbors(8, Omega_h::read(test_pitch_coordinates), mlmd_max, neighbors_ids);
	
			int expected_ids[32] = {2, 3, 1, 4, 0, 5, 6, 7, 
									18, 19, 17, 20, 16, 21, 15, 22,
									9, 10, 8, 11, 7, 12, 6, 13, 
									23, 24, 22, 25, 21, 26, 20, 27};
			
			Omega_h::HostRead<Omega_h::LO> host_neighbors_ids(neighbors_ids);

			size_t total_neighbors = 32;
			for (int i = 0; i < total_neighbors; ++i){
				//printf("i = %d, neighbors_ids[%d] = %d\n", i, i, host_neighbors_ids[i]); 
				CAPTURE(i,host_neighbors_ids[i], expected_ids[i]);
				REQUIRE(host_neighbors_ids[i] == expected_ids[i]);
			}
			
		}


		SECTION("check energy bin neighbors")
		{
			
			Omega_h::Write<Omega_h::LO> neighbors_ids(32, 0, "poloidal bin neighbors");
			 test_energy_bin_neighbors(8, Omega_h::read(test_energy_coordinates), meng_max, 
								emax,  neighbors_ids);
			
			 int expected_ids[32] = {0, 1, 2, 3, 4, 5, 6, 7, 
									25, 26, 24, 27, 23, 28, 22, 29, 
									7, 8, 6, 9, 5, 10, 4, 11, 
									15, 16, 14, 17, 13, 18, 12, 19};
			
			 Omega_h::HostRead<Omega_h::LO> host_neighbors_ids(neighbors_ids);
			size_t total_neighbors = 32;

			for (int i = 0; i < total_neighbors; ++i){
				//printf("i = %d, neighbors_ids[%d] = %d\n", i, i, host_neighbors_ids[i]); 
				CAPTURE(i,host_neighbors_ids[i], expected_ids[i]);
				REQUIRE(host_neighbors_ids[i] == expected_ids[i]);
			}
	
			
		}
}

// Reads the GNET file named by the GNET_TEST_FILE environment variable;
// skipped when it is not set.
TEST_CASE("gnet HDF5 reader", "[gnet][hdf5]")
{
  const char* test_file = std::getenv("GNET_TEST_FILE");
  if (test_file == nullptr) {
    SKIP("GNET_TEST_FILE not set");
  }

  GnetHDF5Reader reader(test_file);
  const size_t maxEnergyBins = reader.getMaxEnergyBins();
  const size_t maxRhoBins = reader.getMaxRhoBins();
  const size_t maxToroidalBins = reader.getMaxToroidalBins();
  const size_t maxPitchBins = reader.getMaxPitchBins();

  REQUIRE(maxEnergyBins > 0);
  REQUIRE(maxRhoBins > 0);
  REQUIRE(maxToroidalBins > 0);
  REQUIRE(maxPitchBins > 0);
  REQUIRE(reader.getMaxEnergyEV() > 0.0);
  REQUIRE(reader.getMaxZetaRad() > 0.0);

  // f5d is (grid, energy, pitch) with grid = sum_rho(4 (rho + 1)) * zeta
  const size_t grid_total = 2 * maxRhoBins * (maxRhoBins + 1) * maxToroidalBins;
  auto f5d = reader.getF5D();
  REQUIRE(f5d.extent(0) == grid_total);
  REQUIRE(f5d.extent(1) == maxEnergyBins);
  REQUIRE(f5d.extent(2) == maxPitchBins);
  REQUIRE(reader.getVolume3D().extent(0) == grid_total);
}
