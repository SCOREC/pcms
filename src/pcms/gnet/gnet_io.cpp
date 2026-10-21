#include <pcms/gnet/gnet_io.hpp>
#include <Omega_h_for.hpp>
#include <hdf5.h>
#include <stdexcept>
#include <type_traits>
#include <vector>

static_assert(std::is_same_v<hid_t, std::int64_t>,
              "gnet_io stores HDF5 ids as std::int64_t");

namespace pcms::gnet
{

namespace
{

std::vector<hsize_t> getDatasetDims(hid_t file_id, const std::string& path)
{
  hid_t dataset_id = H5Dopen(file_id, path.c_str(), H5P_DEFAULT);
  if (dataset_id < 0) {
    throw std::runtime_error("Error opening dataset: " + path);
  }

  hid_t dataspace_id = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(dataspace_id);
  std::vector<hsize_t> dims(rank);
  H5Sget_simple_extent_dims(dataspace_id, dims.data(), NULL);

  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  return dims;
}

void readDataset(hid_t file_id, const std::string& path, hid_t mem_type,
                 void* buffer)
{
  hid_t dataset_id = H5Dopen(file_id, path.c_str(), H5P_DEFAULT);
  if (dataset_id < 0) {
    throw std::runtime_error("Error opening dataset: " + path);
  }
  H5Dread(dataset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  H5Dclose(dataset_id);
}

} // namespace

GnetHDF5Reader::GnetHDF5Reader(const std::string& filename)
{
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    throw std::runtime_error("Error opening HDF5 file: " + filename);
  }

  // Read scalar values
  readDataset(file_id, "/distribution5D/dimensions/mxenergy", H5T_NATIVE_INT,
              &meng_max);
  readDataset(file_id, "/distribution5D/dimensions/mxrho", H5T_NATIVE_INT,
              &mrho_max);
  readDataset(file_id, "/distribution5D/dimensions/mxzeta", H5T_NATIVE_INT,
              &mtor_max);
  readDataset(file_id, "/distribution5D/dimensions/mxpitch", H5T_NATIVE_INT,
              &mlmd_max);
  readDataset(file_id, "/distribution5D/range/engmax_eV", H5T_NATIVE_DOUBLE,
              &emax);
  readDataset(file_id, "/distribution5D/range/zetamax_rad", H5T_NATIVE_DOUBLE,
              &zeta_max);

  // Read 5D Distribution Data
  auto f5d_dims = getDatasetDims(file_id, "/distribution5D/data/f5d");
  if (f5d_dims.size() != 3) {
    throw std::runtime_error("Expected 3D dataset /distribution5D/data/f5d in " +
                             filename);
  }
  host_f5d = Kokkos::View<double***, Kokkos::HostSpace>(
    "f5d data", f5d_dims[0], f5d_dims[1], f5d_dims[2]);
  readDataset(file_id, "/distribution5D/data/f5d", H5T_NATIVE_DOUBLE,
              host_f5d.data());

  // Read Volume Data
  auto volume_dims = getDatasetDims(file_id, "/distribution5D/data/volume3d");
  if (volume_dims.size() != 1) {
    throw std::runtime_error(
      "Expected 1D dataset /distribution5D/data/volume3d in " + filename);
  }
  host_volume =
    Kokkos::View<double*, Kokkos::HostSpace>("volume data", volume_dims[0]);
  readDataset(file_id, "/distribution5D/data/volume3d", H5T_NATIVE_DOUBLE,
              host_volume.data());
}

GnetHDF5Reader::~GnetHDF5Reader()
{
  H5Fclose(file_id);
}

HDF5Writer::HDF5Writer(const std::string& filename)
{
  file_id =
    H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    throw std::runtime_error("Failed to create HDF5 file: " + filename);
  }
}

HDF5Writer::~HDF5Writer()
{
  if (file_id >= 0) {
    H5Fclose(file_id);
  }
}

void HDF5Writer::writef2D(const std::string& dataset_name,
                          const Kokkos::View<double**, Kokkos::HostSpace>& data)
{
  hsize_t dims[2] = {static_cast<hsize_t>(data.extent(0)),
                     static_cast<hsize_t>(data.extent(1))};

  hid_t dataspace = H5Screate_simple(2, dims, nullptr);
  hid_t dataset = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           data.data());

  H5Sclose(dataspace);
  H5Dclose(dataset);
}

void write_gtc_output(const std::string& path, const HostMatrix& markers)
{
  HDF5Writer writer(path);
  writer.writef2D("/gtc_n6d", markers);
}

Omega_h::Reals calculate_number_density(
  const size_t mrho_max, const Omega_h::LOs& mpol_max, const size_t mtor_max,
  const Omega_h::LOs& pol_max_scan_array, const size_t meng_max,
  const size_t mlmd_max, const Kokkos::View<double***, Kokkos::HostSpace>& f5d,
  const Kokkos::View<double*, Kokkos::HostSpace>& volume3d,
  const int total_num_particles)
{

  Omega_h::HostRead<Omega_h::LO> host_mpol_max(mpol_max);
  Omega_h::HostRead<Omega_h::LO> host_pol_max_scan_array(pol_max_scan_array);
  Kokkos::View<double*, Kokkos::HostSpace> host_number_density(
    "stores number density in host view", total_num_particles);
  int count = -1;
  for (size_t mrho = 0; mrho < mrho_max; ++mrho) {
    for (size_t mtheta = 0; mtheta < host_mpol_max[mrho]; ++mtheta) {
      for (size_t mzeta = 0; mzeta < mtor_max; ++mzeta) {
        size_t m3d = (host_pol_max_scan_array[mrho] + mtheta) * mtor_max + mzeta;
        for (size_t meng = 0; meng < meng_max; ++meng) {
          for (size_t mlmd = 0; mlmd < mlmd_max; ++mlmd) {
            count += 1;
            host_number_density(count) = f5d(m3d, meng, mlmd) / volume3d(m3d);
          }
        }
      }
    }
  }

  using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;
  Kokkos::View<double*, DeviceSpace> device_number_density(
    "particle density", total_num_particles);

  // copy from host to actual device view
  Kokkos::deep_copy(device_number_density, host_number_density);

  Omega_h::Write<Omega_h::Real> number_density(total_num_particles, 0,
                                               "stores number density");

  Omega_h::parallel_for(
    total_num_particles,
    OMEGA_H_LAMBDA(int i) { number_density[i] = device_number_density(i); });
  return Omega_h::read(number_density);
}

} // namespace pcms::gnet
