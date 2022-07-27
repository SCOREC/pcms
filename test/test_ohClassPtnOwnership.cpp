#include <Omega_h_mesh.hpp>
#include <iostream>
#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName)
{
  namespace ts = test_support;
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank()
                               ? ts::readClassPartitionFile(cpnFileName)
                               : ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  MPI_Barrier(MPI_COMM_WORLD);
  REDEV_ALWAYS_ASSERT(mesh.nelems()); //all ranks should have elements
  Omega_h::vtk::write_parallel("rdv.vtk", &mesh, 2);
  MPI_Barrier(MPI_COMM_WORLD);
  if(!ohComm->rank()) std::cerr << " before CreateClassificationPartition\n";
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

void client0(MPI_Comm comm, Omega_h::Mesh& mesh)
{

  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Client, comm,
                      redev::ClassPtn{});
}
void client1(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Client, comm,
                      redev::ClassPtn{});
}
void coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{

  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Server, comm,
                      setupServerPartition(mesh, cpn_file));
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 4) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                   "/path/to/partitionFile.cpn\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  const std::string name = "meshVtxIds";
  switch (clientId) {
    case -1: coupler(mpi_comm, mesh, classPartitionFile); break;
    case 0: client0(mpi_comm, mesh); break;
    case 1: client1(mpi_comm, mesh); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  return 0;
}
