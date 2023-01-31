#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include <wdmcpl/omega_h_field.h>
#include <wdmcpl/xgc_field_adapter.h>

using wdmcpl::ConstructRCFromOmegaHMesh;
using wdmcpl::Copy;
using wdmcpl::CouplerClient;
using wdmcpl::CouplerServer;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;
using wdmcpl::GO;
using wdmcpl::Lagrange;
using wdmcpl::make_array_view;
using wdmcpl::OmegaHField;
using wdmcpl::OmegaHFieldAdapter;
using wdmcpl::ReadReverseClassificationVertex;
using wdmcpl::ReverseClassificationVertex;

namespace ts = test_support;

// TODO: we should communicate the geometric ids in the overlap regions.
// is there a way to use the isOverlap functor to do this. This allows for
// maximum flexibility moving forward

void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file, int nphi)
{
  wdmcpl::CouplerServer cpl("xgc_n0_coupling", comm,
                            ts::setupServerPartition(mesh, cpn_file), mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  ReverseClassificationVertex rc;
  std::string numbering = "simNumbering";
  WDMCPL_ALWAYS_ASSERT(mesh.has_tag(0, numbering));

  auto is_overlap = ts::markServerOverlapRegion(
    mesh, partition, [](const int dim, const int id) {
      if (id >= 1 && id <= 2) {
        return 1;
      }
      return 0;
    });
  std::vector<wdmcpl::ConvertibleCoupledField*> potential_fields(nphi);
  std::cerr<<"ADDING FIELDS\n";
  for (int i = 0; i < nphi; ++i) {
    std::stringstream field1_name;
    field1_name << "dpot_plane_" << i;
    std::cerr<<field1_name.str()<<"\n";
    potential_fields[i] = cpl.AddField(field1_name.str(),
                 wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>(
                   field1_name.str(), mesh, is_overlap, numbering),
                 FieldTransferMethod::Copy, // to Omega_h
                 FieldEvaluationMethod::None,
                 FieldTransferMethod::Copy, // from Omega_h
                 FieldEvaluationMethod::None, is_overlap,"deltaf1/");
    std::cerr<<field1_name.str()<<" DONE\n";
  }
  std::cerr<<"SEND/Recv loop started\n";
  while (true) {
    // TODO: do this w/o blocking
    for (auto* field : potential_fields) {
      WDMCPL_ALWAYS_ASSERT(field != nullptr);
      field->Receive();
    }
    for (auto* field : potential_fields) {
      field->Send();
    }
  }
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  int size = world->size();
  if (argc != 4) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << "</path/to/omega_h/mesh> "
                   "</path/to/partitionFile.cpn> "
                   "sml_nphi_total";
    }
    exit(EXIT_FAILURE);
  }

  const auto meshFile = argv[1];
  const auto classPartitionFile = argv[2];
  const int sml_nphi_total = std::atoi(argv[3]);

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  omegah_coupler(mpi_comm, mesh, classPartitionFile, sml_nphi_total);
  return 0;
}
