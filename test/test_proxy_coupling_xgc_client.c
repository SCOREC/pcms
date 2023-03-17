#include <wdmcpl/capi/client.h>
#include <mpi.h>
#include <stdlib.h>
#include <flcl-util-cxx.h>
#include <printf.h>

int8_t in_overlap(int dimension, int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if ((id >= 22 && id <= 34)) {
    if (dimension == 2) {
      return 1;
    } else if (dimension == 1) {
      return 1;
    } else if (dimension == 0) {
      return 1;
    }
  }
  return 0;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  c_kokkos_initialize_without_args();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int rank_participates = (rank ==0);
  // client comm is NULL on any ranks that shouldn't participate in CouplingCommunication
  MPI_Comm client_comm;
  MPI_Comm_split(MPI_COMM_WORLD,rank_participates?0:MPI_UNDEFINED,rank,&client_comm);
  WdmCplClientHandle* client =
    wdmcpl_create_client("proxy_couple", client_comm);
  const char* rc_file = argv[1];
  WdmCplReverseClassificationHandle* rc =
    wdmcpl_load_reverse_classification(rc_file, MPI_COMM_WORLD);
  const int nverts = wdmcpl_reverse_classification_count_verts(rc);
  long int* data = calloc(nverts, sizeof(long int));
  WdmCplFieldAdapterHandle* xgc_adapter = wdmcpl_create_xgc_field_adapter(
    "adapter1", MPI_COMM_WORLD, data, nverts, WDMCPL_LONG_INT, rc, &in_overlap);
  WdmCplFieldHandle* field = wdmcpl_add_field(client, "xgc_gids", xgc_adapter);
  // only set the data on rank 0 so that we can verify that XGC Adapter
  // properly broadcasting to the other ranks.
  if (rank == 0) {
    for (int i = 0; i < nverts; ++i) {
      data[i] = i;
    }
  }
  wdmcpl_begin_send_phase(client);
  wdmcpl_send_field(field);
  wdmcpl_end_send_phase(client);
  wdmcpl_begin_receive_phase(client);
  wdmcpl_receive_field(field);
  wdmcpl_end_receive_phase(client);
  for (int i = 0; i < nverts; ++i) {
    if (data[i] != i) {
      printf("ERROR: data[%d] = %ld, should be %d", i, data[i], i);
      abort();
    }
  }
  for (int i = 0; i < nverts; ++i) {
    data[i] *= 2;
  }
  wdmcpl_begin_send_phase(client);
  wdmcpl_send_field_name(client, "xgc_gids");
  wdmcpl_end_send_phase(client);
  wdmcpl_begin_receive_phase(client);
  wdmcpl_receive_field_name(client, "xgc_gids");
  wdmcpl_end_receive_phase(client);
  for (int i = 0; i < nverts; ++i) {
    if (data[i] != 2 * i) {
      printf("ERROR: data[%d] = %ld, should be %d", i, data[i], 2 * i);
      abort();
    }
  }

  wdmcpl_destroy_field_adapter(xgc_adapter);
  free(data);
  wdmcpl_destroy_reverse_classification(rc);
  wdmcpl_destroy_client(client);
  c_kokkos_finalize();
  MPI_Finalize();
  return 0;
}
