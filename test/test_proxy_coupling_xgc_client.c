#include <wdmcpl/capi/client.h>
#include <wdmcpl/capi/kokkos.h>
#include <mpi.h>
#include <stdlib.h>
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
  wdmcpl_kokkos_initialize_without_args();
  int world_rank = -1, world_size = -1, plane_rank = -1, plane_size = -1,
      client_rank = -1, client_size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const int nplanes = 2;
  if (world_size % nplanes != 0) {
    printf("Number of ranks must be divisible by number of planes\n");
    abort();
  }
  int plane = world_rank % nplanes;
  printf("PLANE %d\n", plane);
  // client comm is NULL on any ranks that shouldn't participate in
  // CouplingCommunication
  MPI_Comm client_comm, plane_comm;
  MPI_Comm_split(MPI_COMM_WORLD, plane, world_rank, &plane_comm);
  MPI_Comm_rank(plane_comm, &plane_rank);
  MPI_Comm_size(plane_comm, &plane_size);
  MPI_Comm_split(MPI_COMM_WORLD, (plane_rank == 0) ? 0 : MPI_UNDEFINED,
                 world_rank, &client_comm);
  if (client_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);
  }
  printf("world: %d %d; plane: %d %d ; client: %d %d\n", world_rank, world_size,
         plane_rank, plane_size, client_rank, client_size);

  WdmCplClientHandle* client =
    wdmcpl_create_client("proxy_couple", client_comm);
  const char* rc_file = argv[1];
  WdmCplReverseClassificationHandle* rc =
    wdmcpl_load_reverse_classification(rc_file, MPI_COMM_WORLD);
  const int nverts = wdmcpl_reverse_classification_count_verts(rc);
  // long int* data[nplanes];
  long int* data = calloc(nverts, sizeof(long int));
  WdmCplFieldHandle* field[nplanes];
  WdmCplFieldAdapterHandle* field_adapters[nplanes];
  for (int i = 0; i < nplanes; ++i) {
    char field_name[100];
    sprintf(field_name, "xgc_gids_plane_%d", i);
    printf("%s\n", field_name);
    int communicating_rank = (i == plane) && (plane_rank == 0);
    if (plane == i) {
      field_adapters[i] = wdmcpl_create_xgc_field_adapter(
        "adapter1", plane_comm, data, nverts, WDMCPL_LONG_INT, rc, &in_overlap);
    } else {
      field_adapters[i] = wdmcpl_create_dummy_field_adapter();
    }
    field[i] = wdmcpl_add_field(client, field_name, field_adapters[i],
                                communicating_rank);
  }
  // only set the data on plane_rank 0 so that we can verify that XGC Adapter
  // properly broadcasting to the other ranks.
  if (plane_rank == 0) {
    for (int i = 0; i < nverts; ++i) {
      data[i] = i;
    }
  }
  wdmcpl_begin_send_phase(client);
  wdmcpl_send_field(field[plane]);
  wdmcpl_end_send_phase(client);
  wdmcpl_begin_receive_phase(client);
  wdmcpl_receive_field(field[plane]);
  wdmcpl_end_receive_phase(client);
  // check data on all ranks. This should be set either from RDV communication
  // or the broadcast in the XGC Field Adapter
  for (int i = 0; i < nverts; ++i) {
    if (data[i] != i) {
      printf("ERROR: data[%d] = %ld, should be %d", i, data[i], i);
      abort();
    }
  }
  // only set the data on plane_rank 0 so that we can verify that XGC Adapter
  // properly broadcasting to the other ranks.
  if (plane_rank == 0) {
    for (int i = 0; i < nverts; ++i) {
      data[i] *= 2;
    }
  }
  wdmcpl_begin_send_phase(client);
  wdmcpl_send_field(field[plane]);
  wdmcpl_end_send_phase(client);
  wdmcpl_begin_receive_phase(client);
  wdmcpl_receive_field(field[plane]);
  wdmcpl_end_receive_phase(client);
  // check data on all ranks. This should be set either from RDV communication
  // or the broadcast in the XGC Field Adapter
  for (int i = 0; i < nverts; ++i) {
    if (data[i] != 2 * i) {
      printf("ERROR: data[%d] = %ld, should be %d", i, data[i], 2 * i);
      abort();
    }
  }
  for (int i = 0; i < nplanes; ++i) {
    wdmcpl_destroy_field_adapter(field_adapters[i]);
  }
  free(data);
  wdmcpl_destroy_reverse_classification(rc);
  wdmcpl_destroy_client(client);
  if(client_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&client_comm);
  }
  MPI_Comm_free(&plane_comm);
  wdmcpl_kokkos_finalize();
  MPI_Finalize();
  return 0;
}
