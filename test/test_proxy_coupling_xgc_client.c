#include <wdmcpl/capi/client.h>
#include <mpi.h>

int8_t in_overlap(int dimension, int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if (dimension == 2 && (id >= 22 && id <= 34)) {
    return 1;
  } else if (dimension == 1 && (id >= 21 && id <= 34)) {
    return 1;
  } else if (dimension == 0 && (id >= 21 && id <= 34)) {
    return 1;
  }
  return 0;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  WdmCplClient* client = wdmcpl_create_client("proxy_couple", MPI_COMM_WORLD);
  const int data_size = 1000;
  double data[data_size];
  WdmCplReverseClassificationHandle* rc =
    wdmcpl_load_reverse_classification("/path/to/rc_file.txt", MPI_COMM_WORLD);
  WdmCplFieldAdapter* xgc_adapter = wdmcpl_create_xgc_field_adapter(
    "adapter1", &data, data_size, WDMCPL_DOUBLE, rc, &in_overlap);
  wdmcpl_add_field(client, "xgc_field", xgc_adapter);
  wdmcpl_send_field(client, "xgc_field");
  wdmcpl_receive_field(client, "xgc_field");

  wdmcpl_destroy_field_adapter(xgc_adapter);
  wdmcpl_destroy_reverse_classification(rc);
  wdmcpl_destroy_client(client);

  MPI_Finalize();
  return 0;
}