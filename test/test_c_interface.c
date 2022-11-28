#include <wdmcpl/capi/client.h>
#include <mpi.h>
#include <stddef.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  WdmCplClient* client = wdmcpl_create_client("proxy_couple", MPI_COMM_WORLD);

  WdmCplFieldAdapter* xgc_adapter = NULL;//wdmcpl_create_xgc_field_adapter(
    //"adapter1",
    //);
  wdmcpl_add_field(client, "xgc_field", xgc_adapter);
  wdmcpl_send_field(client, "xgc_field");
  wdmcpl_receive_field(client, "xgc_field");

  wdmcpl_destroy_client(client);

  MPI_Finalize();
  return 0;
}