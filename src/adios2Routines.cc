#include "adios2Routines.h"

namespace coupler {


  Array2d<double>* receive_density(const std::string cce_folder,
      const adios2_handler &handler, MPI_Comm comm) { 
      adios2::IO io = handler.IO; 
      adios2::Engine engine = handler.eng;
    const std::string name = "gene_density";
    return receive2d_from_ftn<double>(cce_folder,name, io, engine, comm);
  }

  void send_density(const std::string cce_folder, const Array2d<double>* density,
      const adios2_handler &handler, adios2::Variable<double> &send_id) {
    adios2::IO io = handler.IO; 
    adios2::Engine engine = handler.eng;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string fld_name = "cpl_density";
    send2d_from_C(density, cce_folder, fld_name, io, engine, send_id);
    std::cerr << rank <<  ": send " << fld_name <<" done \n";
  }

  Array2d<double>* receive_field(const std::string cce_folder,
      const adios2_handler &handler, MPI_Comm comm) { 
      adios2::IO io = handler.IO; 
      adios2::Engine engine = handler.eng;
    const std::string name = "xgc_field";
    return receive2d_from_ftn<double>(cce_folder,name, io, engine, comm);
  }

  void send_field(const std::string cce_folder, const Array2d<double>* field,
      const adios2_handler &handler, adios2::Variable<double> &send_id) {
    adios2::IO io = handler.IO; 
    adios2::Engine engine = handler.eng;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string fld_name = "cpl_field";
    send2d_from_C(field, cce_folder, fld_name, io, engine, send_id);
    std::cerr << rank <<  ": send " << fld_name <<" done \n";
  }

}
