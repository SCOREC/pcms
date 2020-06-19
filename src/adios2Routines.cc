#include "adios2Routines.h"

namespace coupler {


  Array2d<CV>* receive_density(const std::string cce_folder,
      const adios2_handler &handler,GO my_start[2],GO my_count[2],MPI_Comm comm) { 
      adios2::IO io = handler.IO; 
      adios2::Engine engine = handler.eng;
    const std::string name = "gene_density";
    return receive2d_from_ftn<CV>(cce_folder,name, io, engine,my_start,my_count,comm);
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
      const adios2_handler &handler,GO my_start[2], GO my_count[2], MPI_Comm comm) { 
      adios2::IO io = handler.IO; 
      adios2::Engine engine = handler.eng;
    const std::string name = "xgc_field";
    return receive2d_from_ftn<double>(cce_folder,name, io, engine,my_start,my_count, comm);
  }

  void send_field(const std::string cce_folder, const Array2d<CV>* field,
      const adios2_handler &handler, adios2::Variable<CV> &send_id) {
    adios2::IO io = handler.IO; 
    adios2::Engine engine = handler.eng;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string fld_name = "cpl_field";
    send2d_from_C(field, cce_folder, fld_name, io, engine, send_id);
    std::cerr << rank <<  ": send " << fld_name <<" done \n";
  }
 
 void AdiosProTransFortranCpp2D(LO rankout,LO mypxout,LO mypyout, const LO mypex,const LO mypey,const LO npx,const LO npy)
{
  mypxout = mypey;
  mypyout = mypex; 
  LO npxbar;
  npxbar=npy; 
  rankout = mypyout*npxbar+mypxout; 
}

}
