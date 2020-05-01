#ifndef COUPLING_H_
#define COUPLING_H_

#include <adios2.h>
#include <iostream>
#include <cassert>
#include <Kokkos_Core.hpp> //not used
#include <typeinfo> //not used
#include <string>
#include <fstream>
#include "couplingConstants.h" //not used
#include "couplingTypes.h"

namespace coupler { 
  /** Storage of double precision 2D array data
   *  and associated meta data
   */
  template<class T> 
  class Array2d {
    public:
      Array2d(GO gH, GO gW, GO lH, GO lW, GO start) :
        globH(gH), globW(gW), locH(lH), locW(lW), locFirstCol(start) {
          vals = new T[locH*locW];
      }   
      ~Array2d() {
        globH = globW = locH = locW = 0;
        delete [] vals;
      }   
      T val(long i) const {
        assert(i<(locH*locW));
        return vals[i];
      }   
      T* data() const { return vals; };
      GO globalH() const { return globH; };
      GO globalW() const { return globW; };
      GO localH() const { return locH; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      T* vals;
      GO globH;
      GO globW;
      GO locH;
      GO locW;
      GO locFirstCol;
  };

  template<class T>
  class Array1d {
    public:
      Array1d(GO gW, GO lW, GO start) :
        globW(gW), locW(lW), locFirstCol(start) {
        vals = new T[locW];
      }
      ~Array1d() {
        globW = locW = 0;
        delete [] vals;
      }
      T val(long i) const {
        assert(i<(locW));
        return vals[i];
      }
      T* data() const { return vals; };
      GO globalW() const { return globW; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      T* vals;
      GO globW;
      GO locW;
      GO locFirstCol;
  };
 
  /** Destroy Array1d and Array2d objects
   */
  template<typename T>
    void destroy(Array1d<T>* a) {
      delete a;
    }

  template<typename T>
    void destroy(Array2d<T>* a) {
      delete a;
    }






  /** Sanity check values
   */
  template<typename T>
  void printSomeDensityVals(const Array2d<T>* density) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //asserting the density values received from GENE
    if(!rank)
    {
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank <<  ": first 10 density at "<< i
          << " is "<< density->val(i) <<"\n";
      }
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank << ": first 10 for rank 1 at: [67236]" << " + "<< i
          << " is " << density->val(67236 + i) << "\n";
      }
    }
  
    if(rank == 1)
    {
      for (int i = 0; i < 10; i++)
      {
        int offset = ((density->localW() - 1) * density->localH()) + 67235 - 9;
        std::cerr << rank << ": last 10 for rank 0 at: [67235 - 9]" << " + "<< i
          << " is " << density->val(offset  + i) << "\n";
      }
      int last_ten = (density->localH() * density->localW()) - 10;
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank <<  ": last 10 density at " << last_ten + i << " is "
          << density->val(last_ten + i) <<"\n";
      }
    }
  }
  
  
  template<typename T> 
  Array1d<T>* receive1d_from_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
    const std::string fname = dir + "/" + name + ".bp";
  
    if(!eng){
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "480"}
          });
      eng = read_io.Open(fname, adios2::Mode::Read);
      std::cerr << rank << ": " << name << " engine created\n";
    }
    else{
      std::cerr << rank << ": receive engine already exists \n";
    }
  
    eng.BeginStep();
    adios2::Variable<T> adios_var = read_io.InquireVariable<T>(name);
  
    const auto total_size = adios_var.Shape()[0];
    const auto my_start = (total_size / nprocs) * rank;
    const auto my_count = (total_size / nprocs);
    std::cout << " Reader of rank " << rank << " reading " << my_count
              << " floats starting at element " << my_start << "\n";
  
    const adios2::Dims start{my_start};
    const adios2::Dims count{my_count};
  
    const adios2::Box<adios2::Dims> sel(start, count);
    Array1d<T>* field = new Array1d<T>{total_size, my_count, 
    	my_start};
  
    adios_var.SetSelection(sel);
    eng.Get(adios_var, field->data());
    eng.EndStep();
    return field;
  }
  
  
  
  
  
  /* receive columns (start_col) to (start_col + localW) */
  template<typename T>
  Array2d<T>* receive2d_from_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
    const std::string fname = dir + "/" + name + ".bp";
  
    if(!eng){
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "480"}
          });
      eng = read_io.Open(fname, adios2::Mode::Read);
      std::cerr << rank << ": " << name << " engine created\n";
    }
    else{
      std::cerr << rank << ": receive engine already exists \n";
    }
    eng.BeginStep();
    adios2::Variable<T> adVar = read_io.InquireVariable<T>(name);
  
    const auto ftn_glob_height = adVar.Shape()[0] ; //4
    const auto ftn_glob_width = adVar.Shape()[1]; // 256005
    //fortran to C transpose
    const auto c_glob_height = ftn_glob_width;
    const auto c_glob_width = ftn_glob_height;
  
    GO local_width  =  c_glob_width / nprocs;
    const GO start = rank * local_width;
    if(rank == nprocs - 1) local_width += c_glob_width%nprocs; // 2
  
    fprintf(stderr, "%d 1.0 name %s nprocs %d"
        "c_glob_width %lu c_glob_height %lu local_width %lu start %lu\n",
        rank, name.c_str(), nprocs,
        c_glob_width, c_glob_height, local_width, start);
  
    Array2d<T>* a2d = new Array2d<T>(c_glob_height, c_glob_width,
        c_glob_height, local_width, start);
    const::adios2::Dims my_start({a2d->start_col(), 0});
    assert(a2d->localH() == a2d->globalH());
    const::adios2::Dims my_offset({a2d->localW(), a2d->globalH()});
    const adios2::Box<adios2::Dims> sel(my_start, my_offset);
  
    adVar.SetSelection(sel);
    eng.Get<T>(adVar, a2d->data());
    eng.EndStep();
    std::cerr << rank <<  ": receive " << name << " done \n";
    return a2d;
  }
  
  /* send columns (start_col) to (start_col + localW) */
  template<typename T>
  void send2d_from_C(const Array2d<T>* a2d, const std::string dir,
      const std::string name, adios2::IO &coupling_io,
      adios2::Engine &engine, adios2::Variable<T> &send_id) {
    const::adios2::Dims g_dims({a2d->globalW(), a2d->globalH()});
    const::adios2::Dims g_offset({a2d->start_col(), 0});
    assert(a2d->localH() == a2d->globalH());
    const::adios2::Dims l_dims({a2d->localW(), a2d->globalH()});
  
    const std::string fname = dir + "/" + name + ".bp";
    if (!engine){
      coupling_io.SetEngine("Sst");
      send_id = coupling_io.DefineVariable<T>(name,
          g_dims, g_offset, l_dims);
      engine = coupling_io.Open(fname, adios2::Mode::Write);
    }
    else{
      std::cerr << ": field engine already created \n";
    }
  
    engine.BeginStep();
    engine.Put<T>(send_id, a2d->data());
    engine.EndStep();
  }
  
  /** Receive PreProc values from GENE
   */
  template<typename T>
  Array1d<T>* receive_gene_pproc(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &engine) {
    const std::string name = "gene_pproc";
    return receive1d_from_ftn<T>(cce_folder,name, io, engine);
  }
  

  Array2d<double>* receive_density(const std::string cce_folder,
  		    adios2::IO &io, adios2::Engine &engine);

  void send_density(const std::string cce_folder, const Array2d<double>* density,
      adios2::IO &io, adios2::Engine &engine, adios2::Variable<double> &send_id);

  Array2d<double>* receive_field(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &eng);

  void send_field(const std::string cce_folder, const Array2d<double>* field,
      adios2::IO &io, adios2::Engine &engine, adios2::Variable<double> &send_id); 

  /** Close the Adios2 engine objects
   */
  void close_engines(adios2::Engine engine[], const int i);
}//end namespace coupler

#endif
