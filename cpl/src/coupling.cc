#include<adios2.h>
#include<iostream>
#include<mpi.h>
#include<cassert>
//#include <Kokkos_Core.hpp>
#include <typeinfo>

typedef long unsigned GO;

class Array2d {
  public:
    Array2d(GO gH, GO gW, GO lH, GO lW, GO start) :
      globH(gH), globW(gW), locH(lH), locW(lW), locFirstCol(start) {
        vals = new double[locH*locW];
    }
    ~Array2d() {
      globH = globW = locH = locW = 0;
      delete [] vals;
    }
    double val(long i) const {
      assert(i<(locH*locW));
      return vals[i];
    }
    double* data() const { return vals; };
    GO globalH() const { return globH; };
    GO globalW() const { return globW; };
    GO localH() const { return locH; };
    GO localW() const { return locW; };
    GO start_col() const { return locFirstCol; };
  private:
    double* vals;
    GO globH;
    GO globW;
    GO locH;
    GO locW;
    GO locFirstCol;
};

void printSomeDensityVals(const Array2d* density) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //asserting the density values received from GENE
  if(!rank)
  {
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank <<  ": first 10 density at "<< i << " is "<< density->val(i) <<"\n";
    }
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank << ": first 10 for rank 1 at: [67236]" << " + "<< i << " is " << density->val(67236 + i) << "\n";
    }
  }

  if(rank == 1)
  {
    for (int i = 0; i < 10; i++)
    {
      int offset = ((density->localW() - 1) * density->localH()) + 67235 - 9; //width
      std::cerr << rank << ": last 10 for rank 0 at: [67235 - 9]" << " + "<< i << " is " << density->val(offset  + i) << "\n";
    }
    int last_ten = (density->localH() * density->localW()) - 10;
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank <<  ": last 10 density at " << last_ten + i << " is "<< density->val(last_ten + i) <<"\n";
    }
  }
}

/* receive columns (start_col) to (start_col + localW) */
Array2d* receive2d_from_ftn(const std::string dir, const std::string name, adios2::IO &read_io, adios2::Engine &eng) {
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  const std::string fname = dir + "/" + name + ".bp";

  if(!eng){
    read_io.SetEngine("Sst");
    read_io.SetParameters({{"DataTransport","RDMA"},  {"OpenTimeoutSecs", "480"}});
    eng = read_io.Open(fname, adios2::Mode::Read);
    std::cerr << rank << ": " << name << " engine created\n";
  }
  else{
    std::cerr << rank << ": receive engine already exists \n";
  }
  eng.BeginStep();
  adios2::Variable<double> adVar = read_io.InquireVariable<double>(name);

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

  Array2d* a2d = new Array2d(c_glob_height, c_glob_width,
      c_glob_height, local_width, start);
  const::adios2::Dims my_start({a2d->start_col(), 0});
  assert(a2d->localH() == a2d->globalH());
  const::adios2::Dims my_offset({a2d->localW(), a2d->globalH()});
  const adios2::Box<adios2::Dims> sel(my_start, my_offset);

  adVar.SetSelection(sel);
  eng.Get<double>(adVar, a2d->data());
  eng.EndStep();
  std::cerr << rank <<  ": receive " << name << " done \n";
  return a2d;
}

/* send columns (start_col) to (start_col + localW) */
void send2d_from_C(const Array2d* a2d, const std::string dir, const std::string name, adios2::IO &coupling_io, adios2::Engine &engine, adios2::Variable<double> &send_id) {
  const::adios2::Dims g_dims({a2d->globalW(), a2d->globalH()});
  const::adios2::Dims g_offset({a2d->start_col(), 0});
  assert(a2d->localH() == a2d->globalH());
  const::adios2::Dims l_dims({a2d->localW(), a2d->globalH()});

  const std::string fname = dir + "/" + name + ".bp";
  if (!engine){
    coupling_io.SetEngine("Sst");
    send_id = coupling_io.DefineVariable<double>(name, g_dims, g_offset, l_dims);
    engine = coupling_io.Open(fname, adios2::Mode::Write);
  }
  else{
    std::cerr << ": field engine already created \n";
  }

  engine.BeginStep();
  engine.Put<double>(send_id, a2d->data());
  engine.EndStep();
}

Array2d* receive_density(const std::string cce_folder, adios2::IO &io, adios2::Engine &engine)
{
  const std::string name = "gene_density";
  return receive2d_from_ftn(cce_folder,name, io, engine);
}

void send_density(const std::string cce_folder, const Array2d* density,  adios2::IO &io, adios2::Engine &engine, adios2::Variable<double> &send_id)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string fld_name = "cpl_density";
  send2d_from_C(density, cce_folder, fld_name, io, engine, send_id);
  std::cerr << rank <<  ": send " << fld_name <<" done \n";
}

Array2d* receive_field(const std::string cce_folder,  adios2::IO &io, adios2::Engine &eng)
{
  const std::string name = "xgc_field";
  return receive2d_from_ftn(cce_folder,name, io, eng);
}

void send_field(const std::string cce_folder, const Array2d* field, adios2::IO &io, adios2::Engine &engine, adios2::Variable<double> &send_id)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string fld_name = "cpl_field";
  send2d_from_C(field, cce_folder, fld_name, io, engine, send_id);
  std::cerr << rank <<  ": send " << fld_name <<" done \n";
}

void close_engines(adios2::Engine engine[]){
  for(int i = 0; i < 4; i++)
  {
    engine[i].Close();
  }
}

/*
void exParFor() {
  Kokkos::parallel_for(
      4, KOKKOS_LAMBDA(const int i) {
        printf("Hello from kokkos thread i = %i\n", i);
      });
}
*/

int main(int argc, char **argv){
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//  Kokkos::initialize(argc, argv);
/*  if(!rank) {
    printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());
    exParFor();
  }
*/
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO IO[4];
  adios2::Engine engines[4];
  adios2::Variable<double> send_var[2];
  const std::string cce_folder = "../coupling";
  const int time_step = 1, RK_count = 4;

  IO[0] = adios.DeclareIO("gene_density");
  IO[1] = adios.DeclareIO("cpl_density");
  IO[2] = adios.DeclareIO("xgc_field");
  IO[3] = adios.DeclareIO("cpl_field");
  
  for (int i = 0; i < time_step; i++)
  {
    for (int j = 0; j < RK_count; j++)
    {
      Array2d* density = receive_density(cce_folder, IO[0], engines[0]);
      printSomeDensityVals(density);
      send_density(cce_folder, density, IO[1], engines[1], send_var[0]);
      delete density;

      Array2d* field = receive_field(cce_folder, IO[2], engines[2]);
      send_field(cce_folder, field, IO[3], engines[3], send_var[1]);
      delete field;
    }
  }

  close_engines(engines);
//  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

