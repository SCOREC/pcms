#include<adios2.h>
#include<fstream>
#include<iostream>
#include <mpi.h>
#include<vector>
#include <algorithm> 
#include<numeric>
#include<time.h>

//global variables
bool cce_dpot_index0 ;

int cce_density_model, cce_step, cce_field_step, cce_last_node, cce_comm_density_mode;
int cce_first_surface_coupling, cce_last_surface_coupling, cce_field_model, cce_comm_field_mode;
int cce_npsi, cce_dt, cce_side, cce_alpha, cce_all_surface_number, itime, mype;

unsigned long cce_field_node_number, sml_nphi_total, sml_intpl_mype, cce_node_number, cce_first_node;
std::string cce_folder;

int g_width = 0;
int g_height = 0;
int timestep = 0;


double *dens_ptr = NULL;
double *field_ptr = NULL;

void initialize_coupling();
void finalize_coupling();

void receive_density(int rank, int nprocs);
void send_density(int rank, int nprocs);
void receive_field(int rank, int nprocs);
void send_field(int rank, int flag);


int main(int argc, char **argv){
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  initialize_coupling();

  receive_density(rank, nprocs);

  send_density(rank, nprocs);
 // delete[] dens_ptr;

  receive_field(rank, nprocs);
  send_density(rank, nprocs);

//  send_field(rank, nprocs);

  delete[] field_ptr;
  return 0;

  MPI_Finalize();
}

void initialize_coupling()
{
  cce_alpha = 0.5;
  cce_density_model = 0;
  cce_step = 0; //In case of restart, one might want to change this
  cce_field_step = 0;
  cce_first_node = 10;
  cce_last_node = 0;
  cce_comm_density_mode = 2;
  cce_all_surface_number = 1;
  cce_first_surface_coupling = -1;
  cce_last_surface_coupling = -1;
  cce_field_model = 0;
  cce_comm_field_mode = 0;
  cce_npsi = -1;
  cce_dpot_index0=false;
  cce_dt = -1; 

  // attempt to read the remaining values from 'coupling.in' file - cce_side, cce_first_node, cce_last_node, cce_folder
  std::map<std::string, int> input_map;
  std::map<std::string, int>::iterator iter;

  std::fstream in;
  in.open("coupling.in");
  if (!in)
  {
    std::cout << " Error, could not open the coupling.in file. \n";
  }

  int value;
  std::string key;

  while (in >> key >> value)
  {  
    getline(in, key, '=');
    input_map[key] = value;
  }
  iter = input_map.find("cce_last_node");
  if(iter != input_map.end())
  {
    std::cout << "cce_last_node is = " << iter->second << "\n";
    cce_last_node = iter->second;
  }

  iter = input_map.find("cce_last_node");
  if(iter != input_map.end())
  {
    std::cout << "cce_last_node is = " << iter->second << "\n";
    cce_first_node = iter->second;
  }
  cce_node_number = cce_last_node - cce_first_node + 1;
  cce_node_number = 212817 - 1875 + 1;
  cce_side = 3;
  cce_folder = "../coupling";
}

void receive_density(int rank, int nprocs)
{
  std::string fld_name = "gene_density";

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO dens_io = adios.DeclareIO("gene_density");
  dens_io.SetEngine("Sst");
  dens_io.SetParameters({{"DataTransport","RDMA"},  {"OpenTimeoutSecs", "360"}});

  adios2::Engine engine = dens_io.Open(cce_folder + "/gene_density.bp", adios2::Mode::Read);
  fprintf(stderr,"GENE-to-coupling density engine created by %d\n", rank);

  engine.BeginStep();
  adios2::Variable<double> dens_id = dens_io.InquireVariable<double>(fld_name);
  auto width = dens_id.Shape()[0]; // 32 
  auto height = dens_id.Shape()[1];// 183529

  int count  =  width / nprocs;
  if(rank == nprocs - 1) count += width%nprocs; // 16
  const int start = rank * count;

  fprintf(stderr, "%d 1.0 nprocs %d width %d height %d count %d start %d\n",
      rank, nprocs, width, height, count, start);

  g_width = height;
  g_height = width;

  const::adios2::Dims my_start({start, 0}); 
  const::adios2::Dims my_count({count, height}); 
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  dens_ptr = new double[height * count]; 

  dens_id.SetSelection(sel);
  engine.Get<double>(dens_id, dens_ptr);
  engine.EndStep();

  //asserting the density values received from GENE
  if(!rank)
  {
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank <<  ": first 10 density at "<< i << " is "<< dens_ptr[i] <<"\n";
    }
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank << ": first 10 for rank 1 at: [67236]" << " + "<< i << " is " << dens_ptr[67236 + i] << "\n";
    }
  }

  if(rank == 1)
  {
    for (int i = 0; i < 10; i++)
    {
      int offset = ((count - 1) * height) + 67235 - 9; //width 
      std::cerr << rank << ": last 10 for rank 0 at: [67235 - 9]" << " + "<< i << " is " << dens_ptr[offset  + i] << "\n";
    }
    int last_ten = (height * count) - 10;
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank <<  ": last 10 density at " << last_ten + i << " is "<< dens_ptr[last_ten + i] <<"\n";
    }
  }
  engine.Close(); // this is done at an external step
  std::cerr << rank <<  ": receive density done \n";
}

void send_density(int rank, int nprocs)
{
  int count, start;
  adios2::Variable<double> varid;
  adios2::Engine write_engine; 
  std::string fld_name = "cpl_density";

  //if(!timestep)
//  {
    count = g_height / nprocs ; 
    if(rank == nprocs - 1) count += g_height%nprocs;
    start = rank * count;

    const::adios2::Dims g_dims({g_height, g_width});
    const::adios2::Dims g_offset({start, 0});
    const::adios2::Dims l_dims({count, g_width});
    adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO write_io = adios.DeclareIO("cpl_density");
    write_io.SetEngine("Sst");

   varid = write_io.DefineVariable<double>(fld_name, g_dims, g_offset, l_dims);
   write_engine = write_io.Open(cce_folder + "/cpl_density.bp", adios2::Mode::Write);
  //}
  timestep++;
  write_engine.BeginStep();
  write_engine.Put<double>(varid, dens_ptr);
  write_engine.EndStep();
  //write_engine.Close();
  //if(timestep) write_engine.Close();
  std::cerr << rank <<  ": send density done \n";
}


void receive_field(int rank, int nprocs)
{
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO read_io = adios.DeclareIO("xgc_field");
  read_io.SetEngine("Sst");
  read_io.SetParameters({{"DataTransport","RDMA"},  {"OpenTimeoutSecs", "360"}});

  adios2::Engine read_engine = read_io.Open(cce_folder + "/xgc_field.bp", adios2::Mode::Read);
  std::cout << "XGC-to-coupling field engine created\n";

  read_engine.BeginStep();
  adios2::Variable<double> field_id = read_io.InquireVariable<double>("xgc_field");
  auto height = field_id.Shape()[0] ; //4
  auto width = field_id.Shape()[1]; // 256005

  int count  =  height / nprocs;
  const int start = rank * count;
  if(rank == nprocs - 1) count += height%nprocs; // 2

  fprintf(stderr, "%d 1.0 nprocs %d width %d height %d count %d start %d\n",
      rank, nprocs, width, height, count, start);

  g_width = height;// 4 updating global height to be used for the field send
  g_height = width;// 128002 updating global width
  const::adios2::Dims my_start({start, 0}); 
  const::adios2::Dims my_count({count, width}); 
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  field_ptr = new double[width * count]; 

  field_id.SetSelection(sel);
  read_engine.Get<double>(field_id, field_ptr);
  read_engine.EndStep();
  //  read_engine.Close();
  std::cerr << rank <<  ": receive field done \n";
}


void send_field(int rank, int nprocs)
{
  std::string fld_name = "cpl_field";

  int count  = g_height / nprocs;
  if(rank == nprocs - 1) count += g_height%nprocs;
  int start = rank * count;

  const::adios2::Dims g_dims({g_height, g_width});
  const::adios2::Dims g_offset({start, 0});
  const::adios2::Dims l_dims({count, g_width});

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO coupling_io = adios.DeclareIO("cpl_field");
  coupling_io.SetEngine("Sst");

  auto send_id = coupling_io.DefineVariable<double>(fld_name, g_dims, g_offset, l_dims);
  adios2::Engine engine = coupling_io.Open(cce_folder + "/cpl_field.bp", adios2::Mode::Write);

  engine.BeginStep();
  engine.Put<double>(send_id, field_ptr);
  engine.EndStep();
  engine.Close();
  std::cerr << rank <<  ": send field done \n";
}


