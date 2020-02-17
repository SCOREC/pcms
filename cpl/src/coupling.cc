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

bool first_step=true;

int g_width = 0;
int g_height = 0;

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

  std::cerr << rank <<  ": g_height: "<< g_height << " g_width: "<< g_width <<" \n";
  receive_density(rank, nprocs);
  std::cerr << rank <<  ": receive density done 2.0 \n";

  std::cerr << rank <<  ": g_height: "<< g_height << " g_width: "<< g_width <<" \n";
  send_density(rank, nprocs);
  std::cerr << rank <<  ": send density done 2.1 \n";
  delete[] dens_ptr;

  std::cerr << rank <<  ": g_height: "<< g_height << " g_width: "<< g_width <<" \n";
  receive_field(rank, nprocs);
  std::cerr << rank <<  ": receive field done 2.2 \n";

  std::cerr << rank <<  ": g_height: "<< g_height << " g_width: "<< g_width <<" \n";
  send_field(rank, nprocs);
  std::cerr << rank <<  ": send field done 10.0 \n";

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
  cce_folder = "/global/homes/d/damilare";
}

void receive_density(int rank, int nprocs)
{
  std::string fld_name = "gene_density"; // or data_from_gene??

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
  //This is temporarily for debugging #ADA
  g_width = height;
  g_height = width;
  const::adios2::Dims my_start({start, 0}); //for DebugON
  const::adios2::Dims my_count({count, height}); //for DebugON
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  dens_ptr = new double[height * count]; //contiguously allocate this on the heap from free-list

  dens_id.SetSelection(sel);
  engine.Get<double>(dens_id, dens_ptr);
  engine.EndStep();

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
    std::cerr << 1.30 << std::endl;
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
    std::cerr << 1.40 << std::endl;
    //    bar[1] = &density[0];
    std::cerr << 1.41 << std::endl;
  }
  engine.Close(); // this is done at an external step
}

void send_density(int rank, int nprocs)
{
  int count, start;
  std::string fld_name = "cpl_density";

  std::cerr << "5.0" << std::endl;
  std::cerr << " g_height:" << g_height << " g_width: "<< g_width << std::endl;
  count = g_height / nprocs ; // break the height up
  if(rank == nprocs - 1) count += g_height%nprocs;
  start = rank * count;
  std::cerr << "5.1" << std::endl;
  std::cerr << " start: " << start << " count: "<< count << std::endl;

  const::adios2::Dims g_dims({g_height, g_width});
  const::adios2::Dims g_offset({start, 0});
  const::adios2::Dims l_dims({count, g_width});
  std::cerr << "5.2" << std::endl;

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO write_io = adios.DeclareIO("cpl_density");
  std::cerr << "5.3" << std::endl;
  write_io.SetEngine("Sst");
  std::cerr << "5.4" << std::endl;

  adios2::Variable<double> varid = write_io.DefineVariable<double>(fld_name, g_dims, g_offset, l_dims);
  if(varid) std::cerr << " valid varID" << std::endl;
  std::cerr << "5.5" << std::endl;
  adios2::Engine write_engine = write_io.Open(cce_folder + "/cpl_density.bp", adios2::Mode::Write);
  if(!write_engine) std::cerr << " NULL pointer for writer "<< std::endl;
  std::cerr << "5.6" << std::endl;

  write_engine.BeginStep();
  std::cerr << "5.7" << std::endl;
  write_engine.Put<double>(varid, dens_ptr);
  std::cerr << "5.8" << std::endl;
  write_engine.EndStep();
  std::cerr << "5.9" << std::endl;
  write_engine.Close();
}


void receive_field(int rank, int nprocs)
{
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO read_io = adios.DeclareIO("xgc_field");
  read_io.SetEngine("Sst");
  read_io.SetParameters({{"DataTransport","RDMA"},  {"OpenTimeoutSecs", "360"}});
  std::cerr << "4.2" << std::endl;

  adios2::Engine read_engine = read_io.Open(cce_folder + "/xgc_field.bp", adios2::Mode::Read);
  std::cout << "XGC-to-coupling field engine created\n";

  read_engine.BeginStep();
  std::cerr << "4.3" << std::endl;
  adios2::Variable<double> field_id = read_io.InquireVariable<double>("xgc_field");
  std::cerr << "4.4" << std::endl;
  auto height = field_id.Shape()[0] ; //4
  auto width = field_id.Shape()[1]; // 256005

  //int count  =  width / nprocs;
  int count  =  height / nprocs;
  const int start = rank * count;
//  if(rank == nprocs - 1) count += width%nprocs; // 128002 //128003

  fprintf(stderr, "%d 1.0 nprocs %d width %d height %d count %d start %d\n",
      rank, nprocs, width, height, count, start);
  //This is temporarily for debugging #ADA
  g_width = height;// 4
  g_height = width;// 128002
  const::adios2::Dims my_start({start, 0}); //for DebugON
  std::cerr << "4.5" << std::endl;
  //const::adios2::Dims my_count({count, height}); //for DebugON
  const::adios2::Dims my_count({count, width}); //for DebugON
  std::cerr << "4.6" << std::endl;
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  std::cerr << "4.7" << std::endl;
//  field_ptr = new double[height * count]; 
  field_ptr = new double[width * count]; 
  std::cerr << "4.8" << std::endl;

  field_id.SetSelection(sel);
  std::cerr << "4.9" << std::endl;
  read_engine.Get<double>(field_id, field_ptr);
  std::cerr << "4.9.9" << std::endl;
  read_engine.EndStep();
  read_engine.Close();
}


void send_field(int rank, int nprocs)
{
  std::string fld_name = "cpl_field";

  std::cerr << "9.0"<< std::endl;
  int count  = g_height / nprocs;
  if(rank == nprocs - 1) count += g_height%nprocs;
  int start = rank * count;

  std::cerr << "9.1"<< std::endl;
  const::adios2::Dims g_dims({g_height, g_width});
  const::adios2::Dims g_offset({start, 0});
  const::adios2::Dims l_dims({count, g_width});
  std::cerr << "9.2"<< std::endl;

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  std::cerr << "9.3"<< std::endl;
  adios2::IO coupling_io = adios.DeclareIO("cpl_field");
  std::cerr << "9.4"<< std::endl;
  coupling_io.SetEngine("Sst");

  auto send_id = coupling_io.DefineVariable<double>(fld_name, g_dims, g_offset, l_dims);
  std::cerr << "9.5"<< std::endl;
  adios2::Engine engine = coupling_io.Open(cce_folder + "/cpl_field.bp", adios2::Mode::Write);
  std::cerr << "9.6"<< std::endl;

  engine.BeginStep();
  std::cerr << "9.7"<< std::endl;
  engine.Put<double>(send_id, field_ptr);
  std::cerr << "9.8"<< std::endl;
  engine.EndStep();
  engine.Close();
  std::cerr << "9.9"<< std::endl;
}


