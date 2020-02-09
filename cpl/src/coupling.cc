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
using twod_vec = std::vector<std::vector<double>>;
twod_vec dens_holder = {{0.0},{0.0}};
twod_vec dens_sender = {{0.0},{0.0}};

//global functions;
void initialize_coupling();
void finalize_coupling();

//functions to trim a string
static inline void ltrim(std::string &s);
static inline void rtrim(std::string &s);
static inline void trim(std::string &s);

void receive_density(double* &vec, int rank, int nprocs);
void send_density(const twod_vec &density, int rank, int size);

void receive_field(twod_vec &data_block, int rank, int size);
void send_field(const twod_vec &field, std::vector<double> pot0, int flag);


int main(int argc, char **argv){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  twod_vec dens_vec = {{0.0},{0.0}};

  initialize_coupling();

  double * dens_arr = NULL;
  receive_density(dens_arr, rank, size);
  std::cerr << rank <<  " 2.0 \n";

  /*  send_density(dens_vec, rank, size);
      std::cerr << rank <<  " 3.0 \n";
      */
  /*  if(!rank)
      {
      for (int i = 0; i < 10; i++)
      {
      std::cerr << rank <<  " The values of the first 10 at "<< i <<" is "<< (double)dens_arr[i] <<"\n";
      }

      for (int i = 0; i < 10; i++)
      {
      std::cerr << rank <<  " The values of the last 10 at "<< (183529 * 32)-9+i <<" is "<< (double)dens_arr[(32*183529)-9+i] <<"\n";
      }
      }*/
  //    send field data to GENE

  //  end loop
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

void finalize_coupling()
{

}

static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
        }));
}

static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
        }).base(), s.end());
}

static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

void receive_density(double * &foo, int rank, int nprocs)
{
  std::string fld_name = "gene_density"; // or data_from_gene??

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugOFF);
  adios2::IO dens_io = adios.DeclareIO("density_coupling");
  dens_io.SetEngine("Sst");
  dens_io.SetParameters({{"DataTransport","RDMA"},  {"OpenTimeoutSecs", "360"}});

  trim(cce_folder);
  auto p =  dens_io.Parameters();
  //print out the set parameters for this IO
  if (rank == 0)
  {
    time_t itime = time(NULL);// print wall time
    std::cout << "The wall clock time at Open is "<< ctime(&itime) << std::endl;
    for(auto it = p.begin(); it != p.end(); ++it)
    {
      std::cout<< it->first << " second: " << it->second <<std::endl;
    }
  }

  adios2::Engine engine = dens_io.Open("/global/homes/d/damilare/density.bp", adios2::Mode::Read);
  fprintf(stderr,"GENE-to-coupling density engine created by %d\n", rank);
  fprintf(stderr,"%d 0.6\n", rank);

  engine.BeginStep();
  std::cerr << rank <<  " 0.7 \n";
  adios2::Variable<double> dens_id = dens_io.InquireVariable<double>(fld_name);
  auto height = dens_id.Shape()[0];
  auto width = dens_id.Shape()[1];
  std::cerr << rank <<  " 0.9\n";

  //int count  =  height / nprocs;
  int count  =  width / nprocs;
  if(rank == nprocs - 1) count += width%nprocs;
  //if(rank == nprocs - 1) count += height%nprocs;
  const int start = rank * count;

  fprintf(stderr, "%d 1.0 nprocs %d width %d height %d count %d start %d\n",
      rank, nprocs, width, height, count, start);
  const::adios2::Dims my_start({0, start});
  std::cerr << rank <<  " 1.1 \n";
  const::adios2::Dims my_count({height, count});
  //const::adios2::Dims my_count({width, count});
  std::cerr << rank <<  " 1.2 \n";
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  std::cerr << rank <<  " 1.3 \n";
  //foo = new double[width * count];
  foo = new double[height * count];

  std::cerr << rank <<  " 1.41 \n";
  dens_id.SetSelection(sel);
  std::cerr << rank <<  " 1.5 \n";
  engine.Get<double>(dens_id, foo);
  std::cerr << rank <<  " 1.6 \n";
  engine.EndStep();

  // confirming the values sent
  for (int i = 0; i < 10; i++)
  {
    // the first 10 values for each rank
    std::cerr << rank <<  ": first 10 density at "<< i << " is "<< (double)foo[i] <<"\n";
    // for process zero, get the last values at row 67325 -1
    if(!rank)
    {
      int last_ten = ((67235+1) * height) - 10;
      std::cerr << rank <<  ": last 10 at " <<  last_ten + i << " is "<< (double)foo[last_ten + i] <<"\n";
      int next_ten = ((67235) * height) - 10;
      std::cerr << rank <<  ": previous 10 at " <<  next_ten + i << " is "<< (double)foo[next_ten + i] <<"\n";
    }
  }

//  if(!rank)
  {
	  std::cerr << rank << ": LAST ENTRY: 2,151,552 is " << foo[2151552] << "\n";
  }

  // for rank 1, get the last 10
  int last_ten = (height * count) - 10;
  if(rank == 1)
  {
    for (int i = 0; i < 10; i++)
    {
      std::cerr << rank <<  ": last 10 density at " << last_ten + i << " is "<< (double)foo[last_ten + i] <<"\n";
    }
    std::cerr << rank <<  " 1.7 \n";
    engine.Close();
    std::cerr << rank <<  " 1.8 \n";
    //dens_holder = dens;
  }
}

void send_density(const twod_vec &dens, int rank, int size)
{
  int Nx, Ny, count_x, start_x;
  std::string fld_name = "cpl_density";

  Nx = dens[0].size();  //node_number
  Ny = dens.size();  //nphi
  count_x = Nx / size;
  if(rank == size - 1) count_x += Nx%size;
  start_x = rank * count_x;

  const::adios2::Dims gdims({Ny, Nx});
  const::adios2::Dims goffset({0, start_x});
  const::adios2::Dims ldims({Ny, count_x});

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO read_io = adios.DeclareIO("density_from_coupling");
  read_io.SetEngine("Sst");

  auto varid = read_io.DefineVariable<double>(fld_name, gdims, goffset, ldims);
  trim(cce_folder);
  adios2::Engine read_engine = read_io.Open(cce_folder + "/cpl_density.bp", adios2::Mode::Write);

  read_engine.BeginStep();
  read_engine.Put<double>(varid, (dens.data())->data() );
  read_engine.EndStep();
  read_engine.Close();
  dens_sender = dens;
}


void receive_field(twod_vec &data_block, int rank, int size)
{
  int count, start;
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO send_io = adios.DeclareIO("field_coupling");
  send_io.SetEngine("Sst");
  trim(cce_folder);

  adios2::Engine send_engine = send_io.Open((cce_folder + "/field.bp"), adios2::Mode::Read);
  std::cout << "XGC-to-coupling field engine created\n";
  send_engine.BeginStep();
  adios2::Variable<double> field_id = send_io.InquireVariable<double>("dadat");
  auto height = field_id.Shape()[0] ;
  auto width = field_id.Shape()[1];
  std::cout << "first dim - Incoming variable is of size " << height << std::endl;
  std::cout << "second dim - Incoming variable is of size " << width << std::endl;

  count  = (width / size);
  if(rank == size - 1) count += width%size;
  start = rank * count;

  std::cout << "cField Reader of rank " << rank << " reading " << count
    << " floats starting at element " << start << "\n";

  const adios2::Dims my_start({0, start});
  const adios2::Dims my_count({height, count});
  const adios2::Box<adios2::Dims> sel(my_start, my_count);
  // resize the block to fit 
  for(int i=0; i < height; i++)
  {
    data_block[i].resize(count);// accessing wrong index and sizing wrong, fix
  }

  field_id.SetSelection(sel);
  send_engine.Get<double>(field_id, (data_block.data())->data());
  send_engine.EndStep();
  send_engine.Close();
}


void send_field(const twod_vec &field, int rank, int size)
{
  const std::size_t width = field[0].size();
  const std::size_t height = field.size();

  std::string fld_name = "cpl_field";

  const::adios2::Dims gdims({height, width});
  int  count  = (width / size);
  if(rank == size - 1) count += width%size;
  const::adios2::Dims goffset({0, rank * count});
  const::adios2::Dims ldims({height, count});

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO coupling_io = adios.DeclareIO("field_from_coupling");
  coupling_io.SetEngine("Sst");

  auto field_id = coupling_io.DefineVariable<double>(fld_name, gdims, goffset, ldims);
  trim(cce_folder);
  adios2::Engine engine = coupling_io.Open(cce_folder + "/cpl_field.bp", adios2::Mode::Write);

  engine.BeginStep();
  engine.Put<double>(field_id, (field.data())->data());
  engine.EndStep();
  engine.Close();
}


