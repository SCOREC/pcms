#include<adios2.h>
#include<fstream>
#include<iostream>
//#include<map>
#include <mpi.h>
#include<vector>
#include <algorithm> 
//#include <cctype>
//#include <locale>
//#include<string.h>
//#include<thread>
#include<numeric>

//global variables
bool cce_dpot_index0 ;

int cce_density_model, cce_step, cce_field_step, cce_last_node, cce_comm_density_mode;
int cce_first_surface_coupling, cce_last_surface_coupling, cce_field_model, cce_comm_field_mode;
int cce_npsi, cce_dt, cce_side, cce_alpha, cce_all_surface_number, itime, mype;

unsigned long cce_field_node_number, sml_nphi_total, sml_intpl_mype, cce_node_number, cce_first_node;

std::string cce_folder;
using twoD_vec = std::vector<std::vector<double>>;

// cce_density is set as a global variable
std::vector<double> density;

//global functions;
void initialize_coupling();
void finalize_coupling();

//functions to trim a string
static inline void ltrim(std::string &s);
static inline void rtrim(std::string &s);
static inline void trim(std::string &s);

void send_dens_to_XGC(twoD_vec density, int iphi, int nphi, int block_count, int block_start,int block_end, int comm);
void receive_field_from_XGC(twoD_vec &data_block,int block_start,int block_end, int block_count, int nphi);

void send_field_to_GENE(const twoD_vec field, std::vector<double> pot0, int flag);
void receive_dens_from_GENE();


int main(int argc, char **argv){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  initialize_coupling();
  //  call cce_initialize to get the required cce_values
  // initialize adios2
  //  create a while loop

  //    get density data from GENE


  //    CPL-XGC Writer - send density data to XGC in the same format

  //    get field data from XGC


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
  cce_folder = "../coupling";
  // maxplanes gets its values from ncuts thats defined in an hdf5 file  
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


void send_dens_to_XGC(twoD_vec dens, int iphi, unsigned long nphi, unsigned long block_count, unsigned long block_start, unsigned long block_end, int comm)
{
  int maxplane;
  adios2::Dims gdims, goffset, ldims;
  std::string fld_name = "data";

  gdims = {cce_node_number, nphi};
  goffset = {block_start, 0};
  ldims = {block_count, nphi};

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO cplDensIO = adios.DeclareIO("density_from_coupling");
  cplDensIO.SetEngine("Sst");

  auto bpDens = cplDensIO.DefineVariable<double>(fld_name, gdims, goffset, ldims);
  trim(cce_folder);
  adios2::Engine cplDensWriter = cplDensIO.Open(cce_folder + "/cpl_density.bp", adios2::Mode::Write);
  maxplane = nphi - 1;

  cplDensWriter.BeginStep();
  cplDensWriter.Put<double>(bpDens, (dens.data())->data() );
  cplDensWriter.EndStep();
  cplDensWriter.Close();
}


void receive_field_from_XGC(twoD_vec &data_block, int block_start, int block_end, unsigned long block_count, unsigned long nphi)
{
  int maxplane = nphi - 1;
  adios2::Dims start, count;

  start = {(cce_first_node - 1 + block_start), 0};
  count = {block_count, nphi};

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO xgc_fieldIO = adios.DeclareIO("field_from_XGC");
  xgc_fieldIO.SetEngine("Sst");
  trim(cce_folder);

  adios2::Engine xgc_field_reader = xgc_fieldIO.Open((cce_folder + "/field.bp"), adios2::Mode::Read);
  std::cout << "XGC-to-coupling field engine created\n";
  xgc_field_reader.BeginStep();
  adios2::Variable<double> bp_field = xgc_fieldIO.InquireVariable<double>("dadat");

  auto shape =  bp_field.Shape();
  size_t data_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
  data_block.resize(data_size);
  const adios2::Box<adios2::Dims> sel(start, count);

  bp_field.SetSelection(sel);
  xgc_field_reader.Get<double>(bp_field, (data_block.data())->data());
  xgc_field_reader.EndStep();
  xgc_field_reader.Close();
}


void send_field_to_GENE(const std::vector<double> field, int flag)
{
  int maxplane;
  adios2::Dims gdims, goffset, ldims;
  std::string fld_name = "dadat";

  gdims = {cce_field_node_number, sml_nphi_total};
  goffset = {0, sml_intpl_mype};
  ldims = {cce_field_node_number, 1};

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO gene_fieldIO = adios.DeclareIO("field_to_gene");
  gene_fieldIO.SetEngine("Sst");

  auto bp_gene_field = gene_fieldIO.DefineVariable<double>(fld_name, gdims, goffset, ldims);
  trim(cce_folder);
  adios2::Engine gene_field_writer = gene_fieldIO.Open(cce_folder + "/cpl_field.bp", adios2::Mode::Write);

  gene_field_writer.BeginStep();
  gene_field_writer.Put<double>(bp_gene_field, field.data());
  gene_field_writer.EndStep();
  gene_field_writer.Close();

}


void receive_dens_from_GENE()
{
  adios2::Dims start, count;
  std::vector<double> arrtmp(cce_node_number);
  std::string fld_name = "data"; // or data_from_gene??

  start = {0, sml_intpl_mype};
  count = {cce_node_number, 1};

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO gene_densIO = adios.DeclareIO("density_from_GENE");
  gene_densIO.SetEngine("Sst");

  trim(cce_folder);
  adios2::Engine gene_density_reader = gene_densIO.Open((cce_folder + "/density.bp"), adios2::Mode::Read);
  std::cout << "GENE-to-coupling density engine created\n";

  gene_density_reader.BeginStep();
  adios2::Variable<double> bp_density = gene_densIO.InquireVariable<double>(fld_name);
  auto shape = bp_density.Shape();
  size_t data_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
  density.resize(data_size);

  const adios2::Box<adios2::Dims> sel(start, count);
  bp_density.SetSelection(sel);

  gene_density_reader.Get<double>(bp_density, arrtmp.data());
  gene_density_reader.EndStep();
  gene_density_reader.Close();
  density = arrtmp;
}
