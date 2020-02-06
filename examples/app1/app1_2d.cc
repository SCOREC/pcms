#include<iostream>
#include<vector>
#include<assert.h>
#include<adios2.h>
#include<mpi.h>
#include<numeric>

using twoD_vec = std::vector<std::vector<float>>;


void write_density(const twoD_vec &density, int rank, int size);
void read_field(twoD_vec &dens, int rank, int size);

int main(int argc, char *argv[])
{
  int rank, size, step = 0;
  twoD_vec field = {{0.0},{0.0}};

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);   

  std::vector<float> density = {
    (float)10.0 * rank + 0, (float)10.0 * rank + 1, (float)10.0 * rank + 2,
    (float)10.0 * rank + 3, (float)10.0 * rank + 4, (float)10.0 * rank + 5,
    (float)10.0 * rank + 6, (float)10.0 * rank + 7, (float)10.0 * rank + 8,
    (float)10.0 * rank + 9};

  twoD_vec dens = {density, density};

  while(step <5)
  {
    write_density(dens, rank, size);
    read_field(field, rank, size);
    std::cout << "This is for step "<< step <<std::endl;
    step++;
  }


  MPI_Barrier(MPI_COMM_WORLD);
  if ( !rank )
  {
    assert (field[0][0] == (float)0);
    std::cout << "The asserted first field value is " << field[0][0] << "\n";
    std::cout << "gene proxy done\n";
  }
  MPI_Finalize();
  return 0;
}


void write_density(const twoD_vec &density, int rank, int size)
{
  const std::size_t Nx = density[0].size();
  const std::size_t Ny = density.size();

  try
  {
    adios2::ADIOS gc_adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO gdensIO = gc_adios.DeclareIO("gcIO");
    gdensIO.SetEngine("Sst");
    adios2::Dims gdims, start, count;
    gdims = {Ny, size * Nx};
    start = {0, rank * Nx};
    count = {Ny, Nx};

    auto bp_gdens = gdensIO.DefineVariable<float>("bp_gdens", gdims, start, count);

    adios2::Engine gdensWriter = gdensIO.Open("gdens.bp", adios2::Mode::Write);

    gdensWriter.BeginStep();
    gdensWriter.Put<float>(bp_gdens, (density.data())->data());
    gdensWriter.EndStep();
    gdensWriter.Close();
  }

  catch(std::invalid_argument &e)
  {
    std::cout << "Invalid argumant exception, STOPPING PROGRAM from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch(std::ios_base::failure &e)
  {
    std::cout << "IO system base failure exception, STOPPING PROGRAM from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch(std::exception &e)
  {
    std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
    std::cout << e.what() << "\n";
  }
}

void read_field(twoD_vec &field, int rank, int size)
{
  try
  {
    adios2::ADIOS cg_adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO cfieldIO = cg_adios.DeclareIO("cgIO");
    cfieldIO.SetEngine("Sst");

    adios2::Engine cfieldReader = cfieldIO.Open("xfield.bp", adios2::Mode::Read);
    cfieldReader.BeginStep();

    adios2::Variable<float> bp_cfield = cfieldIO.InquireVariable<float>("bp_xfield");
    auto height = bp_cfield.Shape()[0] ;
    auto width = bp_cfield.Shape()[1];
    std::cout << "first dim - Incoming variable is of size " << height << std::endl;
    std::cout << "second dim - Incoming variable is of size " << width << std::endl;

    auto shape =  bp_cfield.Shape();
    size_t data_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    std::cout << "shapes accumulated gives Incoming variable is of size " << data_size <<std::endl;
    long unsigned int f_count  = (width/size);
    long unsigned int f_start = rank * f_count;
//    if(rank == size - 1) f_count += width%size;

    std::cout << "cField Reader of rank " << rank << " reading " << f_count
      << " floats starting at element " << f_start << "\n";

    const adios2::Dims start{0, f_start};
    const adios2::Dims count{height, f_count};
    const adios2::Box<adios2::Dims> sel(start, count);

    for(int i=0; i < height; i++)
    {
      field[i].resize(width);
    }

    bp_cfield.SetSelection(sel);
    cfieldReader.Get(bp_cfield, (field.data())->data());
    cfieldReader.EndStep();
    cfieldReader.Close();
  }
  catch (std::invalid_argument &e)
  {
    std::cout << "Invalid argument exception, STOPPING PROGRAM from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch (std::ios_base::failure &e)
  {
    std::cout << "IO System base failure exception, STOPPING PROGRAM "
      "from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch (std::exception &e)
  {
    std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
    std::cout << e.what() << "\n";
  }
}
