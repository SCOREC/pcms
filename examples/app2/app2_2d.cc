#include<iostream>
#include<vector>
#include<assert.h>
#include<adios2.h>
#include<mpi.h>
#include<numeric>

using twoD_vec = std::vector<std::vector<float>>; 
void read_density(twoD_vec &density, int rank, int size);
void write_field(const twoD_vec &field, int rank, int size);

int main(int argc, char *argv[])
{
  int rank , size, step = 0;
  twoD_vec density = {{0.0},{0.0}};

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);   

    while(step < 5)
   {
  read_density(density, rank, size);

  // returning the same varaible density across the loop
  write_field(density, rank, size);
  std::cout << "This is for time step " << step << std::endl;
   step++;
   }
  MPI_Barrier(MPI_COMM_WORLD);
  if( !rank )
  {
    assert (density[0][0] == (float)0);
    std::cout << "The asserted first density value is " << density[0][0] << "\n";
    std::cout << "xgc proxy done\n";
  }
  MPI_Finalize();
  return 0;
}

void write_field(const twoD_vec  &field, int rank, int size)
{
  const std::size_t Nx = field[0].size();
  const std::size_t Ny = field.size();

  try
  {
    std::cout << rank << " start writing \n";
    adios2::ADIOS xc_adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO xfieldIO = xc_adios.DeclareIO("xcIO");
    xfieldIO.SetEngine("Sst");
    adios2::Dims gdims, start, count;

    gdims = {Ny, Nx};
    const std::size_t count_x = Nx/size;
    start = {0, rank * count_x};
    count = {Ny, count_x};

    auto bp_xfield = xfieldIO.DefineVariable<float>("bp_xfield", gdims, start, count);
    adios2::Engine xfieldWriter = xfieldIO.Open("xfield.bp", adios2::Mode::Write);

    xfieldWriter.BeginStep();
    xfieldWriter.Put<float>(bp_xfield, (field.data())->data());
    xfieldWriter.EndStep();
    xfieldWriter.Close();
    std::cout << rank << " Done writing \n";
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

void read_density(twoD_vec &dens, int rank, int size)
{

  try
  {
    std::cout << rank << " start reading\n";
    adios2::ADIOS xc_adios(MPI_COMM_WORLD, adios2::DebugON);

    adios2::IO cdensIO = xc_adios.DeclareIO("xcIO");
    cdensIO.SetEngine("Sst");

    adios2::Engine cdensReader = cdensIO.Open("gdens.bp", adios2::Mode::Read);
    cdensReader.BeginStep();
    adios2::Variable<float> bp_cdens = cdensIO.InquireVariable<float>("bp_gdens");
    auto height = bp_cdens.Shape()[0] ;
    auto width = bp_cdens.Shape()[1];
    std::cout << "first dim - Incoming variable is of size " << height << std::endl;
    std::cout << "second dim - Incoming variable with size " << width << std::endl;

    auto shape =  bp_cdens.Shape();
    size_t data_size = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    std::cout << "shapes accumulated variable has size " << data_size <<std::endl;

    long unsigned int f_count  = (width/size);
    long unsigned int f_start = rank * f_count;
    //if(rank == size - 1) f_count += width%size;


    std::cout << "cDensity Reader of rank " << rank << " reading " << f_count
      << " floats starting at element " << f_start << "\n";

    const adios2::Dims start{0, f_start};
    const adios2::Dims count{height, f_count};
    const adios2::Box<adios2::Dims> sel(start, count);

    for(int i=0; i < height; i++)
    {
      dens[i].resize(width);
    }

    bp_cdens.SetSelection(sel);
    cdensReader.Get(bp_cdens, (dens.data())->data());
    cdensReader.EndStep();

    cdensReader.Close();
    std::cout << rank << " done reading\n";
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
