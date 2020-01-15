#include<iostream>
#include<vector>
#include<assert.h>
#include<adios2.h>
#include<mpi.h>

void read_density(std::vector<float> &density, int rank, int size);
void write_field(const std::vector<float> &field, int rank, int size);

int main(int argc, char *argv[])
{
  int rank , size;
  std::vector<float> density;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);   

  read_density(density, rank, size);
  assert (sizeof(density[0])==4);
  printf ("The first density value is %f\n",density[0]);

  // returning the same varaible density across the loop
  write_field(density, rank, size);
  MPI_Barrier(MPI_COMM_WORLD);
  if( !rank )
    std::cout << "xgc proxy done\n";
  MPI_Finalize();
  return 0;
}


void write_field(const std::vector<float> &field, int rank, int size)
{
  const std::size_t Nx = field.size();

  try
  {
	std::cout << rank << " start writing \n";
  	adios2::ADIOS xc_adios(MPI_COMM_WORLD, adios2::DebugON);
  	adios2::IO xfieldIO = xc_adios.DeclareIO("xcIO");
  	xfieldIO.SetEngine("Sst");

  	auto bp_xfield = xfieldIO.DefineVariable<float>("bp_xfield", {size * Nx}, {rank * Nx}, {Nx});

  	adios2::Engine xfieldWriter = xfieldIO.Open("xfield.bp", adios2::Mode::Write);

  	xfieldWriter.BeginStep();
  	xfieldWriter.Put<float>(bp_xfield, field.data());
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

void read_density(std::vector<float> &dens, int rank, int size)
{

  try
  {
	std::cout << rank << " start reading\n";
  	adios2::ADIOS xc_adios(MPI_COMM_WORLD, adios2::DebugON);

  	adios2::IO cdensIO = xc_adios.DeclareIO("xcIO");
  	cdensIO.SetEngine("Sst");

  	adios2::Engine cdensReader = cdensIO.Open("cdens.bp", adios2::Mode::Read);
  	cdensReader.BeginStep();
  	adios2::Variable<float> bp_cdens =
  		cdensIO.InquireVariable<float>("bp_cdens");
  	std::cout << "Incoming variable is of size " << bp_cdens.Shape()[0]
  		<< "\n";
  	const std::size_t total_size = bp_cdens.Shape()[0];
  	const std::size_t my_start = (total_size / size) * rank;
  	const std::size_t my_count = (total_size / size);
  	std::cout << "cDensity Reader of rank " << rank << " reading " << my_count
  		<< " floats starting at element " << my_start << "\n";

  	const adios2::Dims start{my_start};
  	const adios2::Dims count{my_count};

  	const adios2::Box<adios2::Dims> sel(start, count);

  	dens.resize(my_count);

  	bp_cdens.SetSelection(sel);
  	cdensReader.Get(bp_cdens, dens.data());
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
