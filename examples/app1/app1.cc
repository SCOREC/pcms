#include<iostream>
#include<vector>

#include<adios2.h>
#include<mpi.h>

void write_density(const std::vector<float> &density, int rank, int size);
void read_field(std::vector<float> &dens, int rank, int size);

int main(int argc, char *argv[])
{
  int rank , size;
  std::vector<float> field = {0.0};

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);   

  std::vector<float> density = {
    (float)10.0 * rank + 0, (float)10.0 * rank + 1, (float)10.0 * rank + 2,
    (float)10.0 * rank + 3, (float)10.0 * rank + 4, (float)10.0 * rank + 5,
    (float)10.0 * rank + 6, (float)10.0 * rank + 7, (float)10.0 * rank + 8,
    (float)10.0 * rank + 9};

  write_density(density, rank, size);
  read_field(field, rank, size);
 MPI_Finalize();
return 0;
}


void write_density(const std::vector<float> &density, int rank, int size)
{
 const std::size_t Nx = density.size();

  try
  {
    adios2::ADIOS gc_adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO gdensIO = gc_adios.DeclareIO("gcIO");
    gdensIO.SetEngine("Sst");

    auto bp_gdens = gdensIO.DefineVariable<float>("bp_gdens", {size * Nx}, {rank * Nx}, {Nx});

    adios2::Engine gdensWriter = gdensIO.Open("gdens.bp", adios2::Mode::Write);

    gdensWriter.BeginStep();
    gdensWriter.Put<float>(bp_gdens, density.data());
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

void read_field(std::vector<float> &field, int rank, int size)
{
    try
    {
        adios2::ADIOS cg_adios(MPI_COMM_WORLD, adios2::DebugON);
        adios2::IO cfieldIO = cg_adios.DeclareIO("cgIO");
        cfieldIO.SetEngine("Sst");

        adios2::Engine cfieldReader = cfieldIO.Open("cfield.bp", adios2::Mode::Read);
        cfieldReader.BeginStep();

        adios2::Variable<float> bp_cfield = cfieldIO.InquireVariable<float>("bp_cfield");
        std::cout << "Incoming variable is of size " << bp_cfield.Shape()[0] << "\n";

        const std::size_t total_size = bp_cfield.Shape()[0];
        const std::size_t my_start = (total_size / size) * rank;
        const std::size_t my_count = (total_size / size);
        std::cout << "cField Reader of rank " << rank << " reading " << my_count
                  << " floats starting at element " << my_start << "\n";

        const adios2::Dims start{my_start};
        const adios2::Dims count{my_count};

        const adios2::Box<adios2::Dims> sel(start, count);
        field.resize(my_count);

        bp_cfield.SetSelection(sel);
        cfieldReader.Get(bp_cfield, field.data());
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
