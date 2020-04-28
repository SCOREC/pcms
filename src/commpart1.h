#include "coupling.h"
namespace coupler {

class Part1ParalPar3D {
  public:
    LO mype ; // my process rank in mpi_comm_world
    LO mype_x; // my process rank in comm_x
    LO mype_y; // my process rank in comm_y
    LO mype_z; // my process rank in comm_z
    MPI_Comm comm_x;
    MPI_Comm comm_y;
    MPI_Comm comm_z;
    MPI_Comm comm_cart;
    LO  NP; // The total number of processes
    LO npx,npy,npz;
    LO nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
    LO ny0,nyb,lj0,lj1,lj2,lm0,lm1,lm2;
    LO nz0,nzb,lk0,lk1,lk2,ln0,ln1,ln2;
//    LO myli0,mylj1,myl12;  // The indexes of box y after Fourier transform
    int periods[3]={0,1,1};
    double* xcoords=NULL; // The 1d array storing the radial position of all flux surfaces
    double* pzcoords=NULL; // The 1d array storing the poloidal angle of all vertices along the poloidal surface curve.
    double* pzp=NULL; // The 1d array storing the poloial on each process.
    double dz;  // The equal step length along the poloidal flux curve.
 // parameters for creating the magnetic field, density and temperature ground. 
    LO res_fact;
    ~Part1ParalPar3D()
    {
      if(xcoords!=NULL)  delete[] xcoords;
      if(pzcoords!=NULL) delete[] pzcoords;
      if(pzp!=NULL)      delete[] pzp;
    }     
};

void InitPart1ParalPar3D(Part1ParalPar3D& p1pp3d);

template<class T>
void InputfromFile(T* numbers,LO ARRAY_SIZE,std::string filename)
{
    LO count = 0;             // Loop counter variable
    T x;
    LO rank;
    std::ifstream inputFile;        // Input file stream object
    // Open the file.
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
      std::cout<<filename<<'\n';
    }
    inputFile.open(filename);
    if (!inputFile) {
        std::cout << "Unable to open file";
        exit(1); // terminate with error
    }
    // Read the numbers from the file into the array.
    while (inputFile >> x){
       numbers[count]=x;
       count++;
    }
    inputFile.close();
    if(rank==0){
      for (count = 0; count < ARRAY_SIZE; count++){
	  std::cout << numbers[count] << "\n";
      }
    }
}

void InitPart1paral3DInCoupler(Part1ParalPar3D  &p1pp3d);

void CreateSubCommunicators(Part1ParalPar3D  &p1pp3d);

void MpiFreeComm(Part1ParalPar3D  &p1pp3d);
}











