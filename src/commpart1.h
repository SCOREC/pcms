
namespace coupler{
class Part1ParalPar3D {
  public:
    GO mype ; // my process rank in mpi_comm_world
    GO mype_x; // my process rank in comm_x
    GO mype_y; // my process rank in comm_y
    GO mype_z; // my process rank in comm_z
    GO comm_x;
    GO comm_y;
    GO comm_z;
    GO NP; // The total number of processes
    GO npx,nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
    GO npy,ny0,nyb,lj0,lj1,lh2,lm0,lm1,lm2;
    GO npz,nz0 nzb,lk0,lk1,lk2,ln0,ln1,ln2;
    GO mylj0,mylj1,myl12;  // The indexes of box y after Fourier transform
    GO periods[]={0,1,1};
    double* xcoords; // The 1d array storing the radial position of all flux surfaces
    double* pzcoords; // The 1d array storing the poloidal angle of all vertices along the poloidal surface curve.
    double* pzp; // The 1d array storing the poloial on each process.
    double dz;  // The equal step length along the poloidal flux curve.
}

void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d);

}
