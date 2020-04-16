#include <coupling1.h>

namespace coupler {
 
class BoundaryDescr3D{
  GO nzb;
  double*** updenz==NULL; // The upper  boundary buffer on z domain for interpolation and storing the real quantiies resulted from the backward Fourier transform of complex charged density.
  double*** lowdenz==NULL; 
  std::complex<double>*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  quantiies resulted from the forward Fourier transform of electrosttic potential.
  std::complex<double>*** lowpotentz==NULL;
}

void InitBoundaryDescr3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)  
  bdesc.nzb=p1pp3d.nzb;
  updenz=new double**[pipp3d.li0];
  lowdenz=new double**[pipp3d.li0];
  for(int i=0;i<pipp3d.li0;i++){
    updenz[i]=new double*[dp3d.part1lj0];
    lowdenz[i]=new double*[dp3d.part1lj0];
    for(int j=0;j<dp3d.part1lj0;j++)
      updenz[i][j]=new double[bdesc.nzb];
      lowdenz[i][j]=new double[bdesc.nzb];
  } 
  uppotentz=new std::complex<double>**[p3m3d.xboxinds[0][pipp3d.mype_x]];
  lowpotentz=new std::complex<double>**[p3m3d.xboxinds[0][pipp3d.mype_x]];
  for(int i=0;i<p3m3d.xboxinds[0][pipp3d.mype_x]){
    uppotentz[i]=new std::complex<double>*[dp3d.part3lj0];
    lowpotentz[i]=new std::complex<double>*[dp3d.part3lj0]; 
    for(int j=0;j<dp3d.part3lj0;j++)
      uppotentz[i][j]=new std::complex<double>[bdesc.nzb];
      lowpotentz[i][j]=new std::complex<double>[bdesc.nzb];
  }


void zPotentBoundaryBufAssign(MPI_Comm mpitype,BoundaryDescr3D &bdesc,DatasProc3D$ dp3d, Part3Mesh3D& p3m3d,\
     Part1ParalPar3D &p1pp3d){
  if(lowbuf==NULL||upbuf==NULL){
    std::cout<<"ERROR:the boundary buffer of the potential must be allocated beforing invoking this routine."
    std::exit();
  }
  GO li0,lj0,lk0,nzb;
  li0=p3d3d.li0,lj0=p3d3d.j0,nzb=bdesc.nzb;
  if(p1pp3d.npz>1){
    for(int i=0;i<li0;i++){
      lk0=p3m3d.mylk0(i)
      if(lk0<nzb){
         std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension."
         std::exit();
       } 
      for(int j=0;j<lj0;j++){
        if(p1pp3d.periods[2]==1){
          mpisendrecv_aux1D(mpitype,nzb,li0,lj0,lk0,bdesc.lowpotentz[i][j][0],bdesc.uppotentz[i][j][0],\
          dp3d.potentout[i][j][0],p1pp3d);  
         } else {
           std::cout<<"The topology is not right for parallely domain."<<'\n';
           std::exit();
         }
       }  
    }
  } else{
      if(p1pp3d.periods[2]==1){
        for(int i=0;i<li0;i++){
          lk0=p3m3d.mylk0(i)
          if(lk0<nzb){
            std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension."
            std::exit();
          }  
          for(int j=0;j<lj0;j++){
            for(int k=0;k<nzb;k++)
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][k];
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][lk0-nzb+k];
           }     
       }
     } else {
         std::cout<<"The topology is not right for serial y domain."<<'\n';
         std::exit();
     }
   }
 }
}
/*  stuct boundarydescrip{
    public:
      GO lower; //number of boundary points on lower side
      GO upper; //number of boundary points on upper side
      GO count; //how often exchange the direction
      GO n_points; //number of points(incl boundary) in 1st dimension
      GO innerfirst; //index of last inner point, when counting starts with 0 for the most left point
      GO innerlast; //
      GO subarray_size; //size of the subarrays which consists the dimension to exchange
      GO n_upper_subarrays; //number of subarrays in the upper boundary
      GO n_lower_subarrays; //number of subarrays in the lower boundary
      bool mpi_type_allocated; //is true if two mpi datatypes have been set and committed
      GO mpi_lower_datatype; //user defined datatypes for the exchange of lower boundary
      GO mpi_upper_datatype; //user defined datatypes for the exchange of upper boundary
      GO exchange_direction; //
  }

  void InitBoundDescrip(boundarydescrip& bdesc, GO a_suarray_size, \
   GO a_n_subarrays, GO a_n_lower_subarrays,GO a_n_upper_subarrays, \
   GO a_count){
     bdesc.lower = a_subarray_size*a_n_lower_subarrays;
     bdesc.upper = a_subarray_size*a_n_upper_subarrays;
     bdesc.count = a_count;
     bdesc.n_points=a_n_subarrays*a_n_lower_subarrays;
     bdesc.innerfirst=bdesc.n_points-bdesc.upper-1;
     bdesc.subarray_size=a_subarray_size;
     bdesc.n_upper_subarrays=a_n_upper_subarrays;
     bdesc.n_lower_subarrays=a_n_lower_subarrays;
     bdesc.mpi_type_allocated = false;  
  }

  void  MPIExchangeType(boundarydescrip& bdesc, MPI_Comm mpitype){
     MPI_Type_vector(bdesc.count,bdesc.lower,bdesc.n_points,mpitype, \
     &bdesc.mpi_lower_datatype);
     MPI_Type_commit(&bdesc.mpi_lower_datatype);

     MPI_Type_vector(bdesc.count,bdesc.lower,bdesc.n_points,mpitype, \
     &bdesc.mpi_upper_datatype);
     MPI_Type_commit(&bdesc.mpi_upper_datatype);

  } 
*/ 
