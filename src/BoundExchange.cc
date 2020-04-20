#include <coupling1.h>

namespace coupler {

void InitBoundaryDescr3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)  
{
  bdesc.nzb=p1pp3d.nzb;
  bdesc.updenz=new double**[pipp3d.li0];
  bdesc.lowdenz=new double**[pipp3d.li0];
  for(int i=0;i<pipp3d.li0;i++){
    bdesc.updenz[i]=new double*[dp3d.part1lj0];
    bdesc.lowdenz[i]=new double*[dp3d.part1lj0];
    for(int j=0;j<dp3d.part1lj0;j++){
      bdesc.updenz[i][j]=new double[bdesc.nzb];
      bdesc.lowdenz[i][j]=new double[bdesc.nzb];
    }
  } 
  bdesc.uppotentz=new std::complex<double>**[p3m3d.xboxinds[0][pipp3d.mype_x]];
  bdesc.lowpotentz=new std::complex<double>**[p3m3d.xboxinds[0][pipp3d.mype_x]];
  bdesc.upzpart3=new double*[p3m3d.xboxinds[0][pipp3d.mype_x]];
  bdesc.lowzpart3=new double*[p3m3d.xboxinds[0][pipp3d.mype_x]];
  for(int i=0;i<p3m3d.xboxinds[0][pipp3d.mype_x]){
    bdesc.uppotentz[i]=new complex<double>*[p3m3d.lj0];
    bdesc.lowpotentz[i]=new complex<double>*[p3d3d.lj0]; 
    bdesc.upzpart3[i]=new double[bdesc.nzb];
    bdesc.lowzpart3[i]=new double[bdesc.nzb];
    for(int j=0;j<dp3d.part3lj0;j++){
      bdesc.uppotentz[i][j]=new complex<double>[bdesc.nzb];
      bdesc.lowpotentz[i][j]=new complex<double>[bdesc.nzb];
    }
  }
}

void zPotentBoundaryBufAssign(MPI_Datatype mpitype,BoundaryDescr3D &bdesc,DatasProc3D& dp3d, Part3Mesh3D& p3m3d,\
     Part1ParalPar3D &p1pp3d)
{
  if(lowbuf==NULL||upbuf==NULL){
    std::cout<<"ERROR:the boundary buffer of the potential must be allocated beforing invoking this routine."
    std::exit;
  }
  GO li0,lj0,lk0,nzb;
  li0=p3m3d.xboxinds[0][pipp3d.mype_x];
  lj0=p3m3d.lj0;
  nzb=bdesc.nzb;
  if(p1pp3d.npz>1){
    for(int i=0;i<li0;i++){
      lk0=p3m3d.mylk0[i];
      if(lk0<nzb){
         std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension."
         std::exit();
       } 
      if(p1pp3d.periods[2]==1){ 
        mpisendrecv_aux1D(mpitype,p1pp3d.comm_z,nzb,li0,lj0,lk0,bdesc.lowzpart3[i],bdesc.upzpart3[i], \
          p3m3d.pzcoords[i]); 
        if(p1pp3d.comm_z==0) 
        for(int j=0;j<lj0;j++){
          mpisendrecv_aux1D(mpitype,p1pp3d.comm_z,nzb,li0,lj0,lk0,bdesc.lowpotentz[i][j],bdesc.uppotentz[i][j],\
              dp3d.potentout[i][j]); 
        }
      } else {
         std::cout<<"The topology is not right for parallely domain."<<'\n';
         std::exit();
      }
    }  
  } else{
      if(p1pp3d.periods[2]==1){
        for(int i=0;i<li0;i++){
          lk0=p3m3d.mylk0[i]
          if(lk0<nzb){
            std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension."
            std::exit();
          }  
          for(GO k=0;k<nzb-1;k++){
            lowzpart3[i][k]=p3m3d.pzcoords[i][lk0-nzb+k];
            upzpart3[i][k]=p3m3d.pzcoords[i][k];
          }
          for(int j=0;j<lj0;j++){
            for(int k=0;k<nzb;k++){
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][k];
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][lk0-nzb+k];
            }  
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
