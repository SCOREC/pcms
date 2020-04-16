#include <coupling1.h>

namespace coupler {

/*  
  template<class T1,mpitype>
   MPIExchangeType(T1& bdesc, T2){
     MPI_Type_vector(bdesc.count,bdesc.lower,bdesc.n_points,mpitype, \
     &bdesc.mpi_lower_datatype);
     MPI_Type_commit(&bdesc.mpi_lower_datatype);
 
     MPI_Type_vector(bdesc.count,bdesc.lower,bdesc.n_points,mpitype, \
     &bdesc.mpi_upper_datatype);
     MPI_Type_commit(&bdesc.mpi_upper_datatype);

  } 
*/ 
// z is periodic
  template<class T>
void  zDensityBoundaryBufAssign(MPI_Comm mpitype,GO nzb,GO lx,GO ly,GO lz,T*** lowbuf,T*** upbuf,T*** box,Part1ParalPar3D &p1pp3d){  
  if(lowbuf==NULL||upbuf==NULL){
    std::cout<<"ERROR:the boundary buffer must be alloctted beforing involing this routine."
    std::exit();
  }
  if(p1pp3d.npz>1){
    if(lz>=nzb){
      mpisendrecv_aux2D(mpitype,lx,ly,lz,lowbuf,upbuf,box,p1pp3d);   
    } else if{
      std::cout<<"ERROR: nzb is larger than lzb. Large lzb is required."       
      }
  } else{
      if(p1pp3d.periods[2]==1){
        for(GO i=0;i<lx-1;i++){
          for(GO j=0;j<ly-1;j++){
            for(GO k=0;k<nzb;k++){
              lowbuf[i][j][k]=box[i][j][lz-nzb+k];
              upbuf[i][j][k]=box[i][j][k];
            }
          }
        }             
      } else{
          std::cout<<"The topology is not right."<<'\n';
          std::exit();
        }
    } 
}    


template<class T>
void mpisendrecv_aux2D(MPI_Comm mpitype,GO nzb,GO lx,GO ly,GO lz,T*** lowbuf,T*** upbuf,T*** box,Part1ParalPar3D &p1pp3d){
      T* sendbuf=new T[lx*ly*nzb];
      T* recvbuf=new T[lx*ly*nzb];
      MPI_STATUS status;
      int rank_souce,rank dest;
      MPI_Cart_shift(p1pp3d.comm_z,2,1,&rank_souce,&rank_dest);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb;k++)
              sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][k];
       }
      }        
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,100,recvbuf,lx*ly*nzb,mpitype,rank_source,101, \
                  p1pp3d.comm_z,&status);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb;k++)
              upbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
       }            
      }
      MPI_Cart_shift(p1pp3d.comm_z,2,-1,int rank_souce,int rank_dest);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb;k++)
              sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][lz-nzb+k];
       }
     }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,&rank_dest,102,recvbuf,lx*ly*nzb,mpitype,&rank_source,103, \
                  p1pp3d.comm_z,&status);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb;k++)
              lowbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k]
       }
      }      
    }
  delete[] sendbuf;
  delete[] recvbuf;
}
 

template<class T>
void mpisendrecv_aux1D(MPI_Comm mpitype,GO nzb,GO xind,GO yind,GO zind, T* lowbuf,T* upbuf,T* box1d,Part1ParalPar3D &p1pp3d){
      MPI_STATUS status;
      int rank_source, rank_dest;
      MPI_Cart_shift(p1pp3d.comm_z,3,1,&rank_source,&rank_dest);  
      MPI_Sendrecv(box1d[0],nzb,mpitype,rank_dest,100,upbuf[0],nzb,mpitype,rank_source,101, \
            p1pp3d.comm_z,&status);
      MPI_Cart_shift(p1pp3d.comm_z,3,-1,&rank_source,&rank_dest); 
      MPI_Sendrecv(box1d[zind-nzb],nzb,mpitype,rank_dest,102,lowbuf[0],nzb,mpitype,rank_source,103, \
            p1pp3d.comm_z,&status);
}  








}
