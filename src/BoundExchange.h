#include <coupling1.h>

namespace coupler {

class BoundaryDescr3D{
  GO nzb;
  double*** upzpart3;
  double*** lowzpart3;
  double*** updenz==NULL; // The upper  boundary buffer on z domain for interpolation and storing the real quantiies resulted from the backward Fourier transform of complex charged density.
  double*** lowdenz==NULL;
  std::complex<double>*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  quantiies resulted from the forward Fourier transform of electrosttic potential.
  std::complex<double>*** lowpotentz==NULL;
}

void InitBoundaryDescr3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d);

void zPotentBoundaryBufAssign(MPI_Datatype mpitype,BoundaryDescr3D &bdesc,DatasProc3D& dp3d, Part3Mesh3D& p3m3d,\
     Part1ParalPar3D &p1pp3d)



  template<class T>
void  zDensityBoundaryBufAssign(MPI_Datatype mpitype,GO nzb,GO lx,GO ly,GO lz,T*** lowbuf, \
       T*** upbuf,T*** box,Part1ParalPar3D &p1pp3d)
  {  
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
            for(GO k=0;k<nzb-1;k++){
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
void mpisendrecv_aux2D(MPI_Datatype mpitype,MPI_Comm comm,GO nzb,GO lx,GO ly,GO lz,T*** lowbuf,T*** upbuf,T*** box)
{
      T* sendbuf=new T[lx*ly*nzb];
      T* recvbuf=new T[lx*ly*nzb];
      MPI_STATUS status;
      int rank_souce,rank dest;
      MPI_Cart_shift(comm,2,1,&rank_souce,&rank_dest);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb-1;k++)
           sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][k];
       }
      }        
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,100,recvbuf,lx*ly*nzb,mpitype,rank_source,101, \
                  comm,&status);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb-1;k++)
              upbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
       }            
      }
      MPI_Cart_shift(comm,2,-1,int rank_souce,int rank_dest);
      for(GO i=0;i<lx-1;i++){
       for(GO j=0;j<ly-1;j++){
         for(GO k=0;k<nzb-1;k++)
            sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][lz-nzb+k];
       }
     }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,&rank_dest,102,recvbuf,lx*ly*nzb,mpitype,&rank_source,103, \
                  comm,&status);
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
void mpisendrecv_aux1D(MPI_Datatype mpitype,MPI_Comm comm,GO nzb,GO xind,GO yind,GO zind, T* lowbuf,T* upbuf,T* box1d)
{
      MPI_STATUS status;
      int rank_source, rank_dest;
      MPI_Cart_shift(comm,3,1,&rank_source,&rank_dest);  
      MPI_Sendrecv(box1d[0],nzb,mpitype,rank_dest,100,upbuf[0],nzb,mpitype,rank_source,101, \
            comm,&status);
      MPI_Cart_shift(comm,3,-1,&rank_source,&rank_dest); 
      MPI_Sendrecv(box1d[zind-nzb],nzb,mpitype,rank_dest,102,lowbuf[0],nzb,mpitype,rank_source,103, \
            comm,&status);
}  

}
