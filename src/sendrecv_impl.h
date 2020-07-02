#ifndef SENDRECV_IMPL_H
#define SENDRECV_IMPL_H

#include "couplingTypes.h"
#include <mpi.h>
#include <type_traits> // is_same
#include <cassert>

namespace {
//anonymous namespace is not accessible outside this file
//Gerrett - better way to hide this function?

//TODO I think there is a cleaner way - ask Gerrett
template<class T>
MPI_Datatype getMpiType(T foo) {
  MPI_Datatype mpitype;
  //determine the type based on what is being sent
  if( std::is_same<T, double>::value ) {
    mpitype = MPI_DOUBLE;
  } else if ( std::is_same<T, coupler::CV>::value ) { 
    //https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node48.htm
    mpitype = MPI_CXX_DOUBLE_COMPLEX;
  } else if ( std::is_same<T, coupler::GO>::value ) { 
    mpitype = MPI_UNSIGNED_LONG;
  } else if ( std::is_same<T, coupler::LO>::value ) {
    mpitype = MPI_INT;
  } else {
    assert(false);
    fprintf(stderr, "Unknown type in %s... exiting\n", __func__);
    exit(EXIT_FAILURE);
  }
  return mpitype;
}

}

namespace coupler {

template<class T>
void mpisendrecv_aux2D(MPI_Comm comm,LO nzb,LO lx,LO ly,LO lz,
     T*** lowbuf,T*** upbuf,T*** box)
{
      MPI_Datatype mpitype = getMpiType(T());
      T* sendbuf=new T[lx*ly*nzb];
      T* recvbuf=new T[lx*ly*nzb];
      MPI_Status status;
      int rank_source,rank_dest;
      MPI_Cart_shift(comm,0,1,&rank_source,&rank_dest);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++){
           sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][k];
         }
       }
      }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,100,recvbuf,lx*ly*nzb,mpitype,rank_source,100,
                  comm,&status);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++){
              upbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
         }
       }
      }
      MPI_Cart_shift(comm,0,-1, &rank_source, &rank_dest);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++)
            sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][lz-nzb+k];
       }
     }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,102,recvbuf,lx*ly*nzb,mpitype,rank_source,102,
                  comm,&status);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++)
              lowbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
       }
      }
  delete[] sendbuf;
  delete[] recvbuf;
}

template<class T>
void mpisendrecv_aux1D(MPI_Comm comm,LO nzb,LO xind,LO yind,LO zind,
     T* lowbuf,T* upbuf,T* box1d)
{

      MPI_Datatype mpitype = getMpiType(T());
      MPI_Status status;
      int rank_source, rank_dest;
      MPI_Cart_shift(comm,0,1,&rank_source,&rank_dest);
      MPI_Sendrecv(box1d,nzb,mpitype,rank_dest,100,upbuf,nzb,mpitype,rank_source,100,
            comm,&status);
      MPI_Cart_shift(comm,0,-1,&rank_source,&rank_dest);
      MPI_Sendrecv(&box1d[zind-nzb],nzb,mpitype,rank_dest,102,lowbuf,nzb,mpitype,rank_source,102,
            comm,&status);
}

}

#endif
