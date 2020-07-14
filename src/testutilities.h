#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

#include "adios2Routines.h"

namespace coupler {

enum class TestCase : unsigned {
  off,
  t0, // better name?
  invalid
};

// input and output utilities
template<class T>
void OutputtoFile(T* numbers,LO array_size,std::string filename){ 
  std::ofstream outputFile; 
  LO rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  outputFile.open(filename);
  if(!outputFile){
    std::cout<<"unable to open: filename="<<filename<<'\n';
    exit(1);
  }
  for(int i=0;i<array_size;i++){
    outputFile<<numbers[i]<<'\n';
  }
  outputFile.close();  
}


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

// for diagnosis
template<class T>
void printSumm1D(T* array, GO inds1d[2],T sum,
     MPI_Comm &comm, std::string name,LO numiter)
{  LO num=0;
  LO rank;
  MPI_Comm_rank(comm,&rank);
  for(GO i=inds1d[0];i<inds1d[1]+1;i++){
      sum+=array[i];
      num+=1;
//     if(rank==0) std::cout<<i<<" "<<array[i]<<'\n';
    }
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<" "<<num<<'\n';
}



template<class T>
void printSumm2D(T** array, GO inds1d[2],GO inds2d[2],T sum,
     MPI_Comm comm, std::string name,LO numiter)
{ 
  for(GO i=inds1d[0];i<inds1d[1]+1;i++){
    for(GO j=inds2d[0];j<inds2d[1]+1;j++){
      sum+=array[i][j];
    }
  }
  LO rank;
  MPI_Comm_rank(comm,&rank);
  std::cout<<"numiter,rank "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}


template<class T>
void printSumm3D(T*** array, LO inds1d,LO inds2d,LO *inds3d,T sum,
     MPI_Comm comm, std::string name,LO numiter)
{
  for(GO i=0;i<inds1d;i++){
    for(GO j=0;j<inds2d;j++){
      for(GO k=0;k<inds3d[i];k++)      
        sum+=array[i][j][k];
    }
  }
  LO rank;
  MPI_Comm_rank(comm,&rank);
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}

template<class T>
void printSumm4D(T**** array, LO inds1d,LO inds2d,LO inds3d,LO* inds4d,T sum,
     MPI_Comm comm, std::string name,LO numiter)
{
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      for(LO k=0;k<inds3d;k++){
        for(LO l=0;k<inds4d[i];l++)
          sum+=array[i][j][k][l];
      }
    }
  }
  LO rank;
  MPI_Comm_rank(comm,&rank);
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}


}
#endif
