#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

#include "coupling.h"

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


}

#endif
