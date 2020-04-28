#include <coupling.h>

namespace coupler {
 
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









}

