#ifndef TESTUTILITIES_H
#define TESTUTILITIES_H

#include "adios2Routines.h"

namespace coupler {

enum class TestCase : unsigned {
  off,
  t0, // better name?
  invalid
};

enum class CouplingCase : unsigned {
  genexgc,
  gemxgc,
  off
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
  for(GO i=inds1d[0];i<inds1d[1];i++){
      sum+=array[i];
      num+=1;
//     if(rank==0) std::cout<<i<<" "<<array[i]<<'\n';
    }
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<" "<<num<<'\n';
}



template<class T>
void printSumm2D(T** array, LO inds1d,LO inds2d,T sum,
     LO rank, std::string name,LO numiter)
{ 
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      sum+=array[i][j];
    }
  }
  printf("%s,numiter,rank,sum=%2d %2d %E \n",name.c_str(),numiter,rank,sum);  
//std::cout<<"numiter,rank "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}


template<class T>
void printSumm3D(T*** array, LO inds1d,LO inds2d,LO *inds3d,T sum,
     MPI_Comm comm, std::string name,LO numiter)
{
  LO rank;
  MPI_Comm_rank(comm,&rank); 
  for(GO i=0;i<inds1d;i++){
    for(GO j=0;j<inds2d;j++){
      for(GO k=0;k<inds3d[i];k++)      
        sum+=array[i][j][k];
    }
  }
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}

template<class T>
void printSumm4D(T**** array, LO inds1d,LO inds2d,LO inds3d,LO* inds4d,T sum,
     MPI_Comm comm, std::string name,LO numiter)
{
  LO rank;
  MPI_Comm_rank(comm,&rank); 
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      for(LO k=0;k<inds3d;k++){
        for(LO l=0;l<inds4d[i];l++){
          sum+=array[i][j][k][l];
//if(rank==0) std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<'\n';
        }   
      }
    }
  }
  std::cout<<"numiter,ranki "<<name<<"="<<numiter<<" "<<rank<<" "<<sum<<'\n';
}

template<class T>
void printminmax(T*** array, LO inds1d, LO inds2d, LO *inds3d,
     LO rank, std::string name, LO numiter)
{ double minmal=0.0;
  double maxmal=0.0;
  T minele,maxele;
  for(GO i=0;i<inds1d;i++){
    for(GO j=0;j<inds2d;j++){
      for(GO k=0;k<inds3d[i];k++){
        if(std::abs(array[i][j][k])>maxmal) maxele=array[i][j][k];
        if(std::abs(array[i][j][k])<minmal) minele=array[i][j][k];
      }
    }
  }
  std::cout<<name<<" numiter,rank,"<<" minele,maxele="<<numiter<<" "<<rank<<" "<<minele<<" "<<maxele<<'\n';
}

template<class T>
void printminmax_var3d(T*** array, LO inds1d, LO inds2d, LO *inds3d, LO rank, std::string name, 
     LO numiter, bool outpattern = true)
{ double minmal=0.0;
  double maxmal=0.0;
  T minele,maxele;
  LO i1 = 0, j1 = 0, k1 = 0;
  for(GO i=0;i<inds1d;i++){
    for(GO j=0;j<inds2d;j++){
      for(GO k=0;k<inds3d[i];k++){
        if(typeid(T) == typeid(CV)) {
          if(std::abs(array[i][j][k])>maxmal) maxele=array[i][j][k];
          if(std::abs(array[i][j][k])<minmal) minele=array[i][j][k];
        }
        else if(typeid(T) == typeid(double)) {
          if(maxmal < array[i][j][k]) maxmal = array[i][j][k];
          i1 = i;
          j1 = j;
          k1 = k;
          if(minmal > array[i][j][k]) minmal = array[i][j][k];
        }
      }
    }
  }
  if(typeid(T) == typeid(CV)) {
    std::cout<<name<<" numiter,rank,"<<" minele,maxele="<<numiter<<" "<<rank<<" "<<minele<<" "<<maxele<<'\n';
  }
  else if(typeid(T) == typeid(double)) {
    printf("name: %s, numiter: %d, rank: %d,i1: %d, j1: %d, k1: %d,  maxval: %f, minval: %f \n", 
            name, numiter, rank, i1, j1, k1, maxmal, minmal);
    if(!outpattern) {
      LO size;
      LO rank;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
      double *maxarray = new double[size];
      double *minarray = new double[size];

      MPI_Gather(&maxmal, 1, MPI_DOUBLE, maxarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&minmal, 1, MPI_DOUBLE, minarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
      double minmal = 0.0;
      double maxmal = 0.0;
      LO min = 0;
      LO max = 0;
      if(rank == 0) {
	for (LO i=0; i<size; i++ ) {
	  if(minmal > minarray[i]) {minmal = minarray[i]; max = i;}
	  if(maxmal < maxarray[i]) {maxmal = maxarray[i]; min = i;}
	}
	printf("name:%s, min: %d, max: %d, m: %d, tomaxmal: %f, tominmal: %f \n", name, min, max, 
        numiter, maxmal, minmal); 
      }
      free(maxarray);
      free(minarray);
   }
  }
}


template<class T>
void printminmax1d(T* array, LO inds1d, LO rank, std::string name, LO numiter, bool outpattern)
{
  double minmal = 0.0;
  double maxmal = 0.0;
  T minele, maxele;
  LO imax, imin;
  for (LO i=0; i< inds1d; i++) {
    if(typeid(T) == typeid(CV)){
      if (std::abs(array[i]) > maxmal) {
	maxmal = std::abs(array[i]);
	maxele = array[i]; 
	imax = i;
      }
      if(std::abs(array[i] < minmal)) {
	minmal=std::abs(array[i]);
	minele=array[i];        
      }
    }
    else if (typeid(T) == typeid(double)){
      if (array[i] > maxmal) {
	maxmal = array[i];
	imax = i;
      }
      if(array[i] < minmal) {
        minmal=array[i];
        imin = i;
      }
    }    
  }
  if(!outpattern) {
    printf("name:%s, numiter: %d, rank: %d, imax: %d, imin: %d, minmal: %19.13f, maxmal: %19.13f \n",
    name.c_str(),numiter,rank,imax, imin, minmal,maxmal);
//  } 
//  else {
    if (typeid(T) == typeid(double)){
      LO size;
      LO rank;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
      double *maxarray = new double[size];
      double *minarray = new double[size];

      MPI_Gather(&maxmal, 1, MPI_DOUBLE, maxarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&minmal, 1, MPI_DOUBLE, minarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
      double minmal = 0.0;
      double maxmal = 0.0;
      if(rank == 0) {
	for (LO i=0; i<size; i++ ) {
	  if(minmal > minarray[i]) minmal = minarray[i];
	  if(maxmal < maxarray[i]) maxmal = maxarray[i];
	}
	printf("m: %d, name: %s, tomaxmal: %19.13f, tominmal: %19.13f \n", numiter, name, maxmal, minmal); 
      }
      free(maxarray);
      free(minarray);
    }
  }
}

template<class T>
void printminmax2d(T** array, LO inds1d,LO inds2d,LO rank, std::string name,LO numiter)
{ double minmal=0.0;
  double maxmal=0.0;
  T minele,maxele;
  LO i1,j1;
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      if(std::abs(array[i][j])>maxmal){
	maxmal=std::abs(array[i][j]);
	maxele=array[i][j];
	i1=i;
	j1=j;
       }
      if(std::abs(array[i][j])<minmal){
	minmal=std::abs(array[i][j]);
	minele=array[i][j];
      }
    }
  }
//  std::cout<<name<<" numiter,rank,"<<" minele,maxele="<<numiter<<" "<<rank<<" "<<minele<<" "<<maxele<<'\n';
  printf("%s,numiter,rank,minele,maxele=%2d %2d %E %E \n",name.c_str(),numiter,rank,minele,maxele);
  std::cout<<name<<" "<<rank<<" "<<i1<<" "<<j1<<'\n';
}

template<class T>
void printminmax3d(T*** array, LO inds1d,LO inds2d,LO inds3d,LO rank, std::string name, 
     LO numiter, bool outpattern = true)
{ double minmal=0.0;
  double maxmal=0.0;
  T minele,maxele;
  LO i1,j1,k1;
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      for(LO k=0;k<inds3d;k++){
        if(typeid(T) == typeid(CV)) {
	  if(std::abs(array[i][j][k])>maxmal){
	    maxmal=std::abs(array[i][j][k]);
	    maxele=array[i][j][k];
	    i1=i;
	    j1=j;
	    k1=k;
	   }
	  if(std::abs(array[i][j][k])<minmal){
	    minmal=std::abs(array[i][j][k]);
	    minele=array[i][j][k];
	  }
        }
        else if(typeid(T) == typeid(double)) {
          if(array[i][j][k] > maxmal) maxmal = array[i][j][k];
          if(array[i][j][k] < minmal) minmal = array[i][j][k];
        } 
     }
   }
 }

// printf("numiter: %d, rank: %d, 3Dminval: %f, 3dmaxval: %f \n", numiter, rank, minmal, maxmal);
//  std::cout<<name<<" numiter,rank,"<<" 3Dminval,3Dmaxval="<<numiter<<" "<<rank<<" "<<minmal<<" "<<maxmal<<'\n';
//  std::cout<<name<<" "<<rank<<" "<<i1<<" "<<j1<<" "<<k1<<'\n';
  if(!outpattern) {
    if (typeid(T) == typeid(double)){
      LO size;
      LO rank;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      double *maxarray = new double[size];
      double *minarray = new double[size];

      MPI_Gather(&maxmal, 1, MPI_DOUBLE, maxarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&minmal, 1, MPI_DOUBLE, minarray, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      double minmal = 0.0;
      double maxmal = 0.0;
      if(rank == 0) {
        for (LO i=0; i<size; i++ ) {
          if(minmal > minarray[i]) minmal = minarray[i];
          if(maxmal < maxarray[i]) maxmal = maxarray[i];
        }
        printf("rank: %d, name: %s, m: %d, to3Dmaxmal: %19.18f, to3Dminmal: %19.18f \n", rank, name, numiter, maxmal, minmal);
      }
      free(maxarray);
      free(minarray);
    }
  }
}


template<class T>
void printminmax4d(T**** array, LO inds1d,LO inds2d,LO inds3d,LO inds4d,
     MPI_Comm comm, std::string name,LO numiter)
{ double minmal=0.0;
  double maxmal=0.0;
  T minele,maxele;
  LO i1,j1,k1,l1;
  for(LO i=0;i<inds1d;i++){
    for(LO j=0;j<inds2d;j++){
      for(LO k=0;k<inds3d;k++){
        for(LO l=0;l<inds4d;l++){
	  if(std::abs(array[i][j][k][l])>maxmal){
	    maxmal=std::abs(array[i][j][k][l]);
	    maxele=array[i][j][k][l];  
            i1=i; 
            j1=j;
            k1=k;
            l1=l;
	   }
	  if(std::abs(array[i][j][k][l])<minmal){ 
	    minmal=std::abs(array[i][j][k][l]);
	    minele=array[i][j][k][l];
	  }
	}
      }
    }
  }
  LO rank;
  MPI_Comm_rank(comm,&rank);
  std::cout<<name<<" numiter,rank,"<<" minele,maxele="<<numiter<<" "<<rank<<" "<<minele<<" "<<maxele<<'\n';
  std::cout<<name<<" "<<rank<<" "<<i1<<" "<<j1<<" "<<k1<<" "<<l1<<'\n';
}



/*for reorder the datas between different collectives*/
template<class T>
void reshuffleforward(T* array,const LO nstart,const LO vertnum)
{
  T* tmp=new T[vertnum];
  for(LO i=0;i<vertnum-nstart;i++)
    tmp[i]=array[nstart+i];
  for(LO j=vertnum-nstart;j<vertnum;j++)
    tmp[j]=array[j-vertnum+nstart];
  for(LO k=0;k<vertnum;k++)
    array[k]=tmp[k];
  delete[] tmp;
}


template<class T>
void reshufflebackward(T* array,const LO nstart,const LO vertnum)
{
  T* tmp=new T[vertnum];
  for(LO i=vertnum-nstart;i<vertnum;i++)
    tmp[i-vertnum+nstart]=array[i];
  for(LO j=0;j<vertnum-nstart;j++)
    tmp[j+nstart]=array[j];
  for(LO k=0;k<vertnum;k++)
    array[k]=tmp[k];
  delete[] tmp;
}

std::string filename(std::string const& s);


}
#endif  /*testutilities.h*/
