#include <coupling.h>

namespace coupler {

#define pi_ 3.1415926

  class Part1ParalPar3D

/*  template<class T> // T can be int, long int,or double
     Aarray1D{
      T* data   // data is 1d array

     };

  template<class T> // T can be int, long int, or double
    Array2D{
     T* data  // data is 2d array
    };
*/
  template<class T>  //T is the Aarray1D; num is an integer
    receive_field1D(T* array1D,const std::string cce_folder,
    const std::string name, GO num, MPI_Comm Comm){

  }

 template<class T> //T is the Arrary1D
   send_field1D(T* array1D, const std::string cce_folder, 
    const std::string name, GO num, MPI_Comm Comm){

 } 

  template<class T>  //T is the Aarray2D,num[2]
    receive_field2D(T* array2D,const std::string cce_folder,
    const std::string name,GO num, MPI_Comm Comm){


  }

 template<class T> //T is the Arrary2D, num[2]
   send_field2D(T* array2D, const std::string cce_folder, 
    const std::string name, GO num,MPI_Comm Comm){

 } 

 template<class T>
   receive_field1D_serial(T* array1D, const std::string cce_folder,const std::string name, GO num)

void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d)


}
