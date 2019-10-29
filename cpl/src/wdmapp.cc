//#include <mpi.h>
#include<stdio.h>
#include"cpl_init.h"
#include<iostream>
#include<sstream> //for the stream run function call
//#include<cstdlib>
//#include<cstring>
#include<unistd.h>
#include<sys/wait.h>
#include<sys/stat.h>
#include<fstream>
#include<string>
//#include<string.h>

// TASK_1: create a new process to run the binary serially
// TASK_2: Using MPI, duplicate the number of process - currently run as mpirun -np executable called within benchmark
// TASK_3: Automate the reading of the number of processes needed
// TASK_4: Check the run scripts for the existing coupling and draw up something similar in cpp
// Task_5: Adopt the process for both GENE & XGC

struct run_coupling  
{
  private :
    const std::string XGC_lib/*, GENE_lib*/;
    struct stat buffer;

  public: 
    run_coupling(const std::string& X_lib /*, const std::string& G_lib*/)
      : XGC_lib(X_lib)/*, GENE_lib(G_lib) */ {}
    //gain access to the binaries via this functions
    const std::string& getNameX() const { return XGC_lib;}  
    // const std::string& getNameG() const { return GENE_lib;}  

    // check that the XGC library exists
    inline bool check_lib(const std::string& lib){
      return(stat (lib.c_str(), &buffer) == 0);
    }
    
    // MPI run is exactly as the standard run, calling MPIRUN
    //this was to create a single process to run the binary
    int ls(const char *dir){
      int pid, status;
      if (pid=fork()){
        waitpid(pid, &status, 0);
      } else {
        const char executable[] = "bin/ls";
        execl(executable, executable, dir, NULL);
      }
      return status;
    }


    int run_XGC(const char *dir){
      int pid, status;
      if (pid=fork()){
        waitpid(pid, &status, 0);
      } else {
        const char executable[] = "bin/ls";
        execl(executable, executable, dir, NULL);
      }
      return status;
    }


    int string_run()
    {
      std::stringstream stream;
      stream<< XGC_lib
        <<" "
        <<"argument";
 	return system(stream.str().c_str());
    }
};

int main(int argc, char **argv){
//	MPI_Init(NULL, NULL);
 //Ignore the use of OOP for this first iteration
  std::string MPI_RUN = "mpirun ";
  std::string PROC = "64 ";
  std::string XPATH = "../../WDM_XGC/xgc-es-coupling";
  std::string str = MPI_RUN + "-np " + PROC + XPATH;
  const char *command = str.c_str(); 
  system(command);
  system("../../benchmark/XGC/bin/coupling");
/*  run_coupling beta(lib1);
  beta.check_lib(lib1);
  std::cout << "ls'ing /" << std::endl;
  std::cout << "returned: " << beta.ls("/") << std::endl;*/
  return 0;

// MPI_Finalize();
}
