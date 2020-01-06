//#include <mpi.h>
#include"cpl_init.h"
#include<sstream> //for the stream run function call
#include<unistd.h>
#include<sys/wait.h>
#include<sys/stat.h>
#include<string>

// TASK_1: Execute both codes using the system commands 
// TASK_2: validate OOP implementation is necessary i.e. considering process level coordination
// TASK_3: Compare both codes

struct couple  
{
  private :
    const std::string XGC_lib;
    const int N_X, n_X;
    struct stat buf;

  public: 
    couple(const std::string& X_lib, const int& nodes_X, const int& threads_X )
      : XGC_lib(X_lib), N_X(nodes_X), n_X(threads_X) {}

    //access binaries 
    const std::string& getNameX() const { return XGC_lib;}  

    // check that the XGC library exists
    inline bool check_libs(){ // if checking const char* file-> stat(file,&buf)
      return(stat (XGC_lib.c_str(), &buf) == 0);// use "continue" when checking for GENE
    }
    
    // standard srun
    int run_X()
    {
      std::stringstream stream;
      stream << "srun -N " << N_X  <<" -n " << n_X <<" "  << XGC_lib <<" > XGC.output 2> XGC.error";
 	return system(stream.str().c_str());
    }
    
    //this creates a single process to run the binary
    int run_X_2(const char *dir){
      int pid, status;
      if (pid=fork()){
        waitpid(pid, &status, 0);
      } else {
        const char* executable = XGC_lib.c_str();
        std::string str = "srun -N " + N_X;//+ "-n " + n_X + " " +  XGC_lib + " > XGC.output 2> XGC.error";
        execl(executable, str.c_str(), (char*)0);
      }
      return status;
    }

};

int main(int argc, char **argv){
//	MPI_Init(NULL, NULL);

  std::string SRUN = "srun ", PROC_X = "16", PROC_G = "16", THREAD_X = "512", THREAD_G="384";
  // full proc_x => 256 & thread_x => 6144
  
  std::string X_PATH = "~/XGC_adios2/xgc-es", X_OPT = " > XGC.output 2> XGC.error";
  std::string G_PATH = "~/gene_Adios2_cori/bin/gene_cori", G_OPT = " > GENE.output 2> GENE.error";

  std::string str = SRUN + "-N " + PROC_X + "-n " + THREAD_X + " "+ X_PATH + X_OPT;
  const char *command = str.c_str(); 
  system(command);

  //system("wait");

  str = SRUN + "-N " + PROC_G + "-n " + THREAD_G + G_PATH + G_OPT;
  command = str.c_str(); 
  system(command);
//	                  -----> XGC 
//                       | 
// -----> BENCHMARK ---->
// (coupling binary loc) |
// 	                  -----> GENE

  return 0;

// MPI_Finalize();
}
