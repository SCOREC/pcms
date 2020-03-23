## Build on SCOREC RHEL7

### Setup

Clone the repo and create a build directory:

```
git clone git@github.com:SCOREC/wdmapp_coupling.git
mkdir build-wdmCoupler-rhel7
```

Run the following commands to select the compiler and set environment variables
for dependencies.  *This needs to be done every time you create a new shell and
build within it.*

```
module purge
module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
module use /opt/scorec/spack/dev/lmod/linux-rhel7-x86_64/Core
module load gcc mpich adios2
```

### Build

```
cd build-wdmCoupler-rhel7
cmake ../wdmapp_coupling/cpl -DCMAKE_CXX_COMPILER=mpicxx
make
```

If all goes well you will have a `cpl` binary in the `src` directory.

## Build on NERSC Cori

### Setup

Clone the repo and create a build directory:

```
git clone git@github.com:SCOREC/wdmapp_coupling.git
mkdir build-wdmCoupler-intel-cori
```

Run the following commands to select the compiler and set environment variables
for dependencies.  *This needs to be done every time you create a new shell and
build within it.*

```
a2=/project/projectdirs/m499/Software/adios2/DEFAULT/cori_haswell/DEFAULT
module load cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$a2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$a2/lib64
```

### Build

```
cd build-wdmCoupler-intel-cori
cmake ../wdmapp_coupling/cpl -DCMAKE_CXX_COMPILER=CC
make
```

If all goes well you will have a `cpl` binary in the `src` directory.

### Run

To run the `cpl` binary the `module load` and `export` commands need to be added
to your script that calls `srun`.  For example,

`runCoriCpl.sh`

```
#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=5
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell

a2=/project/projectdirs/m499/Software/adios2/DEFAULT/cori_haswell/DEFAULT
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$a2/lib64

srun /path/to/build-wdmCoupler-intel-cori/src/cpl
```

Where this script is submitted with `sbatch runCoriCpl.sh`.
