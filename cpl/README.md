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
cmake ../wdmapp_coupling/cpl
make # there will be a bunch of warnings
```

If all goes well you will have a `cpl` binary in the `src` directory.
