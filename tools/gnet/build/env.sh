# Environment for building gnet2gtc. Edit this file, then:
#   source env.sh
#   ./build-deps.sh
#   ./build-pcms.sh

# -----------------------------------------------------------------------------
# 1. Compilers and libraries
#    Load these before sourcing this file (module names differ per system):
#      GCC with C++20 support   (tested 12.3)
#      MPI, mpicc/mpicxx        (tested MPICH 4.1.1)
#      CUDA, nvcc               (tested 12.1.1)
#      CMake >= 3.20            (tested 3.26)
#      HDF5 C library           (tested 1.14.1)
#      git
#    For example:
#
# module load gcc/12.3.0 mpich/4.1.1 cuda/12.1.1 cmake/3.26.3 hdf5/1.14.1
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2. GPU architecture (Kokkos name), see ../README.md for the full table:
#    VOLTA70 (V100), AMPERE80 (A100), AMPERE86 (RTX 30xx, A40),
#    ADA89 (RTX 40xx, L4, L40), HOPPER90 (H100)
# -----------------------------------------------------------------------------
export GPU_ARCH="${GPU_ARCH:-AMPERE80}"

# -----------------------------------------------------------------------------
# 3. Where to put everything
# -----------------------------------------------------------------------------
export GNET_ROOT="${GNET_ROOT:-$HOME/gnet2gtc}"
export SRC_DIR="$GNET_ROOT/src"          # git checkouts
export BUILD_DIR="$GNET_ROOT/build"      # build trees
export INSTALL_DIR="$GNET_ROOT/install"  # everything is installed here
export JOBS="${JOBS:-8}"                 # parallel build jobs

# -----------------------------------------------------------------------------
# 4. pcms source with gnet2gtc
# -----------------------------------------------------------------------------
export PCMS_REPO="${PCMS_REPO:-https://github.com/SCOREC/pcms.git}"
export PCMS_BRANCH="${PCMS_BRANCH:-gnet-gtc-coupling}"

# set to ON to also build the pcms unit tests (needs Catch2, longer build)
export BUILD_TESTS="${BUILD_TESTS:-OFF}"

# -----------------------------------------------------------------------------
# Derived settings (no need to edit)
# -----------------------------------------------------------------------------
for tool in mpicc mpicxx nvcc cmake git; do
  command -v "$tool" >/dev/null || echo "WARNING: $tool not found, load it first (see section 1)"
done
export NVCC_WRAPPER_DEFAULT_COMPILER="$(which mpicxx 2>/dev/null)"
export CMAKE_PREFIX_PATH="$INSTALL_DIR${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export PATH="$INSTALL_DIR/bin:$PATH"
