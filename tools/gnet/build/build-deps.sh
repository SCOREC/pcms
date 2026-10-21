#!/usr/bin/env bash
# Builds the pcms dependencies needed by gnet2gtc into $INSTALL_DIR.
# Run after: source env.sh
set -euo pipefail

: "${INSTALL_DIR:?source env.sh first}"

# Tested versions
KOKKOS_REF=4764a9a79                # develop, reports 4.5.99
KOKKOS_KERNELS_REF=4.6.01
OMEGA_H_REF=3cd28bcd                # scorec-v11.0.2 + 7
ADIOS2_REF=v2.10.2
REDEV_REF=1e4c876                   # v4.3.1 + 50
PERFSTUBS_REF=2676497
CATCH2_REF=v3.8.0

RELEASE_FLAGS=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG"
               -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" -DBUILD_SHARED_LIBS=ON)

# clone <url> <dir> <ref>: clone once, then check out the tested ref
clone() {
  local url=$1 dir=$SRC_DIR/$2 ref=$3
  if [ ! -d "$dir/.git" ]; then
    git clone "$url" "$dir"
  fi
  git -C "$dir" fetch --tags origin
  if ! git -C "$dir" rev-parse --verify --quiet "$ref^{commit}" >/dev/null; then
    echo "ERROR: $ref not found in $url" >&2
    exit 1
  fi
  git -C "$dir" checkout --detach "$ref"
}

# build <name> <cmake args...>: configure, build and install one package
build() {
  local name=$1
  shift
  echo "=== building $name"
  cmake -S "$SRC_DIR/$name" -B "$BUILD_DIR/$name" "$@"
  cmake --build "$BUILD_DIR/$name" -j "$JOBS" --target install
}

mkdir -p "$SRC_DIR" "$BUILD_DIR" "$INSTALL_DIR"

# Kokkos: CUDA through nvcc_wrapper, chosen automatically by Kokkos
clone https://github.com/kokkos/kokkos.git kokkos "$KOKKOS_REF"
build kokkos "${RELEASE_FLAGS[@]}" \
  -DCMAKE_CXX_STANDARD=17 \
  -DKokkos_ARCH_"$GPU_ARCH"=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_CUDA_LAMBDA=ON \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=ON

clone https://github.com/kokkos/kokkos-kernels.git kokkos-kernels "$KOKKOS_KERNELS_REF"
build kokkos-kernels "${RELEASE_FLAGS[@]}" \
  -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
  -DCMAKE_C_COMPILER="$(which mpicc)"

clone https://github.com/SCOREC/omega_h.git omega_h "$OMEGA_H_REF"
build omega_h "${RELEASE_FLAGS[@]}" \
  -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
  -DCMAKE_C_COMPILER="$(which mpicc)" \
  -DOmega_h_USE_MPI=ON \
  -DOmega_h_USE_Kokkos=ON \
  -DBUILD_TESTING=OFF

clone https://github.com/ornladios/ADIOS2.git ADIOS2 "$ADIOS2_REF"
build ADIOS2 "${RELEASE_FLAGS[@]}" \
  -DADIOS2_USE_MPI=ON \
  -DADIOS2_USE_Fortran=OFF \
  -DADIOS2_USE_Python=OFF \
  -DADIOS2_USE_HDF5=OFF \
  -DADIOS2_USE_CUDA=OFF \
  -DADIOS2_BUILD_EXAMPLES=OFF \
  -DBUILD_TESTING=OFF

clone https://github.com/UO-OACISS/perfstubs.git perfstubs "$PERFSTUBS_REF"
build perfstubs "${RELEASE_FLAGS[@]}"

clone https://github.com/SCOREC/redev.git redev "$REDEV_REF"
build redev "${RELEASE_FLAGS[@]}" \
  -DCMAKE_CXX_COMPILER="$(which mpicxx)" \
  -DCMAKE_C_COMPILER="$(which mpicc)" \
  -Dperfstubs_DIR="$INSTALL_DIR/lib/cmake" \
  -DBUILD_TESTING=OFF

if [ "$BUILD_TESTS" = "ON" ]; then
  clone https://github.com/catchorg/Catch2.git Catch2 "$CATCH2_REF"
  build Catch2 "${RELEASE_FLAGS[@]}" -DBUILD_TESTING=OFF
fi

echo "=== dependencies installed in $INSTALL_DIR"
