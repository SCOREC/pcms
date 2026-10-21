#!/usr/bin/env bash
# Builds pcms with gnet2gtc and installs it into $INSTALL_DIR.
# Run after: source env.sh && ./build-deps.sh
set -euo pipefail

: "${INSTALL_DIR:?source env.sh first}"

if [ ! -d "$SRC_DIR/pcms/.git" ]; then
  git clone "$PCMS_REPO" "$SRC_DIR/pcms"
fi
git -C "$SRC_DIR/pcms" fetch origin
git -C "$SRC_DIR/pcms" checkout "$PCMS_BRANCH"
git -C "$SRC_DIR/pcms" pull --ff-only origin "$PCMS_BRANCH"

# the pcms test suite needs its data repo at configure time
TEST_FLAGS=(-DBUILD_TESTING="$BUILD_TESTS")
if [ "$BUILD_TESTS" = "ON" ]; then
  if [ ! -d "$SRC_DIR/pcms_testcases/.git" ]; then
    git clone --depth 1 https://github.com/jacobmerson/pcms_testcases.git \
      "$SRC_DIR/pcms_testcases"
  fi
  TEST_FLAGS+=(-DPCMS_TEST_DATA_DIR="$SRC_DIR/pcms_testcases")
fi

echo "=== building pcms"
cmake -S "$SRC_DIR/pcms" -B "$BUILD_DIR/pcms" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
  -DCMAKE_C_COMPILER="$(which mpicc)" \
  -DBUILD_SHARED_LIBS=ON \
  -Dperfstubs_DIR="$INSTALL_DIR/lib/cmake" \
  -DPCMS_ENABLE_GNET=ON \
  -DPCMS_ENABLE_OMEGA_H=ON \
  -DPCMS_ENABLE_MESHFIELDS=OFF \
  -DPCMS_ENABLE_PETSC=OFF \
  -DPCMS_ENABLE_XGC=OFF \
  -DPCMS_ENABLE_C=OFF \
  -DPCMS_ENABLE_Fortran=OFF \
  -DPCMS_ENABLE_Python=OFF \
  -DPCMS_ENABLE_SPDLOG=OFF \
  "${TEST_FLAGS[@]}"
cmake --build "$BUILD_DIR/pcms" -j "$JOBS" --target install

echo "=== gnet2gtc installed: $INSTALL_DIR/bin/gnet2gtc"
