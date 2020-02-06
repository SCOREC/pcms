#/bin/bash -x

# edit these params
#CXX=mpicxx
CXX=CC

cmake $1 \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$PWD/install \
#      ${CONFIG_PARAMS}

