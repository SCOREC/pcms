# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-jammy:develop AS builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment && \
set -o noclobber \
&&  (echo spack: \
&&   echo '  # add package specs to the `specs` list' \
&&   echo '  specs:' \
&&   echo '  - fftw' \
&&   echo '  - pcms@develop+python~tests %gcc ^mpich ^omega-h+trilinos ^trilinos@16.1.0+exodus ^netcdf-c@4.8.1+mpi' \
&&   echo '  - py-numpy' \
&&   echo '  view: /opt/views/view' \
&&   echo '  concretizer:' \
&&   echo '    unify: true' \
&&   echo '  packages:' \
&&   echo '    llvm:' \
&&   echo '      externals:' \
&&   echo '      - spec: llvm@18.1.3+clang~flang~lld~lldb' \
&&   echo '        prefix: /usr' \
&&   echo '        extra_attributes:' \
&&   echo '          compilers:' \
&&   echo '            c: /usr/bin/clang-18' \
&&   echo '            cxx: /usr/bin/clang++-18' \
&&   echo '    gcc:' \
&&   echo '      externals:' \
&&   echo '      - spec: gcc@13.3.0 languages:='"'"'c,c++,fortran'"'"'' \
&&   echo '        prefix: /usr' \
&&   echo '        extra_attributes:' \
&&   echo '          compilers:' \
&&   echo '            c: /usr/bin/gcc' \
&&   echo '            cxx: /usr/bin/g++' \
&&   echo '            fortran: /usr/bin/gfortran' \
&&   echo '  config:' \
&&   echo '    install_tree:' \
&&   echo '      root: /opt/software' \
&&   echo '    build_jobs: 4') > /opt/spack-environment/spack.yaml

RUN spack repo add https://github.com/spack/spack-packages.git

# Clone and add the custom pcms spack repo
RUN git clone https://github.com/jacobmerson/pcms-spack.git /opt/pcms-spack && \
    spack repo add /opt/pcms-spack/spack_repo/pcms

# Patch omega-h package.py to add the SEACASExodus flag
RUN sed -i 's/args.append("-DOmega_h_USE_Trilinos:BOOL=ON")/args.append("-DOmega_h_USE_Trilinos:BOOL=ON")\n            args.append("-DOmega_h_USE_SEACASExodus:BOOL=ON")/' \
    $(spack location -p omega-h)/package.py

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh


# Bare OS image to run the installed executables
FROM ubuntu:22.04

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view

RUN echo "export PYTHONPATH=$(dirname $(ls $(spack location -i pcms)/lib/python*/site-packages/pcms*.so)):\$PYTHONPATH" >> /etc/bash.bashrc

ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

