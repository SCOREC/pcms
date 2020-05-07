## contents

src/
-README.md - this file
-BoundaryDescr3D.[h|cc] - BoundaryDescr3D class defines and supports
                          process boundary data exchange
-dataprocess.[h|cc] - DatasProc3D class supports interpolation and fourier
                      transform between meshes and fields
-dataprocess\_impl.h - templated functions supporting dataprocess.h
-fourierdataproc.cc - DatasProc3D fourier transform functions
-interpo.cc - ? interpolation functions
-sendrev\_imp.h - templated functions for mpi sendrecv
-CMakeLists.txt - cmake
-commpart1.cc - gene mesh class
-coupling.cc - adios2 send/recv functions and array2d class
-coupling.h - header for coupling.cc
-importpart3mesh.cc - xgc mesh class
-testutilities.h - enum for test cases and IO helper functions

test/
-cpl.cc - driver
-test\_init.cc - test constructors

