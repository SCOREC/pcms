#ifndef COUPLING_H_
#define COUPLING_H_

#include <adios2.h>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include <Kokkos_Core.hpp>
#include <typeinfo>
#include <fftw3.h>

#include <couplingConstants.h>

namespace coupler {
  class Part1ParalPar3D
  void InitPart1ParalPar3D(Part1ParalPar3D  &p1pp3d);

  /** GO = global ordinate to count/number
   *  quantities over the entire domain
   */
  typedef long unsigned GO;

  /** LO = local ordinate to count/number
   *  quantities over a sub-domain
   */
  typedef unsigned LO;

  /** Storage of double precision 2D array data
   *  and associated meta data
   */
  class Array2d;  //defined in coupling.cc

  /** Destroy an Array2d object
   */
  void destroy(Array2d* a);

  /** Receive PreProc values from GENE
   */
  Array2d* receive_gene_pproc(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &engine);

  /** Receive density from GENE
   */
  Array2d* receive_density(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &engine);

  /** Sanity check values
   */
  void printSomeDensityVals(const Array2d* density);

  /** Send density to XGC
   */
  void send_density(const std::string cce_folder,
      const Array2d* density,
      adios2::IO &io, adios2::Engine &engine,
      adios2::Variable<double> &send_id);

  /** Receive field from XGC
   */
  Array2d* receive_field(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &eng);

  /** Send field to XGC
   */
  void send_field(const std::string cce_folder, const Array2d* field,
      adios2::IO &io, adios2::Engine &engine,
      adios2::Variable<double> &send_id);

  /** Close the Adios2 engine objects
   */
  void close_engines(adios2::Engine engine[]);
}//end namespace coupler

#endif
