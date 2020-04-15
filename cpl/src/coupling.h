#ifndef COUPLING_H_
#define COUPLING_H_

#include <adios2.h>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include <Kokkos_Core.hpp>
#include <typeinfo>
#include <fftw3.h>

#include "couplingConstants.h"

namespace coupler {
  class Part1ParalPar3D;
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
  template<class T> class Array2d {
    public:
      Array2d(GO gH, GO gW, GO lH, GO lW, GO start) :
        globH(gH), globW(gW), locH(lH), locW(lW), locFirstCol(start) {
          vals = new double[locH*locW];
      }   
      ~Array2d() {
        globH = globW = locH = locW = 0;
        delete [] vals;
      }   
      double val(long i) const {
        assert(i<(locH*locW));
        return vals[i];
      }   
      double* data() const { return vals; };
      GO globalH() const { return globH; };
      GO globalW() const { return globW; };
      GO localH() const { return locH; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      double* vals;
      GO globH;
      GO globW;
      GO locH;
      GO locW;
      GO locFirstCol;
  };

  template<class T>
  class Array1d {
    public:
      Array1d(GO gW, GO lW, GO start) :
        globW(gW), locW(lW), locFirstCol(start) {
        vals = new T[locW];
      }
      ~Array1d() {
        globW = locW = 0;
        delete [] vals;
      }
      T val(long i) const {
        assert(i<(locW));
        return vals[i];
      }
      T* data() const { return vals; };
      GO globalW() const { return globW; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      T* vals;
      GO globW;
      GO locW;
      GO locFirstCol;
  };
 
  /** Destroy an Array2d object
   */
  void destroy(Array1d<class T>* a);

  void destroy(Array2d<class T>* a);

  /** Receive PreProc values from GENE
   */
  Array1d<class T>* receive_gene_pproc(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &engine);

  /** Receive density from GENE
   */
  Array2d<class T>* receive_density(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &engine);

  /** Sanity check values
   */
  void printSomeDensityVals(const Array2d<class T>* density);

  /** Send density to XGC
   */
  void send_density(const std::string cce_folder,
      const Array2d<class T>* density,
      adios2::IO &io, adios2::Engine &engine,
      adios2::Variable<double> &send_id);

  /** Receive field from XGC
   */
  Array2d<class T>* receive_field(const std::string cce_folder,
      adios2::IO &io, adios2::Engine &eng);

  /** Send field to XGC
   */
  void send_field(const std::string cce_folder, const Array2d<class T>* field,
      adios2::IO &io, adios2::Engine &engine,
      adios2::Variable<double> &send_id);

  /** Close the Adios2 engine objects
   */
  void close_engines(adios2::Engine engine[], const int i);
}//end namespace coupler

#endif
