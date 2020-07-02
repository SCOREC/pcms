#ifndef COUPLING_TYPES_H
#define COUPLING_TYPES_H

#include <complex>

namespace coupler { 
  /** GO = global ordinate to count/number
   *  quantities over the entire domain
   *  TODO should be long int to match LO
   */
  typedef unsigned long GO;
  /** LO = local ordinate to count/number
   *  quantities over a sub-domain
   */
  typedef int LO;
  /** CV = complex value 
   */
  typedef std::complex<double> CV;
}

#endif
