#ifndef COUPLING_TYPES_H
#define COUPLING_TYPES_H

#include <complex>
#include <vector>

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

  typedef std::vector<std::vector<std::vector<LO> > > vecint3d;
  typedef std::vector<std::vector<LO> > vecint2d;
  typedef std::vector<LO> vecint1d;
}

#endif
