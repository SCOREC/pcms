#ifndef WDM_COUPLING_INTERPOLATION_H
#define WDM_COUPLING_INTERPOLATION_H
template <int order,typename T=double>
constexpr auto lagrange_interpolation(T* x, T* y, T xout ) noexcept{
  T yout = 0;
  for (int i=0; i<(order+1);++i) {
    double pj = 1;
    for(int j=0; j<(order+1); ++j) {
      if(j!=i){
        pj*= (xout-x[j])/(x[i]-x[j]);
      }
    }
    yout += y[i]*pj;
  }
  return yout;
}

#endif // WDM_COUPLING_INTERPOLATION_H
