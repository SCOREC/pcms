#ifndef MULTIDIMARRAY_HPP
#define MULTIDIMARRAY_HPP

#include <Kokkos_Core.hpp>


using RealMatView = Kokkos::View<double**, Kokkos::LayoutRight>;  
using IntMatView = Kokkos::View<int**, Kokkos::LayoutRight>;  
using RealVecView = Kokkos::View<double*, Kokkos::LayoutRight>;
using IntVecView = Kokkos::View<int*, Kokkos::LayoutRight>;

using HostIntVecView = Kokkos::View<int*, Kokkos::HostSpace>;

KOKKOS_INLINE_FUNCTION
int calculateIndex(const IntVecView& dimensions, const int* indices) {
   
   int dim = dimensions.extent(0);
   int index = 0;
    int multiplier = 1;
    for (int i = dim - 1; i >= 0; --i) {
        index += indices[i] * multiplier;
        multiplier *= dimensions(i);
    }
    return index;
}


int calculateTotalSize(const HostIntVecView& dimensions){
   int dim = dimensions.extent(0);
   int size = 1;
   for (int i = 0; i < dim; ++i){
	size *= dimensions(i);
   }
   return size;
}
#endif 
