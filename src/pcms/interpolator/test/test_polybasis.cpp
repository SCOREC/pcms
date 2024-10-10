#include "pcms/interpolator/MLSInterpolation.hpp"

int dim = 3;
int degree = 2;
MatViewType d_slice_length("slice array",degree, dim);

auto host_slice_length = Kokkos::create_mirror_view(d_slice_length);

Kokkos::deep_copy(host_slice_length, 0)

Kokkos::deep_copy(d_slice_length, host_slice_length);


Kokkos::parallel_for(1, KOKKOS_LAMBDA (int k) {
   int m = k;
    basisSliceLengths(d_slice_length);
});

Kokkos::deep_copy(h_slice_length, d_slice_length);

for (int i = 0; i < host_slice_length.extent(0); ++i){
    for (int j = 0; j < host_slice_length.extent(1); ++j){
	std::cout << host_slice_length(i, j) << " ";
    }

	std::cout << "\n";
}
