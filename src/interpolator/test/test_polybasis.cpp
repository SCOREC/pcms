#include "MLSInterpolation.hpp"

int dim = 3;
int degree = 2;
MatViewType d_slice_length("slice array",degree, dim);

auto host_slice_length = Kokkos::create_mirror_view(d_slice_length);

Kokkos::deep_copy(d_slice_length, host_slice_length);

basisSLiceLengths(host_slice_length);

for (i = 0; i < host_slice_length.extent(0); ++i){
    for (j = 0; j < host_slice_length.extent(1); ++j){
	std::cout << host_slice_length(i, j) << " ";
    }

	std::cout << "\n";
}


