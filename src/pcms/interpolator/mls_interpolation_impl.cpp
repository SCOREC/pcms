#include "mls_interpolation_impl.hpp"

namespace pcms
{

namespace detail
{
void calculate_basis_slice_lengths(IntHostMatView& array)
{
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int j = 0; j < dim; ++j) {
    array(0, j) = 1;
  }

  for (int i = 0; i < degree; ++i) {
    array(i, 0) = 1;
  }

  for (int i = 1; i < degree; ++i) {
    for (int j = 1; j < dim; ++j) {
      array(i, j) = array(i, j - 1) + array(i - 1, j);
    }
  }
}

int calculate_basis_vector_size(const IntHostMatView& array)
{
  int sum = 1;
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int i = 0; i < degree; ++i) {
    for (int j = 0; j < dim; ++j) {
      sum += array(i, j);
    }
  }

  return sum;
}

int calculate_scratch_shared_size()
{

  IntVecView shmem_each_team("stores the size required for each team",
                             nvertices_target);
  Kokkos::parallel_for(
    "calculate the size required for scratch for each team", nvertices_target,
    KOKKOS_LAMBDA(const int i) {
      int start_ptr = support.supports_ptr[i];
      int end_ptr = support.supports_ptr[i + 1];
      int nsupports = end_ptr - start_ptr;

      size_t total_shared_size = 0;

      total_shared_size += ScratchMatView::shmem_size(basis_size, basis_size);
      total_shared_size += ScratchMatView::shmem_size(basis_size, nsupports);
      total_shared_size += ScratchMatView::shmem_size(nsupports, basis_size);
      total_shared_size += ScratchVecView::shmem_size(basis_size);
      total_shared_size += ScratchVecView::shmem_size(nsupports) * 3;
      total_shared_size += ScratchMatView::shmem_size(nsupports, 2);
      total_shared_size += ScratchMatView::shmem_size(nsupports, 1);
      shmem_each_team(i) = total_shared_size;
    });

  // namespace KE = Kokkos::Experimental;
  // auto shared_size = KE::max_element(Kokkos::DefaultExecutionSpace(),
  // shmem_each_team); printf("shared size = %d", shared_size);
  int shared_size = 0;
  Kokkos::parallel_reduce(
    "find_max", nvertices_target,
    KOKKOS_LAMBDA(const int i, int& max_val_temp) {
      if (shmem_each_team(i) > max_val_temp) {
        max_val_temp = shmem_each_team(i);
      }
    },
    Kokkos::Max<int>(shared_size));

  return shared_size;
}

} // namespace detail
} // namespace pcms
