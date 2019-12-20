// Filename: csort.cu
// nvcc -c -arch sm_13 csort.cu

#include <thrust/device_vector.h>

#include <thrust/device_vector.h>

#include <thrust/sort.h>


extern "C" {

//Sort for integer arrays

void sort_int_wrapper( int *data, int N)

{

// Wrap raw pointer with a device_ptr

thrust::device_ptr <int> dev_ptr(data);

// Use device_ptr in Thrust sort algorithm

thrust::sort(dev_ptr, dev_ptr+N);

}


//Sort for float arrays

void sort_float_wrapper( float *data, int N)

{

thrust::device_ptr <float> dev_ptr(data);

thrust::sort(dev_ptr, dev_ptr+N);

}


//Sort for double arrays

void sort_double_wrapper( double *data, int N)

{

thrust::device_ptr <double> dev_ptr(data);

thrust::sort(dev_ptr, dev_ptr+N);

}


}



