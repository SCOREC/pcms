// Filename: cscan.cu
// nvcc -c -arch sm_13 cscan.cu

#include <thrust/device_vector.h>

#include <thrust/device_vector.h>

#include <thrust/scan.h>


extern "C" {

//scan for integer arrays

void scan_int_wrapper( int *data, int N)

{

// Wrap raw pointer with a device_ptr

thrust::device_ptr <int> dev_ptr(data);

// Use device_ptr in Thrust scan algorithm

thrust::inclusive_scan(dev_ptr, dev_ptr+N,dev_ptr);

}


//scan for float arrays

void scan_float_wrapper( float *data, int N)

{

thrust::device_ptr <float> dev_ptr(data);

thrust::inclusive_scan(dev_ptr, dev_ptr+N,dev_ptr);

}


//scan for double arrays

void scan_double_wrapper( double *data, int N)

{

thrust::device_ptr <double> dev_ptr(data);

thrust::inclusive_scan(dev_ptr, dev_ptr+N,dev_ptr);

}


}



