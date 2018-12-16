#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <cstdarg>

#include <limits.h>
#include <windows.h>
#include <tchar.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "cuda_runtime.h"
#include "vector_types.h"

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

#ifdef __DRIVER_TYPES_H__
#ifndef DEVICE_RESET
#define DEVICE_RESET cudaDeviceReset();
#endif
#else
#ifndef DEVICE_RESET
#define DEVICE_RESET
#endif
#endif

#ifdef __DRIVER_TYPES_H__
static const char *_cudaGetErrorEnum(cudaError_t error) {
    return cudaGetErrorName(error);
}
#endif

template <typename T> void check(T result, char const *const func, const char *const file,
    int const line) {
    if (result) {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
            static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
        DEVICE_RESET
            // Make sure we call CUDA Device Reset before exiting
            exit(EXIT_FAILURE);
    }
}


__global__ void _kernelTest(double* output, int lyz, int lz)
{
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;

    //Imagine that we have many many work to do here, not just assign a 1 to it
    output[ix * lyz + iy * lz + iz] = 0.001;
}


extern "C" {
    void _cKernelCallConstantFunction(double* output)
    {
        dim3 dblock(2, 2, 2);
        dim3 dthread(2, 2, 2);
        _kernelTest << <dblock, dthread >> > (output, 16, 4);
    }
}

int main()
{
    double * pTable;
    cudaMalloc((void**)&pTable, sizeof(double) * 64);
    _cKernelCallConstantFunction(pTable);
    
    thrust::device_ptr<double> dp(pTable);
    thrust::device_vector<double> d_x(dp, dp + 64);

    double sum = thrust::reduce(d_x.begin(), d_x.end(), 0.0/*This is important! the return type is same as this*/, thrust::plus<double>());
    printf("result is %f\n", (float)sum);
    return 0;
}
