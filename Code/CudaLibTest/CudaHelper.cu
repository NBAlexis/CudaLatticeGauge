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

__global__ void _kernelInitial(int* output, int lyz, int lz)
{
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;

    output[ix * lyz + iy * lz + iz] = 1;
}


__global__ void _kernelPluseOne(int* output, int lyz, int lz)
{
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;

    output[ix * lyz + iy * lz + iz] += 1;
}

__global__ void _kernelMultiplyTwo(int* output, int lyz, int lz)
{
    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int iz = threadIdx.z + blockDim.z * blockIdx.z;

    output[ix * lyz + iy * lz + iz] *= 2;
}


extern "C" {
    void _cKernelInitial(int* output)
    {
        dim3 dblock(2, 2, 1);
        dim3 dthread(8, 8, 1);
        _kernelInitial << <dblock, dthread >> > (output, 16, 1);
    }

    void _cKernelPlusOne(int* output)
    {
        dim3 dblock(2, 2, 1);
        dim3 dthread(8, 8, 1);
        _kernelPluseOne << <dblock, dthread >> > (output, 16, 1);
    }

    void _cKernelMultiplyTwo(int* output)
    {
        dim3 dblock(2, 2, 1);
        dim3 dthread(8, 8, 1);
        _kernelMultiplyTwo << <dblock, dthread >> > (output, 16, 1);
    }
}

int main()
{
    int * pTable;
    cudaMalloc((void**)&pTable, sizeof(int) * 256);
    _cKernelInitial(pTable);
    _cKernelPlusOne(pTable);
    _cKernelMultiplyTwo(pTable);
    _cKernelPlusOne(pTable);
    _cKernelMultiplyTwo(pTable);

    int outData[256];
    cudaMemcpy(outData, pTable, sizeof(int) * 256, cudaMemcpyDeviceToHost);
    
    printf("res=\n");
    for (int i = 0; i < 16; ++i)
    {
        for (int j = 0; j < 16; ++j)
        {
            printf("%d ", outData[i * 16 + j]);
        }
        printf("\n");
    }

    cudaFree(pTable);
    return 0;
}
