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

__global__ void _kInitialArray(int* thearray)
{
    int iX = threadIdx.x + blockDim.x * blockIdx.x;
    int iY = threadIdx.y + blockDim.y * blockIdx.y;
    int iZ = threadIdx.z + blockDim.z * blockIdx.z;

    thearray[iX * 16 + iY * 4 + iZ] = iX * 16 + iY * 4 + iZ;
}

__constant__ class B* _pB;

extern "C" {
    void _cInitialArray(int* thearray)
    {
        dim3 block(1, 1, 1);
        dim3 th(4, 4, 4);

        _kInitialArray << <block, th >> > (thearray);
        checkCudaErrors(cudaGetLastError());
    }
}

class B
{
public:
    B()
    {
        checkCudaErrors(cudaMalloc((void**)&m_pDevicePtr, sizeof(int) * 64));
        _cInitialArray(m_pDevicePtr);
    }
    ~B()
    {
        cudaFree(m_pDevicePtr);
    }
    __device__ int GetNumber(int index)
    {
        m_pDevicePtr[index] = m_pDevicePtr[index] + 1;
        return m_pDevicePtr[index];
    }
    __device__ __inline__ static int GetNumberS(int* ary, int index)
    {
        ary[index] = ary[index] + 1;
        return ary[index];
    }

    int* m_pDevicePtr;
};

__global__ void _kAddArray(int* thearray1, int* thearray2)
{
    int iX = threadIdx.x + blockDim.x * blockIdx.x;
    int iY = threadIdx.y + blockDim.y * blockIdx.y;
    int iZ = threadIdx.z + blockDim.z * blockIdx.z;

    thearray1[iX * 16 + iY * 4 + iZ] = thearray1[iX * 16 + iY * 4 + iZ] + B::GetNumberS(thearray2, iX * 16 + iY * 4 + iZ);
}

__global__ void _kAddArrayB(int* thearray1, B* pB)
{
    int iX = threadIdx.x + blockDim.x * blockIdx.x;
    int iY = threadIdx.y + blockDim.y * blockIdx.y;
    int iZ = threadIdx.z + blockDim.z * blockIdx.z;

    thearray1[iX * 16 + iY * 4 + iZ] = thearray1[iX * 16 + iY * 4 + iZ] + B::GetNumberS(pB->m_pDevicePtr, iX * 16 + iY * 4 + iZ);
}

extern "C" {
    void _cAddArray(int* thearray1, int* thearray2)
    {
        dim3 block(1, 1, 1);
        dim3 th(4, 4, 4);

        _kAddArray << <block, th >> > (thearray1, thearray2);
        checkCudaErrors(cudaGetLastError());
    }

    void _cAddArrayB(int* thearray1, B* pB)
    {
        dim3 block(1, 1, 1);
        dim3 th(4, 4, 4);

        _kAddArrayB << <block, th >> > (thearray1, pB);
        checkCudaErrors(cudaGetLastError());
    }

}

class A
{
public:
    A()
    {
        checkCudaErrors(cudaMalloc((void**)&m_pDevicePtr, sizeof(int) * 64));
        _cInitialArray(m_pDevicePtr);
    }
    ~A()
    {
        checkCudaErrors(cudaFree(m_pDevicePtr));
    }

    void Add(int* toAdd)
    {
        _cAddArray(m_pDevicePtr, toAdd);
    }

    void Add(B* toAdd)
    {
        _cAddArrayB(m_pDevicePtr, toAdd);
    }

    int* m_pDevicePtr;
};



int main(int argc, char * argv[])
{
    B* pB = new B();
    A* pA = new A();
    pA->Add(pB->m_pDevicePtr);

    int* res = (int*)malloc(sizeof(int) * 64);
    checkCudaErrors(cudaMemcpy(res, pA->m_pDevicePtr, sizeof(int) * 64, cudaMemcpyDeviceToHost));
    printf("----------- A=");
    for (int i = 0; i < 8; ++i)
    {
        printf("\n");
        for (int j = 0; j < 8; ++j)
            printf("res %d=%d  ", i * 8 + j, res[i * 8 + j]);
    }
    printf("\n");
    checkCudaErrors(cudaMemcpy(res, pB->m_pDevicePtr, sizeof(int) * 64, cudaMemcpyDeviceToHost));
    printf("----------- B=");
    for (int i = 0; i < 8; ++i)
    {
        printf("\n");
        for (int j = 0; j < 8; ++j)
            printf("res %d=%d  ", i * 8 + j, res[i * 8 + j]);
    }
    printf("\n");
    B* pB2 = new B();
    A* pA2 = new A();
    B* dpB2;
    cudaMalloc((void**)&dpB2, sizeof(B*));
    cudaMemcpy(dpB2, pB2, sizeof(B), cudaMemcpyHostToDevice);
    //cudaMemcpyToSymbol(_pB, pB2, sizeof(B*));
    pA2->Add(dpB2);
    checkCudaErrors(cudaMemcpy(res, pA2->m_pDevicePtr, sizeof(int) * 64, cudaMemcpyDeviceToHost));
    printf("----------- A2=");
    for (int i = 0; i < 8; ++i)
    {
        printf("\n");
        for (int j = 0; j < 8; ++j)
            printf("res %d=%d  ", i * 8 + j, res[i * 8 + j]);
    }
    printf("\n");
    checkCudaErrors(cudaMemcpy(res, pB2->m_pDevicePtr, sizeof(int) * 64, cudaMemcpyDeviceToHost));
    printf("----------- B2=");
    for (int i = 0; i < 8; ++i)
    {
        printf("\n");
        for (int j = 0; j < 8; ++j)
            printf("res %d=%d  ", i * 8 + j, res[i * 8 + j]);
    }
    printf("\n");
    delete pA;
    delete pB;
    delete pA2;
    delete pB2;

    return 0;
}