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

#include "cuComplex.h"

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

#define _TESTLENGTH_ (65536)
#define _ONE_LINE_START_ (12345)
#define _ONE_LINE_END_ (12360)

struct alignas(8) ClassA
{
public:
    __device__ ClassA(): m_iValue(0), m_bySomeTag(0) { ; }
    __device__ ClassA(int v) : m_iValue(v), m_bySomeTag(static_cast<byte>(v & 0x11)) { ; }
    __device__ void Add(int v)
    {
        m_iValue += v;
        m_bySomeTag = static_cast<byte>(m_iValue & 0x11);
    }
    __device__ void DebugPrint() const
    {
        printf("v=%d,%d; ", m_iValue, static_cast<int>(m_bySomeTag));
    }
    int m_iValue;
    byte m_bySomeTag;

    byte _nouse[3];
};

class ClassB
{
public:
    __device__ ClassB() { ; }
    __device__ virtual void InitialArray(ClassA * pArray) = 0;
};

class ClassC : public ClassB
{
public:
    __device__ ClassC() : ClassB() { ; }
    __device__ virtual void InitialArray(ClassA * pArray)
    {
        for (int i = 0; i < _TESTLENGTH_; ++i)
        {
            pArray[i] = ClassA(i);
        }
    }
};

__constant__ ClassB* __pWorker;

__global__ void InitialConstant(ClassB** ppWorker)
{
    (*ppWorker) = (ClassB*)new ClassC();
}

__global__ void PrintArray(ClassA * pArray)
{
    for (int i = _ONE_LINE_START_; i < _ONE_LINE_END_; ++i)
    {
        pArray[i].DebugPrint();
    }
    printf("\n");
}

__global__ void TestA1(ClassA * pArray)
{
    __pWorker->InitialArray(pArray);
}

__global__ void TestA2()
{
    ClassA array[_TESTLENGTH_];
    __pWorker->InitialArray(array);
    for (int i = _ONE_LINE_START_; i < _ONE_LINE_END_; ++i)
    {
        array[i].DebugPrint();
    }
    printf("\n");
}

int main()
{
    //=======================
    //Initial the constant
    ClassB** ppDeviceWorker;
    cudaMalloc((void**)&ppDeviceWorker, sizeof(ClassB*));
    InitialConstant << <1, 1 >> > (ppDeviceWorker);
    cudaMemcpyToSymbol(__pWorker, ppDeviceWorker, sizeof(ClassB*));
    cudaFree(ppDeviceWorker);

    //=======================
    //Test with malloc buffer
    ClassA* pDeviceBuffer;
    cudaMalloc((void**)&pDeviceBuffer, sizeof(ClassA) * _TESTLENGTH_);
    TestA1 << <1, 1 >> > (pDeviceBuffer);
    PrintArray << <1, 1 >> > (pDeviceBuffer);

    //=======================
    //Test with device array
    //TestA2 << <1, 1 >> > ();

    //=======================
    //Free
    cudaFree(pDeviceBuffer);
}
