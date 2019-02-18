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
#include <chrono>
#include "cuda_runtime.h"
#include "vector_types.h"

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#define _kLoop (1 << 8)
#define _kThreadPerBlock (1 << 10)

//#define _kLoop 1
//#define _kThreadPerBlock 32

__global__ void _kernelInitial(float* arr)
{
    arr[threadIdx.x + blockIdx.x * blockDim.x] = cosf(threadIdx.x + blockIdx.x * blockDim.x);
}

__global__ void _kernelTestSum1(float* arr, float x)
{
    arr[threadIdx.x + blockIdx.x * blockDim.x] += cosf(x);
}

__global__ void _kernelTestSum2(float* arr, float x)
{
    for (int i = 0; i < _kLoop; ++i)
    {
        arr[threadIdx.x * _kLoop + i] += cosf(x);
    }
}

__global__ void _kernelPrint(float* arr)
{
    for (int i = 123; i < 128; ++i)
    {
        printf("[%d]=%f,", i, arr[i]);
    }
}

__global__ void _kernelReduceFloat(float* arr, int iJump, int iMax)
{
    //for length 16 array
    //for jump = 1, this is 1->0, 3->2, 5->4, 7->6, 9->10, 11->10, 13->12, 15->14 
    //for jump = 2, this is 2->0, 6->4, 10->8, 14->12 
    //for jump = 4, this is 4->0, 12->8 
    //for jump = 8, this is 8->0, and is finished.

    //id target = idx * (jump << 1)
    //id from = target + jump
    int iIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (iJump << 1) + iJump;
    if (iIdFrom < iMax)
    {
        arr[iIdFrom - iJump] += arr[iIdFrom];
    }
}

int FindLargestPowerOf2(int iIn)
{
    int iRet = 0;
    while ((1 << iRet) < iIn)
    {
        ++iRet;
    }
    return iRet;
}

float Reduce(float* deviceBuffer, int iLength)
{
    int iRequiredDim = (iLength + 1) >> 1;
    int iPower = FindLargestPowerOf2(iRequiredDim);
    for (int i = 0; i <= iPower; ++i)
    {
        int iJump = 1 << i;
        int iThreadNeeded = 1 << (iPower - i);
        int iBlock = iThreadNeeded > 1024 ? iThreadNeeded / 1024 : 1;
        int iThread = iThreadNeeded > 1024 ? 1024 : iThreadNeeded;
        _kernelReduceFloat << <iBlock, iThread >> > (deviceBuffer, iJump, iLength);
    }
    float result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(float) * 1, cudaMemcpyDeviceToHost);
    return result[0];
}

int main()
{
    float* d1;
    float* d2;
    cudaMalloc((void**)&d1, sizeof(float) * _kLoop * _kThreadPerBlock);
    cudaMalloc((void**)&d2, sizeof(float) * _kLoop * _kThreadPerBlock);
    _kernelInitial << <_kLoop, _kThreadPerBlock >> > (d1);
    _kernelInitial << <_kLoop, _kThreadPerBlock >> > (d2);

    UINT t1 = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    float res1 = Reduce(d2, _kLoop * _kThreadPerBlock);

    cudaDeviceSynchronize();
    UINT t2 = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    thrust::device_ptr<float> dp(d1);
    thrust::device_vector<float> d_x(dp, dp + _kLoop * _kThreadPerBlock);
    float res2 = thrust::reduce(d_x.begin(), d_x.end(), 0.0f, thrust::plus<float>());

    cudaDeviceSynchronize();
    UINT t3 = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    printf("res1=%f\n", res1);
    printf("res2=%f\n", res2);

    printf("\n t1 = %d, t2 = %d\n", t2 - t1, t3 - t2);

    cudaFree(d1);
    cudaFree(d2);

    return 0;
}
