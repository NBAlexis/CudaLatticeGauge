//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for some common CUDA usage
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CUDAHELP_H_
#define _CUDAHELP_H_

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
#include <string>
#include "cuda_runtime.h"
#include "vector_types.h"
#include "cuda.h"
#include "cuComplex.h"
#include <cufft.h>
#include <cufftXt.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

#ifdef __DRIVER_TYPES_H__
static const char *_cudaGetErrorEnum(cudaError_t error) {
    return cudaGetErrorName(error);
}
#endif

#ifdef __DRIVER_TYPES_H__
#ifndef DEVICE_RESET
#define DEVICE_RESET cudaDeviceReset();
#endif
#else
#ifndef DEVICE_RESET
#define DEVICE_RESET
#endif
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

void StartTimer(unsigned long long & uiStart)
{
    uiStart = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
}

float StopTimer(unsigned long long uiStart)
{
    return ((std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count()) - uiStart) * 0.001f;
}

__device__ static __inline__ cuComplex cuCsqrtf(const cuComplex& c)
{
    float fRadius = sqrt(c.x*c.x+c.y*c.y);
    float fCosA = c.x / fRadius;
    cuComplex out;
    out.x = sqrt(0.5f * fRadius * (fCosA + 1.0f));
    out.y = sqrt(0.5f * fRadius * (1.0f - fCosA));
    if (c.y < 0.0f)
    {
        out.y = -1.0f * out.y;
    }

    return out;
}

inline void ClearLog()
{
    std::ofstream ofs;
    ofs.open("CudaHelper.log", std::ofstream::out);
    ofs.close();
}

inline void SaveLog(const char* formmat, ...)
{
    std::ofstream ofs;
    ofs.open("CudaHelper.log", std::ofstream::out | std::ofstream::app);
    static char cBuff[4096];

    va_list arg;
    {
        va_start(arg, formmat);
        vsnprintf_s(cBuff, 4095, formmat, arg);
        va_end(arg);
    }
    ofs << cBuff;
    ofs.close();
}

#endif //#ifndef _CUDAHELP_H_

//=============================================================================
// END OF FILE
//=============================================================================