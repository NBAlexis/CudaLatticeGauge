//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for some common CUDA usage
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

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

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

#ifndef _CUDAHELP_H_
#define _CUDAHELP_H_

template <typename T> void check(T result, char const *const func, const char *const file,
    int const line) {
    if (result) {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
            static_cast<unsigned int>(result), cudaGetErrorName(result), func);
        cudaDeviceReset();
            // Make sure we call CUDA Device Reset before exiting
            exit(EXIT_FAILURE);
    }
}

#endif //#ifndef _CUDAHELP_H_

//=============================================================================
// END OF FILE
//=============================================================================