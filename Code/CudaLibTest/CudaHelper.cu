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

class A1
{
public:
    A1() {}
    ~A1() { printf("~A1"); }
};

class A2
{
public:
    A2() {}
    virtual ~A2() { printf("~A2"); }
};

class B1 : public A1
{
public:
    B1() {}
    ~B1() { printf("~B1"); }
};

class B2 : public A2
{
public:
    B2() {}
    ~B2() { printf("~B2"); }
};

class C1 : public B1
{
public:
    C1() {}
    ~C1() { printf("~C1"); }
};

class C2 : public B2
{
public:
    C2() {}
    ~C2() { printf("~C2"); }
};

int main()
{
    A1 * pB1 = (A1*)new C1();
    delete pB1;

    A2 * pB2 = (A2*)new C2();
    delete pB2;

    B1 * pB3 = (B1*)new C1();
    delete pB3;

    B2 * pB4 = (B2*)new C2();
    delete pB4;
}