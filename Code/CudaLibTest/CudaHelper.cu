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

class ClassA
{
public:
    enum { kArrayLength = 20, };
    ClassA()
    {
        memset(m_pPtr, 0, sizeof(int)* kArrayLength);
    }
    __device__ void Add(int v)
    {
        m_pPtr[v] = m_pPtr[v] + v;
    }
    __device__ void DebugPrint() const
    {
        for (int i = 0; i < kArrayLength; ++i)
            printf("%d,", m_pPtr[i]);
        printf("\n");
    }
    int m_pPtr[kArrayLength];
};

__global__ void _kernelAdd(ClassA a)
{
    a.Add(threadIdx.x);
    __syncthreads();
    a.DebugPrint();
}

int main()
{
    ClassA a;
    _kernelAdd << <1, ClassA::kArrayLength >> > (a);
    return 0;
}
