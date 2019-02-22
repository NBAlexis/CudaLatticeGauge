#include "CudaHelper.h"

class ClassA
{
public:
    __device__ ClassA() {  }
    int a;
};

__constant__ ClassA a[1];

int main()
{
    return 0;
}
