#include "CudaHelper.h"

int main()
{
    float* d1;

    size_t free, total;

    printf("\n");

    cudaMemGetInfo(&free, &total);

    printf("%d KB free of total %d KB\n", free / 1024, total / 1024);

    checkCudaErrors(cudaMalloc((void**)&d1, sizeof(float) * (1 << 20) * 1024));

    printf("\n");

    cudaMemGetInfo(&free, &total);

    printf("%d KB free of total %d KB\n", free / 1024, total / 1024);

    checkCudaErrors(cudaFree(d1));

    return 0;
}
