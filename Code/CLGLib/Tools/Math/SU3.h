//=============================================================================
// FILENAME : SU3.h
// 
// DESCRIPTION:
// This is helper functions for calculate SU3 matrix
//
// The SU3 Matrix is
// 0 1 2
// 3 4 5
// 6 7 8
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _SU3_H_
#define _SU3_H_

#define __LINE_MUL(a, b, c, d, ee, ff) cuCaddf(cuCaddf(cuCmulf(left[a], right[d]), cuCmulf(left[b], right[ee])), cuCmulf(left[c], right[ff]))

__BEGIN_NAMESPACE

/**
* res = left * right
*/
__device__ static __inline__ void _deviceSU3Multply(cuComplex* res, const cuComplex* left, const cuComplex * right)
{
    res[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
    res[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
    res[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

    res[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
    res[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
    res[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

    res[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
    res[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
    res[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
}

/**
* left = left * right
*/
__device__ static __inline__ void _deviceSU3Multply(cuComplex* left, const cuComplex * right)
{
    cuComplex res[9];
    res[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
    res[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
    res[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

    res[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
    res[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
    res[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

    res[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
    res[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
    res[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
    checkCudaErrors(cudaMemcpy(left, res, sizeof(cuComplex) * 9, cudaMemcpyDeviceToDevice));
}

/**
* Re[Tr[U]]
*/
__device__ static __inline__ FLOAT _deviceSU3ReTr(const cuComplex * u)
{
    return cuCrealf(cuCmulf(cuCmulf(u[0], u[4]), u[8]));
}

/**
* Conjugate[Transpose[U]]
*/
__device__ static __inline__ FLOAT _deviceSU3Dagger(cuComplex * u)
{
    cuComplex res[9];
    res[0] = cuConjf(u[0]);
    res[1] = cuConjf(u[3]);
    res[2] = cuConjf(u[6]);

    res[3] = cuConjf(u[1]);
    res[4] = cuConjf(u[4]);
    res[5] = cuConjf(u[7]);

    res[6] = cuConjf(u[2]);
    res[7] = cuConjf(u[5]);
    res[8] = cuConjf(u[8]);

    checkCudaErrors(cudaMemcpy(u, res, sizeof(cuComplex) * 9, cudaMemcpyDeviceToDevice));
}

__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
