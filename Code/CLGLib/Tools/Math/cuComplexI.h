//=============================================================================
// FILENAME : cuComplexI.h
// 
// DESCRIPTION:
// This is a complex with integers
//
//
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================

#ifndef _CUCOMPLEXI_H_
#define _CUCOMPLEXI_H_

__BEGIN_NAMESPACE

typedef int2 cuComplexI;

__device__ static __inline__ INT cuCrealI(const cuComplexI& x) { return x.x; }

__device__ static __inline__ INT cuCimagI(const cuComplexI& x) { return x.y; }

__device__ static __inline__ cuComplexI make_cuComplexI(INT r, INT i)
{
    cuComplexI res;
    res.x = r;
    res.y = i;
    return res;
}

__device__ static __inline__ cuComplexI cuConjI(const cuComplexI& x) { return make_cuComplexI(cuCrealI(x), -cuCimagI(x)); }

__device__ static __inline__ cuComplexI cuCaddI(const cuComplexI& a, const cuComplexI& b) { return make_cuComplexI(a.x + b.x, a.y + b.y); }
__device__ static __inline__ cuComplex cuCaddI(const cuComplex& a, const cuComplexI& b) { return make_cuComplex(a.x + b.x, a.y + b.y); }
__device__ static __inline__ cuComplex cuCaddI(const cuComplexI& a, const cuComplex& b) { return make_cuComplex(a.x + b.x, a.y + b.y); }

__device__ static __inline__ cuComplexI cuCsubI(const cuComplexI& a, const cuComplexI& b) { return make_cuComplexI(a.x - b.x, a.y - b.y); }
__device__ static __inline__ cuComplex cuCsubI(const cuComplex& a, const cuComplexI& b) { return make_cuComplex(a.x - b.x, a.y - b.y); }
__device__ static __inline__ cuComplex cuCsubI(const cuComplexI& a, const cuComplex& b) { return make_cuComplex(a.x - b.x, a.y - b.y); }

__device__ static __inline__ cuComplexI cuCmulI(const cuComplexI& a, const cuComplexI& b) { return make_cuComplexI(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
__device__ static __inline__ cuComplex cuCmulI(const cuComplex& a, const cuComplexI& b) { return make_cuComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
__device__ static __inline__ cuComplex cuCmulI(const cuComplexI& a, const cuComplex& b) { return make_cuComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

__device__ static __inline__ FLOAT cuCabsI(const cuComplexI& x) { return sqrt(x.x * x.x + x.y * x.y); }

__device__ static __inline__ cuComplex cuCItoF(const cuComplexI& x) { return make_cuComplex(x.x, x.y); }

__END_NAMESPACE

#endif //#ifndef _CUCOMPLEXI_H_

//=============================================================================
// END OF FILE
//=============================================================================
