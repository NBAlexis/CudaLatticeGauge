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

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

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
    __device__ static __inline__ _Complex cuCFaddI(const _Complex& a, const cuComplexI& b)    { return _make_cuComplex(a.x + b.x, a.y + b.y); }
    __device__ static __inline__ _Complex cuCIaddF(const cuComplexI& a, const _Complex& b)    { return _make_cuComplex(a.x + b.x, a.y + b.y); }

    __device__ static __inline__ cuComplexI cuCsubI(const cuComplexI& a, const cuComplexI& b) { return make_cuComplexI(a.x - b.x, a.y - b.y); }
    __device__ static __inline__ _Complex cuCFsubI(const _Complex& a, const cuComplexI& b)    { return _make_cuComplex(a.x - b.x, a.y - b.y); }
    __device__ static __inline__ _Complex cuCIsubF(const cuComplexI& a, const _Complex& b)    { return _make_cuComplex(a.x - b.x, a.y - b.y); }

    __device__ static __inline__ cuComplexI cuCmulI(const cuComplexI& a, const cuComplexI& b) { return make_cuComplexI(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
    __device__ static __inline__ _Complex cuCFmulI(const _Complex& a, const cuComplexI& b)    { return _make_cuComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
    __device__ static __inline__ _Complex cuCImulF(const cuComplexI& a, const _Complex& b)    { return _make_cuComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

    __device__ static __inline__ Real cuCabsI(const cuComplexI& x) { return _sqrt((Real)(x.x * x.x + x.y * x.y)); }

    __device__ static __inline__ _Complex cuCItoF(const cuComplexI& x) { return _make_cuComplex((Real)x.x, (Real)x.y); }

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _CUCOMPLEXI_H_

//=============================================================================
// END OF FILE
//=============================================================================
