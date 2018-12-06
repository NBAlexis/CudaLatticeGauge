//=============================================================================
// FILENAME : CudaComplexFunction.h
// 
// DESCRIPTION:
// Add some function to cuComplex where it does not have yet
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================

#ifndef _CUDACOMPLEXFUNCTION_H_
#define _CUDACOMPLEXFUNCTION_H_

__BEGIN_NAMESPACE

/**
* arg(c)
*/
__device__ static __inline__ FLOAT cuCargf(const cuComplex& c)
{
    return atan2f(cuCimagf(c), cuCrealf(c));
}

/**
* |c|^2
*/
__device__ static __inline__ FLOAT cuCabsSqf(const cuComplex& c)
{
    return c.x * c.x + c.y * c.y;
}

/**
* c^p
*/
__device__ static __inline__ cuComplex cuCpowerf(const cuComplex& c, FLOAT p)
{
    FLOAT fArg = cuCargf(c) * p;
    FLOAT fAbs = powf(cuCabsf(c), p);
    return make_cuComplex(cosf(fArg) * fAbs, sinf(fArg) * fAbs);
}

/**
* exp(c)
*/
__device__ static __inline__ cuComplex cuCexpf(const cuComplex& c)
{
    FLOAT factor = exp(c.x);
    return make_cuComplex(factor * cosf(c.y), factor * sinf(c.y));
}

/**
* sqrt(c)
*/
__device__ static __inline__ cuComplex cuCsqrt(const cuComplex& c)
{
    FLOAT fRadius = cuCabsf(c);
    FLOAT fCosA = c.x / cuCabsf(c);
    cuComplex out;
    out.x = sqrtf(fRadius * (fCosA + 1.0f) / 2.0f);
    out.y = sqrtf(fRadius * (1.0f - fCosA) / 2.0f);
    // signbit should be false if x.y is negative
    if (signbit(c.y))
        out.y *= -1.0f;

    return out;
}

__device__ static __inline__ cuComplex cuCaddf(const cuComplex& x, FLOAT y)
{
    return make_cuComplex(cuCrealf(x) + y, cuCimagf(x));
}
__device__ static __inline__ cuComplex cuCdivf(const cuComplex& x, FLOAT y)
{
    return make_cuComplex(cuCrealf(x) / y, cuCimagf(x) / y);
}
__device__ static __inline__ cuComplex cuCmulf(const cuComplex& x, FLOAT y)
{
    return make_cuComplex(cuCrealf(x) * y, cuCimagf(x) * y);
}
__device__ static __inline__ cuComplex cuCsubf(const cuComplex& x, FLOAT y)
{
    return make_cuComplex(cuCrealf(x) - y, cuCimagf(x));
}

__END_NAMESPACE

#endif //#ifndef _CUDACOMPLEXFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
