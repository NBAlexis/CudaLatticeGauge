//=============================================================================
// FILENAME : CudaComplexFunction.h
// 
// DESCRIPTION:
// Add some function to _Complex where it does not have yet
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
__device__ static __inline__ Real cuCargf(const _Complex& c)
{
    return atan2f(_cuCimagf(c), _cuCrealf(c));
}

/**
* |c|^2
*/
__device__ static __inline__ Real cuCabsSqf(const _Complex& c)
{
    return c.x * c.x + c.y * c.y;
}

/**
* c^p
*/
__device__ static __inline__ _Complex cuCpowerf(const _Complex& c, Real p)
{
    Real fArg = cuCargf(c) * p;
    Real fAbs = _pow(_cuCabsf(c), p);
    return _make_cuComplex(_cos(fArg) * fAbs, _sin(fArg) * fAbs);
}

/**
* exp(c)
*/
__device__ static __inline__ _Complex cuCexpf(const _Complex& c)
{
    Real factor = _exp(c.x);
    return _make_cuComplex(factor * _cos(c.y), factor * _sin(c.y));
}

/**
* _sqrt(c)
*/
__device__ static __inline__ _Complex cuCsqrtf(const _Complex& c)
{
    Real fRadius = _cuCabsf(c);
    Real fCosA = c.x / _cuCabsf(c);
    _Complex out;
    out.x = _sqrt(fRadius * (fCosA + 1.0f) / 2.0f);
    out.y = _sqrt(fRadius * (1.0f - fCosA) / 2.0f);
    // signbit should be false if x.y is negative
    if (signbit(c.y))
        out.y *= -1.0f;

    return out;
}

__device__ static __inline__ _Complex cuCaddf(const _Complex& x, Real y)
{
    return _make_cuComplex(_cuCrealf(x) + y, _cuCimagf(x));
}
__device__ static __inline__ _Complex cuCaddf(Real y, const _Complex& x)
{
    return _make_cuComplex(_cuCrealf(x) + y, _cuCimagf(x));
}

__device__ static __inline__ _Complex cuCdivf(const _Complex& x, Real y)
{
    return _make_cuComplex(_cuCrealf(x) / y, _cuCimagf(x) / y);
}

__device__ static __inline__ _Complex cuCmulf(const _Complex& x, Real y)
{
    return _make_cuComplex(_cuCrealf(x) * y, _cuCimagf(x) * y);
}
__device__ static __inline__ _Complex cuCmulf(Real y, const _Complex& x)
{
    return _make_cuComplex(_cuCrealf(x) * y, _cuCimagf(x) * y);
}

__device__ static __inline__ _Complex cuCsubf(const _Complex& x, Real y)
{
    return _make_cuComplex(_cuCrealf(x) - y, _cuCimagf(x));
}
__device__ static __inline__ _Complex cuCsubf(Real y, const _Complex& x)
{
    return _make_cuComplex(y - _cuCrealf(x), -_cuCimagf(x));
}

__END_NAMESPACE

#endif //#ifndef _CUDACOMPLEXFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
