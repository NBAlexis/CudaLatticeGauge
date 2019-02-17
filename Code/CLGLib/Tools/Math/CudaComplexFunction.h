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

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    /**
    * arg(c)
    */
    __device__ static __inline__ Real __cuCargf(const _Complex& c)
    {
        return atan2f(_cuCimagf(c), _cuCrealf(c));
    }

    /**
    * |c|^2
    */
    __device__ static __inline__ Real __cuCabsSqf(const _Complex& c)
    {
        return __add(__mul(c.x, c.x), __mul(c.y, c.y));
    }

    /**
    * c^p
    */
    __device__ static __inline__ _Complex __cuCpowerf(const _Complex& c, Real p)
    {
        Real fArg = __mul(__cuCargf(c), p);
        Real fAbs = _pow(_cuCabsf(c), p);
        return _make_cuComplex(__mul(_cos(fArg), fAbs), __mul(_sin(fArg), fAbs));
    }

    /**
    * exp(c)
    */
    __device__ static __inline__ _Complex __cuCexpf(const _Complex& c)
    {
        Real factor = _exp(c.x);
        return _make_cuComplex(factor * _cos(c.y), factor * _sin(c.y));
    }

    /**
    * _sqrt(c)
    */
    __device__ static __inline__ _Complex __cuCsqrtf(const _Complex& c)
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

    __device__ static __inline__ _Complex cuCaddf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(__add(x.x, y), x.y);
    }
    __device__ static __inline__ _Complex cuCaddf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(__add(x.x, y), _cuCimagf(x));
    }

    __device__ static __inline__ _Complex cuCdivf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(__div(x.x, y), __div(x.y, y));
    }

    __device__ static __inline__ _Complex cuCmulf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(__mul(x.x, y), __mul(x.y, y));
    }
    __device__ static __inline__ _Complex cuCmulf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(__mul(x.x, y), __mul(x.y, y));
    }

    __device__ static __inline__ _Complex cuCsubf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(__sub(x.x, y), x.y);
    }
    __device__ static __inline__ _Complex cuCsubf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(__sub(y, x.x), x.y);
    }

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _CUDACOMPLEXFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
