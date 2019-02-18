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
        return c.x *c.x + c.y * c.y;
    }

    /**
    * c^p
    */
    __device__ static __inline__ _Complex __cuCpowerf(const _Complex& c, Real p)
    {
        Real fArg = __cuCargf(c) * p;
        Real fAbs = _pow(_cuCabsf(c), p);
        return _make_cuComplex(_cos(fArg) * fAbs, _sin(fArg) * fAbs);
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
        Real fCosA = __div(c.x, _cuCabsf(c));
        _Complex out;
        out.x = _sqrt(F(0.5) * fRadius * (fCosA + F(1.0)));
        out.y = _sqrt(F(0.5) * fRadius * (F(1.0) - fCosA));
        // signbit should be false if x.y is negative
        if (signbit(c.y))
            out.y *= -F(1.0);

        return out;
    }

    __device__ static __inline__ _Complex cuCaddf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(x.x + y, x.y);
    }
    __device__ static __inline__ _Complex cuCaddf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(x.x + y, x.y);
    }

    __device__ static __inline__ _Complex cuCdivf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(__div(x.x, y), __div(x.y, y));
    }

    __device__ static __inline__ _Complex cuCmulf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }
    __device__ static __inline__ _Complex cuCmulf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }

    __device__ static __inline__ _Complex cuCsubf_cr(const _Complex& x, Real y)
    {
        return _make_cuComplex(x.x - y, x.y);
    }
    __device__ static __inline__ _Complex cuCsubf_rc(Real y, const _Complex& x)
    {
        return _make_cuComplex(y - x.x, x.y);
    }

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _CUDACOMPLEXFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
