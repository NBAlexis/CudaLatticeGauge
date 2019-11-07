//=============================================================================
// FILENAME : CudaComplexFunction.h
// 
// DESCRIPTION:
// Add some function to CLGComplex where it does not have yet
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
    __device__ __host__ static __inline__ Real __cuCargf(const CLGComplex& c)
    {
        return _atan2(_cuCimagf(c), _cuCrealf(c));
    }

    /**
    * |c|^2
    */
    __device__ static __inline__ Real __cuCabsSqf(const CLGComplex& c)
    {
        return c.x *c.x + c.y * c.y;
    }

    /**
    * c^p
    */
    __device__ static __inline__ CLGComplex __cuCpowerf(const CLGComplex& c, Real p)
    {
        const Real fArg = __cuCargf(c) * p;
        const Real fAbs = _pow(_cuCabsf(c), p);
        return _make_cuComplex(_cos(fArg) * fAbs, _sin(fArg) * fAbs);
    }

    /**
    * exp(c)
    */
    __device__ static __inline__ CLGComplex __cuCexpf(const CLGComplex& c)
    {
        const Real factor = _exp(c.x);
        return _make_cuComplex(factor * _cos(c.y), factor * _sin(c.y));
    }

    /**
     * log(c)
     */
    __device__ static __inline__ CLGComplex __cuClogf(const CLGComplex& c)
    {
        const Real fArg = __cuCargf(c);
        return _make_cuComplex(_log(_cuCabsf(c)), fArg > PI ? fArg - PI2 : fArg);
    }

    /**
    * _sqrt(c)
    */
    __device__ static __inline__ CLGComplex __cuCsqrtf(const CLGComplex& c)
    {
        const Real fRadius = _cuCabsf(c);
        const Real fCosA = __div(c.x, fRadius);
        CLGComplex out;
        out.x = _sqrt(F(0.5) * fRadius * (fCosA + F(1.0)));
        out.y = _sqrt(F(0.5) * fRadius * (F(1.0) - fCosA));
        // signbit should be false if x.y is negative
        //if (signbit(c.y))
        //    out.y *= -F(1.0);
        if (c.y < F(0.0)) //same as Mathematica
            out.y *= -F(1.0);

        return out;
    }

    __device__ static __inline__ CLGComplex cuCaddf_cr(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x + y, x.y);
    }
    __device__ static __inline__ CLGComplex cuCaddf_rc(Real y, const CLGComplex& x)
    {
        return _make_cuComplex(x.x + y, x.y);
    }

    __device__ static __inline__ CLGComplex cuCdivf_cr(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(__div(x.x, y), __div(x.y, y));
    }

    __host__ static __inline__ CLGComplex cuCdivf_cr_host(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x / y, x.y / y);
    }

    __device__ __host__ static __inline__ CLGComplex cuCmulf_cr(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }
    __device__ static __inline__ CLGComplex cuCmulf_rc(Real y, const CLGComplex& x)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }

    __device__ static __inline__ CLGComplex cuCsubf_cr(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x - y, x.y);
    }
    __device__ static __inline__ CLGComplex cuCsubf_rc(Real y, const CLGComplex& x)
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
