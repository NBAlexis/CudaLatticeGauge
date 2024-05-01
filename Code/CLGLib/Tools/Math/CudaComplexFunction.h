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
    __device__ __host__ static __inline__ Real __cuCabsSqf(const CLGComplex& c)
    {
        return c.x *c.x + c.y * c.y;
    }

#if !_CLG_DOUBLEFLOAT
    __device__ __host__ static __inline__ DOUBLE __cuCabsSqfd(const cuDoubleComplex& c)
    {
        return c.x * c.x + c.y * c.y;
    }

    __device__ __host__ static __inline__ DOUBLE __cuCabsSqd(const cuDoubleComplex& c)
    {
        return c.x * c.x + c.y * c.y;
    }
#endif

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
        const Real beforesqrtx = F(0.5) * fRadius * (fCosA + F(1.0));
        const Real beforesqrty = F(0.5) * fRadius * (F(1.0) - fCosA);
        out.x = beforesqrtx > _CLG_FLT_MIN ? _sqrt(beforesqrtx) : beforesqrtx;
        out.y = beforesqrty > _CLG_FLT_MIN ? _sqrt(beforesqrty) : beforesqrty;
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
#if !_CLG_DOUBLEFLOAT
    __host__ static __inline__ cuDoubleComplex cuCdivf_cd_host(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(x.x / y, x.y / y);
    }
    __device__ static __inline__ cuDoubleComplex cuCdivf_cd(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(__div(x.x, y), __div(x.y, y));
    }
#endif

    __host__ static __inline__ CLGComplex cuCdivf_cr_host(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x / y, x.y / y);
    }

    __device__ __host__ static __inline__ CLGComplex cuCmulf_cr(const CLGComplex& x, Real y)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }

#if !_CLG_DOUBLEFLOAT

    __device__ __host__ static __inline__ cuDoubleComplex cuCmulf_cd(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(x.x * y, x.y * y);
    }
#endif

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

__host__ __device__ __inline__
UBOOL operator==(const cuComplex& c1, const cuComplex& c2)
{
    return c1.x == c2.x && c1.y == c2.y;
}

__host__ __device__ __inline__
UBOOL operator==(const cuDoubleComplex& c1, const cuDoubleComplex& c2)
{
    return c1.x == c2.x && c1.y == c2.y;
}

__host__ __device__ __inline__
UBOOL operator!=(const cuComplex& c1, const cuComplex& c2)
{
    return !(c1 == c2);
}

__host__ __device__ __inline__
UBOOL operator!=(const cuDoubleComplex& c1, const cuDoubleComplex& c2)
{
    return !(c1 == c2);
}

#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _CUDACOMPLEXFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
