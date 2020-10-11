//=============================================================================
// FILENAME : CLGFloat.h
// 
// DESCRIPTION:
// Add precsion for floats
//
// REVISION:
//  [12/15/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CLGFLOAT_H_
#define _CLGFLOAT_H_

#if _CLG_DOUBLEFLOAT

#define _sqrt sqrt
#define _log log
#define _exp exp
#define _pow pow
#define _sin sin
#define _cos cos
#define __div(a, b) ((a) / (b))
#define __rcp(a) (F(1.0) / (a))
#define _hostlog log
#define _hostlog10 log10
#define _hostexp exp
#define _hostsqrt sqrt

#define _atan2 atan2
#define _make_cuComplex make_cuDoubleComplex
#define _cuCaddf cuCadd
#define _cuCmulf cuCmul
#define _cuCsubf cuCsub
#define _cuConjf cuConj
#define _cuCrealf cuCreal
#define _cuCimagf cuCimag
#define _cuCabsf cuCabs
#define _cuCdivf cuCdiv
#define F(v) v
#if defined(__cplusplus) && defined(__CUDACC__)
#define _floor2int __double2int_rd
#define _round2int __double2int_rn
#else
#define _floor2int(a) static_cast<INT>(floor(a))
#define _round2int(a) static_cast<INT>(round(a))
#endif

#else

#if defined(__cplusplus) && defined(__CUDACC__)
#define _sqrt __fsqrt_rn
#define _log __logf
#define _exp __expf
#define _pow __powf
#define _sin __sinf
#define _cos __cosf
#define __div __fdividef
#define __rcp __frcp_rn
#else
//the __function is Intrinsic Functions which can be only used in device
#define _sqrt sqrtf
#define _sqrtd sqrt
#define _log logf
#define _exp expf
#define _pow powf
#define _sin sinf
#define _cos cosf
#define __div(a, b) ((a) / (b))
#define __rcp(a) (F(1.0) / (a))
#endif

#define _hostlog logf
#define _hostlog10 log10f
#define _hostlogd log
#define _hostlog10d log10
#define _hostexp expf
#define _hostsqrt sqrtf
#define _hostsqrtd sqrt
#define _hostexpd exp

#define _atan2 atan2f
#define _make_cuComplex make_cuComplex
#define _cuCaddf cuCaddf
#define _cuCmulf cuCmulf
#define _cuCsubf cuCsubf
#define _cuConjf cuConjf
#define _cuCrealf cuCrealf
#define _cuCimagf cuCimagf
#define _cuCabsf cuCabsf
#define _cuCdivf cuCdivf
#define F(v) v##f
#if defined(__cplusplus) && defined(__CUDACC__)
#define _floor2int __float2int_rd
#define _round2int __float2int_rn
#else
#define _floor2int(a) static_cast<INT>(floor(a))
#define _round2int(a) static_cast<INT>(round(a))
#endif

#endif

#if _CLG_DOUBLEFLOAT

#define _CLG_FLT_MIN_ 1E-50   //When smaller than this, sqrt, divide will become nan

#define _CLG_FLT_DECIMAL_DIG  17                      // # of decimal digits of rounding precision
#define _CLG_FLT_DIG          15                      // # of decimal digits of precision
#define _CLG_FLT_EPSILON      2.2204460492503131e-016 // smallest such that 1.0+DBL_EPSILON != 1.0
#define _CLG_FLT_HAS_SUBNORM  1                       // type does support subnormal numbers
#define _CLG_FLT_MANT_DIG     53                      // # of bits in mantissa
#define _CLG_FLT_MAX          1.7976931348623158e+308 // max value
#define _CLG_FLT_MAX_10_EXP   308                     // max decimal exponent
#define _CLG_FLT_MAX_EXP      1024                    // max binary exponent
#define _CLG_FLT_MIN          2.2250738585072014e-308 // min positive value
#define _CLG_FLT_MIN_10_EXP   (-307)                  // min decimal exponent
#define _CLG_FLT_MIN_EXP      (-1021)                 // min binary exponent
#define _CLG_FLT_RADIX        2                       // exponent radix
#define _CLG_FLT_TRUE_MIN     4.9406564584124654e-324 // min positive value

#else

//They are not defined in GCC, so we define them explicitly
#define _CLG_FLT_MIN_ 1E-22F   //When smaller than this, sqrt, divide will become nan

#define _CLG_FLT_DECIMAL_DIG  9                       // # of decimal digits of rounding precision
#define _CLG_FLT_DIG          6                       // # of decimal digits of precision
#define _CLG_FLT_EPSILON      1.192092896e-07F        // smallest such that 1.0+FLT_EPSILON != 1.0
#define _CLG_FLT_HAS_SUBNORM  1                       // type does support subnormal numbers
#define _CLG_FLT_GUARD        0
#define _CLG_FLT_MANT_DIG     24                      // # of bits in mantissa
#define _CLG_FLT_MAX          3.402823466e+38F        // max value
#define _CLG_FLT_MAX_10_EXP   38                      // max decimal exponent
#define _CLG_FLT_MAX_EXP      128                     // max binary exponent
#define _CLG_FLT_MIN          1.175494351e-38F        // min normalized positive value
#define _CLG_FLT_MIN_10_EXP   (-37)                   // min decimal exponent
#define _CLG_FLT_MIN_EXP      (-125)                  // min binary exponent
#define _CLG_FLT_NORMALIZE    0
#define _CLG_FLT_RADIX        2                       // exponent radix
#define _CLG_FLT_TRUE_MIN     1.401298464e-45F        // min positive value

#endif

#define _zeroc (_make_cuComplex(F(0.0), F(0.0)))
#define _onec (_make_cuComplex(F(1.0), F(0.0)))
#define _imgc (_make_cuComplex(F(0.0), F(1.0)))
//Those are constants we are using

//save some constant memory of cuda?
#define PI (F(3.141592653589))
#define PISQ (F(9.8696044010893586188344909998761511353137))

// - 1/4294967296UL
#define AM (F(0.00000000023283064365386963))
// - _sqrt(2)
#define SQRT2 (F(1.4142135623730951))
// - 1 / _sqrt(2), or _sqrt(2)/2
#define InvSqrt2 (F(0.7071067811865475))
// - 2.0f * PI
#define PI2 (F(6.283185307179586))

// 1.0f / _sqrt(3)
#define InvSqrt3 (F(0.5773502691896258))
// 2.0f / _sqrt(3)
#define InvSqrt3_2 (F(1.1547005383792517))

#define OneOver3 (F(0.33333333333333333333333333333333))
#define OneOver6 (F(0.16666666666666666666666666666667))
#define OneOver12 (F(0.0833333333333333333333333333333))
#define OneOver24 (F(0.04166666666666666666666666666667))
#define OneOver48 (F(0.02083333333333333333333333333333))
#define FiveOver3 (F(1.66666666666666666666666666666666))

//typically, 0.3-0.5 - arXiv:002.4232
#define OmelyanLambda2 (F(0.38636665500756728))

//11 Stage, from 10.1016/s0010-4655(02)00754-3
//Velocity version, Variant 8
#define _11Stage_Rho (F(0.2539785108410595))
#define _11Stage_Theta (F(-0.03230286765269967))
#define _11Stage_VarTheta (F(0.08398315262876693))
#define _11Stage_Lambda (F(0.6822365335719091))

#define _11Stage_VarThetaP (F(0.2750081212332419))
#define _11Stage_LambdaP (F(-0.1347950099106792))
#define _11Stage_RhoP (F(-0.08442961950707149))
#define _11Stage_ThetaP (F(0.3549000571574260))

static __host__ __device__ __inline__ cuDoubleComplex _cToDouble(const cuComplex& c)
{
    return make_cuDoubleComplex(c.x, c.y);
}
static __host__ __device__ __inline__ cuComplex _cToFloat(const cuDoubleComplex& c)
{
    return make_cuComplex(static_cast<FLOAT>(c.x), static_cast<FLOAT>(c.y));
}

__BEGIN_NAMESPACE

//NOTE, _Complex is already a keyword in GCC
#if _CLG_DOUBLEFLOAT

typedef double Real;
typedef cuDoubleComplex CLGComplex;

#else

typedef float Real;
typedef cuComplex CLGComplex;

#endif

__END_NAMESPACE

#endif//#ifndef _CLGFLOAT_H_

//=============================================================================
// END OF FILE
//=============================================================================