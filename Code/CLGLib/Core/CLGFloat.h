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
#define _hostexp exp
#define _hostsqrt sqrt

#define _make_cuComplex make_cuDoubleComplex
#define _cuCaddf cuCadd
#define _cuCmulf cuCmul
#define _cuCsubf cuCsub
#define _cuConjf cuConj
#define _cuCrealf cuCreal
#define _cuCimagf cuCimag
#define _cuCabsf cuCimag
#define F(v) v

#else

//the __function is Intrinsic Functions which can be only used in device
#define _sqrt sqrtf
#define _log logf
#define _exp expf
#define _pow powf
#define _sin sinf
#define _cos cosf
#define _hostexp expf
#define _hostsqrt sqrtf

#define _make_cuComplex make_cuComplex
#define _cuCaddf cuCaddf
#define _cuCmulf cuCmulf
#define _cuCsubf cuCsubf
#define _cuConjf cuConjf
#define _cuCrealf cuCrealf
#define _cuCimagf cuCimagf
#define _cuCabsf cuCabsf
#define F(v) v##f

#endif

__BEGIN_NAMESPACE

#if _CLG_DOUBLEFLOAT

typedef double Real;
typedef cuDoubleComplex _Complex;

#else

typedef float Real;
typedef cuComplex _Complex;

#endif

__END_NAMESPACE

#endif//#ifndef _CLGFLOAT_H_

//=============================================================================
// END OF FILE
//=============================================================================