//=============================================================================
// FILENAME : DeviceInlineTemplate.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [mm/dd/yy]
//  [07/03/2024 nbale]
//=============================================================================
#include "DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "Tools/Math/DeviceInlineTemplate2.h"
#include "Tools/Math/ZN.h"
#include "Tools/Math/DN.h"

#ifndef _DEVICEINLINETEMPLATE_H_
#define _DEVICEINLINETEMPLATE_H_

__BEGIN_NAMESPACE

#pragma region element-wise

template<typename T> __device__ __inline__ T _makeId() = delete;

template<> __device__ __inline__ INT _makeId<INT>() { return 1; }
template<> __device__ __inline__ UINT _makeId<UINT>() { return 1; }
template<> __device__ __inline__ BYTE _makeId<BYTE>() { return 1; }
template<> __device__ __inline__ Real _makeId<Real>() { return F(1.0); }
template<> __device__ __inline__ CLGComplex _makeId<CLGComplex>() { return _onec; }

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ FLOAT _makeId<FLOAT>() { return 1.0f; }
template<> __device__ __inline__ cuComplex _makeId<cuComplex>()
{
    return make_cuComplex(1.0f, 0.0f);
}
#else
template<> __device__ __inline__ DOUBLE _makeId<DOUBLE>() { return 1.0; }
template<> __device__ __inline__ cuDoubleComplex _makeId<cuDoubleComplex>()
{
    return make_cuDoubleComplex(1.0, 0.0);
}
#endif


template<> __device__ __inline__ deviceSU2 _makeId<deviceSU2>()
{
    return deviceSU2::makeSU2Id();
}

template<> __device__ __inline__ deviceSU3 _makeId<deviceSU3>()
{
    return deviceSU3::makeSU3Id();
}

template<> __device__ __inline__ deviceSU2Vector _makeId<deviceSU2Vector>()
{
    return deviceSU2Vector::makeOneSU2Vector();
}

template<> __device__ __inline__ deviceSU3Vector _makeId<deviceSU3Vector>()
{
    return deviceSU3Vector::makeOneSU3Vector();
}

template<typename T> __device__ __inline__ void _Id(T& v) = delete;

template<> __device__ __inline__ void _Id<INT>(INT &v) { v = 1; }
template<> __device__ __inline__ void _Id<UINT>(UINT& v) { v = 1; }
template<> __device__ __inline__ void _Id<BYTE>(BYTE& v) { v = 1; }
template<> __device__ __inline__ void _Id<Real>(Real& v) { v = F(1.0); }
template<> __device__ __inline__ void _Id<CLGComplex>(CLGComplex& v) { v.x = F(1.0); v.y = F(0.0); }
#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ void _Id<FLOAT>(FLOAT& v) { v = 1.0f; }
template<> __device__ __inline__ void _Id<cuComplex>(cuComplex& v)
{
    v.x = 1.0f; v.y = 0.0f;
}
#else
template<> __device__ __inline__ void _Id<DOUBLE>(DOUBLE& v) { v = 1.0; }
template<> __device__ __inline__ void _Id<cuDoubleComplex>(cuDoubleComplex& v)
{
    v.x = 1.0; v.y = 0.0;
}
#endif
template<> __device__ __inline__ void _Id<deviceSU2>(deviceSU2& v) { v.Id(); }
template<> __device__ __inline__ void _Id<deviceSU3>(deviceSU3& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceSUN<N, NoE>& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceSLNC<N, NoE>& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceUN<N, NoE>& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceON<N, NoE>& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceSON<N, NoE>& v) { v.Id(); }
template<INT N> __device__ __inline__ void _Id(deviceZN<N>& v) { v.Id(); }
template<> __device__ __inline__ void _Id<deviceSU2Vector>(deviceSU2Vector& v) { v.Id(); }
template<> __device__ __inline__ void _Id<deviceSU3Vector>(deviceSU3Vector& v) { v.Id(); }
template<INT N, INT NoE> __device__ __inline__ void _Id(deviceSUNVector<N, NoE>& v) { v.Id(); }
template<> __device__ __inline__ void _Id<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& v) { v.Id(); }

template<typename T> __device__ __inline__ T _makeZero() = delete;
template<typename T> __host__ __inline__ T _makeZeroHost() = delete;
template<> __host__ __inline__ FLOAT _makeZeroHost<FLOAT>() { return 0.0f; }
template<> __host__ __inline__ DOUBLE _makeZeroHost<DOUBLE>() { return 0.0; }
template<> __host__ __inline__ cuComplex _makeZeroHost<cuComplex>() { return make_cuComplex(0.0f, 0.0f); }
template<> __host__ __inline__ cuDoubleComplex _makeZeroHost<cuDoubleComplex>() { return make_cuDoubleComplex(0.0, 0.0); }

template<> __device__ __inline__ INT _makeZero<INT>() { return 0; }
template<> __device__ __inline__ UINT _makeZero<UINT>() { return 0; }
template<> __device__ __inline__ BYTE _makeZero<BYTE>() { return 0; }
template<> __device__ __inline__ Real _makeZero<Real>() { return F(0.0); }
template<> __device__ __inline__ CLGComplex _makeZero<CLGComplex>() { return _zeroc; }

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ FLOAT _makeZero<FLOAT>() { return 0.0f; }
template<> __device__ __inline__ cuComplex _makeZero<cuComplex>()
{
    return make_cuComplex(0.0f, 0.0f);
}
#else
template<> __device__ __inline__ DOUBLE _makeZero<DOUBLE>() { return 0.0; }
template<> __device__ __inline__ cuDoubleComplex _makeZero<cuDoubleComplex>()
{
    return make_cuDoubleComplex(0.0, 0.0);
}
#endif

template<> __device__ __inline__ deviceSU2 _makeZero<deviceSU2>()
{
    return deviceSU2::makeSU2Zero();
}

template<> __device__ __inline__ deviceSU3 _makeZero<deviceSU3>()
{
    return deviceSU3::makeSU3Zero();
}

template<> __device__ __inline__ deviceSU2Vector _makeZero<deviceSU2Vector>()
{
    return deviceSU2Vector::makeZeroSU2Vector();
}

template<> __device__ __inline__ deviceSU3Vector _makeZero<deviceSU3Vector>()
{
    return deviceSU3Vector::makeZeroSU3Vector();
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _makeZero<deviceWilsonVectorSU3>()
{
    return deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
}

template<typename T> __device__ __inline__ void _Zero(T& v) = delete;

template<> __device__ __inline__ void _Zero<INT>(INT& v) { v = 0; }
template<> __device__ __inline__ void _Zero<UINT>(UINT& v) { v = 0; }
template<> __device__ __inline__ void _Zero<BYTE>(BYTE& v) { v = 0; }
template<> __device__ __inline__ void _Zero<Real>(Real& v) { v = F(0.0); }
template<> __device__ __inline__ void _Zero<CLGComplex>(CLGComplex& v) { v.x = F(0.0); v.y = F(0.0); }
#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ void _Zero<FLOAT>(FLOAT& v) { v = 0.0f; }
template<> __device__ __inline__ void _Zero<cuComplex>(cuComplex& v)
{
    v.x = 0.0f; v.y = 0.0f;
}
#else
template<> __device__ __inline__ void _Zero<DOUBLE>(DOUBLE& v) { v = 0.0; }
template<> __device__ __inline__ void _Zero<cuDoubleComplex>(cuDoubleComplex& v)
{
    v.x = 0.0; v.y = 0.0;
}
#endif
template<> __device__ __inline__ void _Zero<deviceSU2>(deviceSU2& v) { v.Zero(); }
template<> __device__ __inline__ void _Zero<deviceSU3>(deviceSU3& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceSUN<N, NoE>& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceSLNC<N, NoE>& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceUN<N, NoE>& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceON<N, NoE>& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceSON<N, NoE>& v) { v.Zero(); }
template<INT N> __device__ __inline__ void _Zero(deviceZN<N>& v) { v.Zero(); }
template<> __device__ __inline__ void _Zero<deviceSU2Vector>(deviceSU2Vector& v) { v.Zero(); }
template<> __device__ __inline__ void _Zero<deviceSU3Vector>(deviceSU3Vector& v) { v.Zero(); }
template<INT N, INT NoE> __device__ __inline__ void _Zero(deviceSUNVector<N, NoE>& v) { v.Zero(); }
template<> __device__ __inline__ void _Zero<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& v) { v.Zero(); }


template<typename T> __device__ __inline__ T _makeAsK(UINT k) = delete;
template<typename T> __device__ __inline__ void _SetAsK(T& v, UINT k) = delete;

template<INT N> __device__ __inline__ void _SetAsK(deviceZN<N>& v, UINT k) { v.SetFromIndex(k); }
template<INT N> __device__ __inline__ void _SetAsK(deviceDN<N>& v, UINT k) { v.SetFromIndex(k); }


template<typename TMatrix, typename TVector> __device__ __inline__ TMatrix _makeContract(const TVector& left, const TVector& right) = delete;

template<> __device__ __inline__ Real _makeContract<Real, Real>(const Real& left, const Real& right)
{
    return left * right;
}

template<> __device__ __inline__ CLGComplex _makeContract<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(_cuConjf(left), right);
}

template<> __device__ __inline__ deviceSU2 _makeContract<deviceSU2, deviceSU2Vector>(const deviceSU2Vector& left, const deviceSU2Vector& right)
{
    return deviceSU2::makeSU2ContractV(left, right);
}

template<> __device__ __inline__ deviceSU3 _makeContract<deviceSU3, deviceSU3Vector>(const deviceSU3Vector& left, const deviceSU3Vector& right)
{
    return deviceSU3::makeSU3ContractV(left, right);
}

template<> __device__ __inline__ deviceSU3 _makeContract<deviceSU3, deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
{
    return deviceSU3::makeSU3Contract(left, right);
}

//This is to make white noise, so for gauge, it is random generator
template<typename T> __device__ __inline__ T _makeGaussian(UINT fatIdx) = delete;

template<> __device__ __inline__ Real _makeGaussian<Real>(UINT fatIdx)
{
    return _deviceRandomGaussFSqrt2(fatIdx);
}

template<> __device__ __inline__ CLGComplex _makeGaussian<CLGComplex>(UINT fatIdx)
{
    const Real r1 = _deviceRandomGaussFSqrt2(fatIdx);
    return _make_cuComplex(F(0.0), r1);
}

template<> __device__ __inline__ deviceSU2Vector _makeGaussian<deviceSU2Vector>(UINT fatIdx)
{
    return deviceSU2Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU3Vector _makeGaussian<deviceSU3Vector>(UINT fatIdx)
{
    return deviceSU3Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU2 _makeGaussian<deviceSU2>(UINT fatIdx)
{
    return deviceSU2::makeSU2RandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU3 _makeGaussian<deviceSU3>(UINT fatIdx)
{
    return deviceSU3::makeSU3RandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _makeGaussian<deviceWilsonVectorSU3>(UINT fatIdx)
{
    return deviceWilsonVectorSU3::makeRandomGaussian(fatIdx);
}

template<typename T> __device__ __inline__ T _makeSumGenerator(Real factor) = delete;

template<> __device__ __inline__ CLGComplex _makeSumGenerator(Real factor)
{
    return _make_cuComplex(F(0.0), factor);
}

template<> __device__ __inline__ deviceSU2 _makeSumGenerator(Real factor)
{
    return deviceSU2::makeSU2SumGenerator(factor);
}

template<> __device__ __inline__ deviceSU3 _makeSumGenerator(Real factor)
{
    return deviceSU3::makeSU3SumGenerator(factor);
}

template<typename T> __device__ __inline__ T _makeZ4(UINT fatIdx) = delete;

template<> __device__ __inline__ Real _makeZ4<Real>(UINT fatIdx)
{
    //not supported
    return F(1.0);
}

template<> __device__ __inline__ CLGComplex _makeZ4<CLGComplex>(UINT fatIdx)
{
    return _deviceRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU2Vector _makeZ4<deviceSU2Vector>(UINT fatIdx)
{
    return deviceSU2Vector::makeRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU3Vector _makeZ4<deviceSU3Vector>(UINT fatIdx)
{
    return deviceSU3Vector::makeRandomZ4(fatIdx);
}


template<typename T> __device__ __inline__ T _makeRandom(UINT fatIdx) = delete;

template<> __device__ __inline__ Real _makeRandom<Real>(UINT fatIdx)
{
    return _deviceRandomF(fatIdx);
}

template<> __device__ __inline__ CLGComplex _makeRandom<CLGComplex>(UINT fatIdx)
{
    const Real fArg = _deviceRandomF(fatIdx) * PI2;
    return _make_cuComplex(_cos(fArg), -_sin(fArg));
}

template<> __device__ __inline__ deviceSU2 _makeRandom<deviceSU2>(UINT fatIdx)
{
    return deviceSU2::makeSU2Random(fatIdx);
}

template<> __device__ __inline__ deviceSU3 _makeRandom<deviceSU3>(UINT fatIdx)
{
    return deviceSU3::makeSU3Random(fatIdx);
}

template<> __device__ __inline__ deviceSU2Vector _makeRandom<deviceSU2Vector>(UINT fatIdx) { return deviceSU2Vector::makeRandom(fatIdx); }
template<> __device__ __inline__ deviceSU3Vector _makeRandom<deviceSU3Vector>(UINT fatIdx) { return deviceSU3Vector::makeRandom(fatIdx); }

template<> __device__ __inline__ deviceWilsonVectorSU3 _makeRandom<deviceWilsonVectorSU3>(UINT fatIdx) { return deviceWilsonVectorSU3::makeRandom(fatIdx); }

/**
* spin index is only for Wilson Dirac fermion, for staggered fermion, put whatever you like it will has no effect
*/
template<typename T> __device__ __inline__ T _makeColorVector(BYTE spin, BYTE colorIdx) = delete;

template<> __device__ __inline__ Real _makeColorVector<Real>(BYTE spin, BYTE colorIdx)
{
    return F(1.0);
}

template<> __device__ __inline__ CLGComplex _makeColorVector<CLGComplex>(BYTE spin, BYTE colorIdx)
{
    return _onec;
}

template<> __device__ __inline__ deviceSU2Vector _makeColorVector<deviceSU2Vector>(BYTE spin, BYTE colorIdx)
{
    if (colorIdx >= 2)
    {
        return _makeId<deviceSU2Vector>();
    }
    return deviceSU2Vector::makeOneSU2VectorColor(colorIdx);
}

template<> __device__ __inline__ deviceSU3Vector _makeColorVector<deviceSU3Vector>(BYTE spin, BYTE colorIdx)
{
    if (colorIdx >= 3)
    {
        return _makeId<deviceSU3Vector>();
    }
    return deviceSU3Vector::makeOneSU3VectorColor(colorIdx);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _makeColorVector<deviceWilsonVectorSU3>(BYTE spin, BYTE colorIdx)
{
    return deviceWilsonVectorSU3::makeOneWilsonVectorSU3SpinColor(spin, colorIdx);
}

#define _make_impl(n) \
template<> __device__ __inline__ deviceSU##n _makeId<deviceSU##n>() { return deviceSU##n::makeSUNId(); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeId<deviceSU##n##Vector>() { return deviceSU##n##Vector::makeOneSUNVector(); } \
template<> __device__ __inline__ deviceSU##n _makeZero<deviceSU##n>() { return deviceSU##n::makeSUNZero(); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeZero<deviceSU##n##Vector>() { return deviceSU##n##Vector::makeZeroSUNVector(); } \
template<> __device__ __inline__ deviceSU##n _makeContract<deviceSU##n, deviceSU##n##Vector>(const deviceSU##n##Vector& left, const deviceSU##n##Vector& right) { return deviceSU##n::makeSUNContractV(left, right); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeGaussian<deviceSU##n##Vector>(UINT fatIdx) { return deviceSU##n##Vector::makeRandomGaussian(fatIdx); } \
template<> __device__ __inline__ deviceSU##n _makeGaussian<deviceSU##n>(UINT fatIdx) { return deviceSU##n::makeSUNRandomGenerator(fatIdx); } \
template<> __device__ __inline__ deviceSU##n _makeSumGenerator(Real factor) { return deviceSU##n::makeSUNSumGenerator(factor); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeZ4<deviceSU##n##Vector>(UINT fatIdx) { return deviceSU##n##Vector::makeRandomZ4(fatIdx); } \
template<> __device__ __inline__ deviceSU##n _makeRandom<deviceSU##n>(UINT fatIdx) { return deviceSU##n::makeSUNRandom(fatIdx); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeRandom<deviceSU##n##Vector>(UINT fatIdx) { return deviceSU##n##Vector::makeRandom(fatIdx); } \
template<> __device__ __inline__ deviceSU##n##Vector _makeColorVector<deviceSU##n##Vector>(BYTE spinIdx, BYTE colorIdx) { if (colorIdx >= n) { return _makeId<deviceSU##n##Vector>(); } return deviceSU##n##Vector::makeOneSUNVectorColor(colorIdx); } \


#define _make_all_imp(n) _DEF_F_N(n, _make_impl)
_make_all_imp(_MAX_SUN)

#define _make_slnc_impl(n) \
template<> __device__ __inline__ deviceSL##n##C _makeId<deviceSL##n##C>() { return deviceSL##n##C::makeSLNCId(); } \
template<> __device__ __inline__ deviceSL##n##C _makeZero<deviceSL##n##C>() { return deviceSL##n##C::makeSLNCZero(); } \
template<> __device__ __inline__ deviceSL##n##C _makeContract<deviceSL##n##C, deviceSU##n##Vector>(const deviceSU##n##Vector& left, const deviceSU##n##Vector& right) { return deviceSL##n##C::makeSLNCContractV(left, right); } \
template<> __device__ __inline__ deviceSL##n##C _makeGaussian<deviceSL##n##C>(UINT fatIdx) { return deviceSL##n##C::makeSLNCRandomGenerator(fatIdx); } \
template<> __device__ __inline__ deviceSL##n##C _makeSumGenerator(Real factor) { return deviceSL##n##C::makeSLNCSumGenerator(factor); } \
template<> __device__ __inline__ deviceSL##n##C _makeRandom<deviceSL##n##C>(UINT fatIdx) { return deviceSL##n##C::makeSLNCRandom(fatIdx); }

#define _make_all_imp_slnc(n) _DEF_FSLNC_N(n, _make_slnc_impl)
_make_all_imp_slnc(_MAX_SLNC)

#define _make_un_impl(n) \
template<> __device__ __inline__ deviceU##n _makeId<deviceU##n>() { return deviceU##n::makeUNId(); } \
template<> __device__ __inline__ deviceU##n _makeZero<deviceU##n>() { return deviceU##n::makeUNZero(); } \
template<> __device__ __inline__ deviceU##n _makeContract<deviceU##n, deviceSU##n##Vector>(const deviceSU##n##Vector& left, const deviceSU##n##Vector& right) { return deviceU##n::makeUNContractV(left, right); } \
template<> __device__ __inline__ deviceU##n _makeGaussian<deviceU##n>(UINT fatIdx) { return deviceU##n::makeUNRandomGenerator(fatIdx); } \
template<> __device__ __inline__ deviceU##n _makeSumGenerator(Real factor) { return deviceU##n::makeUNSumGenerator(factor); } \
template<> __device__ __inline__ deviceU##n _makeRandom<deviceU##n>(UINT fatIdx) { return deviceU##n::makeUNRandom(fatIdx); }

#define _make_all_imp_un(n) _DEF_FUN_N(n, _make_un_impl)
_make_all_imp_un(_MAX_UN)

#define _make_on_impl(n) \
template<> __device__ __inline__ deviceO##n _makeId<deviceO##n>() { return deviceO##n::makeONId(); } \
template<> __device__ __inline__ deviceO##n _makeZero<deviceO##n>() { return deviceO##n::makeONZero(); } \
template<> __device__ __inline__ deviceO##n _makeContract<deviceO##n, deviceSU##n##Vector>(const deviceSU##n##Vector& left, const deviceSU##n##Vector& right) { return deviceO##n::makeONContractV(left, right); } \
template<> __device__ __inline__ deviceO##n _makeGaussian<deviceO##n>(UINT fatIdx) { return deviceO##n::makeONRandomGenerator(fatIdx); } \
template<> __device__ __inline__ deviceO##n _makeSumGenerator(Real factor) { return deviceO##n::makeONSumGenerator(factor); } \
template<> __device__ __inline__ deviceO##n _makeRandom<deviceO##n>(UINT fatIdx) { return deviceO##n::makeONRandom(fatIdx); }

#define _make_all_imp_on(n) _DEF_FON_N(n, _make_on_impl)
_make_all_imp_on(_MAX_ON)

#define _make_son_impl(n) \
template<> __device__ __inline__ deviceSO##n _makeId<deviceSO##n>() { return deviceSO##n::makeSONId(); } \
template<> __device__ __inline__ deviceSO##n _makeZero<deviceSO##n>() { return deviceSO##n::makeSONZero(); } \
template<> __device__ __inline__ deviceSO##n _makeContract<deviceSO##n, deviceSU##n##Vector>(const deviceSU##n##Vector& left, const deviceSU##n##Vector& right) { return deviceSO##n::makeSONContractV(left, right); } \
template<> __device__ __inline__ deviceSO##n _makeGaussian<deviceSO##n>(UINT fatIdx) { return deviceSO##n::makeSONRandomGenerator(fatIdx); } \
template<> __device__ __inline__ deviceSO##n _makeSumGenerator(Real factor) { return deviceSO##n::makeSONSumGenerator(factor); } \
template<> __device__ __inline__ deviceSO##n _makeRandom<deviceSO##n>(UINT fatIdx) { return deviceSO##n::makeSONRandom(fatIdx); }

#define _make_all_imp_son(n) _DEF_FSON_N(n, _make_son_impl)
_make_all_imp_son(_MAX_SON)

template<typename T> __device__ __inline__ void _dagger(T& element) = delete;
template<typename T> __device__ __inline__ T _daggerC(const T& element) = delete;

template<> __device__ __inline__ void _dagger<Real>(Real& element)
{
    //do nothing
}

template<> __device__ __inline__ void _dagger<CLGComplex>(CLGComplex& element)
{
    element.y = -element.y;
}

template<> __device__ __inline__ void _dagger<deviceSU2Vector>(deviceSU2Vector& element)
{
    element.Conjugate();
}

template<> __device__ __inline__ void _dagger<deviceSU3Vector>(deviceSU3Vector& element)
{
    element.Conjugate();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceSUNVector<N, NoE>& element)
{
    element.Conjugate();
}

template<> __device__ __inline__ void _dagger<deviceSU2>(deviceSU2& element)
{
    element.Dagger();
}

template<> __device__ __inline__ void _dagger<deviceSU3>(deviceSU3& element)
{
    element.Dagger();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceSUN<N, NoE>& element)
{
    element.Dagger();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceSLNC<N, NoE>& element)
{
    element.Dagger();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceUN<N, NoE>& element)
{
    element.Dagger();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceON<N, NoE>& element)
{
    element.Transpose();
}

template<INT N, INT NoE> __device__ __inline__ void _dagger(deviceSON<N, NoE>& element)
{
    element.Transpose();
}

template<INT N> __device__ __inline__ void _dagger(deviceZN<N>& element)
{
    element.Dagger();
}

template<> __device__ __inline__ void _dagger<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& element)
{
    element.Conjugate();
}

template<> __device__ __inline__ Real _daggerC<Real>(const Real& element)
{
    return element;
}

template<> __device__ __inline__ CLGComplex _daggerC<CLGComplex>(const CLGComplex& element)
{
    return _cuConjf(element);
}

template<> __device__ __inline__ deviceSU2Vector _daggerC<deviceSU2Vector>(const deviceSU2Vector& element)
{
    return element.ConjugateC();
}

template<> __device__ __inline__ deviceSU3Vector _daggerC<deviceSU3Vector>(const deviceSU3Vector& element)
{
    return element.ConjugateC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSUNVector<N, NoE> _daggerC(const deviceSUNVector<N, NoE>& element)
{
    return element.ConjugateC();
}

template<> __device__ __inline__ deviceSU2 _daggerC<deviceSU2>(const deviceSU2& element)
{
    return element.DaggerC();
}

template<> __device__ __inline__ deviceSU3 _daggerC<deviceSU3>(const deviceSU3& element)
{
    return element.DaggerC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSUN<N, NoE> _daggerC(const deviceSUN<N, NoE>& element)
{
    return element.DaggerC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSLNC<N, NoE> _daggerC(const deviceSLNC<N, NoE>& element)
{
    return element.DaggerC();
}

template<INT N, INT NoE> __device__ __inline__ deviceUN<N, NoE> _daggerC(const deviceUN<N, NoE>& element)
{
    return element.DaggerC();
}

template<INT N, INT NoE> __device__ __inline__ deviceON<N, NoE> _daggerC(const deviceON<N, NoE>& element)
{
    return element.TransposeC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSON<N, NoE> _daggerC(const deviceSON<N, NoE>& element)
{
    return element.TransposeC();
}

template<INT N> __device__ __inline__ deviceZN<N> _daggerC(const deviceZN<N>& element)
{
    return element.DaggerC();
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _daggerC<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& element)
{
    return element.ConjugateC();
}

template<typename T> __device__ __inline__ void _oppo(T& element) = delete;
template<typename T> __device__ __inline__ T _oppoC(const T& element) = delete;

template<> __device__ __inline__ void _oppo<Real>(Real& element)
{
    element = -element;
}

template<> __device__ __inline__ void _oppo<CLGComplex>(CLGComplex& element)
{
    element.x = -element.x;
    element.y = -element.y;
}

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ void _oppo<FLOAT>(FLOAT& element)
{
    element = -element;
}
template<> __device__ __inline__ void _oppo<cuComplex>(cuComplex& element)
{
    element.x = -element.x;
    element.y = -element.y;
}
#else
template<> __device__ __inline__ void _oppo<DOUBLE>(DOUBLE& element)
{
    element = -element;
}
template<> __device__ __inline__ void _oppo<cuDoubleComplex>(cuDoubleComplex& element)
{
    element.x = -element.x;
    element.y = -element.y;
}
#endif

template<> __device__ __inline__ void _oppo<deviceSU2Vector>(deviceSU2Vector& element)
{
    element.Opposite();
}

template<> __device__ __inline__ void _oppo<deviceSU3Vector>(deviceSU3Vector& element)
{
    element.Opposite();
}

template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceSUNVector<N, NoE>& element)
{
    element.Opposite();
}

template<> __device__ __inline__ void _oppo<deviceSU2>(deviceSU2& element)
{
    element.Opposite();
}

template<> __device__ __inline__ void _oppo<deviceSU3>(deviceSU3& element)
{
    element.Opposite();
}

template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceSUN<N, NoE>& element)
{
    element.Opposite();
}

template<INT N> __device__ __inline__ void _oppo(deviceZN<N>& element)
{
    element.Opposite();
}

template<> __device__ __inline__ void _oppo<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& element)
{
    element.Opposite();
}

template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceSLNC<N, NoE>& element) { element.MulReal(F(-1.0)); }
template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceUN<N, NoE>& element) { element.MulReal(F(-1.0)); }
template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceON<N, NoE>& element) { element.MulReal(F(-1.0)); }
template<INT N, INT NoE> __device__ __inline__ void _oppo(deviceSON<N, NoE>& element) { element.MulReal(F(-1.0)); }

template<> __device__ __inline__ Real _oppoC<Real>(const Real& element)
{
    return -element;
}

template<> __device__ __inline__ CLGComplex _oppoC<CLGComplex>(const CLGComplex& element)
{
    return _make_cuComplex(-element.x, -element.y);
}

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ FLOAT _oppoC<FLOAT>(const FLOAT& element)
{
    return -element;
}
template<> __device__ __inline__ cuComplex _oppoC<cuComplex>(const cuComplex& element)
{
    return make_cuComplex(-element.x, -element.y);
}
#else
template<> __device__ __inline__ DOUBLE _oppoC<DOUBLE>(const DOUBLE& element)
{
    return -element;
}
template<> __device__ __inline__ cuDoubleComplex _oppoC<cuDoubleComplex>(const cuDoubleComplex& element)
{
    return make_cuDoubleComplex(-element.x, -element.y);
}
#endif

template<> __device__ __inline__ deviceSU2Vector _oppoC<deviceSU2Vector>(const deviceSU2Vector& element)
{
    return element.OppositeC();
}

template<> __device__ __inline__ deviceSU3Vector _oppoC<deviceSU3Vector>(const deviceSU3Vector& element)
{
    return element.OppositeC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSUNVector<N, NoE> _oppoC(const deviceSUNVector<N, NoE>& element)
{
    return element.OppositeC();
}

template<> __device__ __inline__ deviceSU2 _oppoC<deviceSU2>(const deviceSU2& element)
{
    return element.OppositeC();
}

template<> __device__ __inline__ deviceSU3 _oppoC<deviceSU3>(const deviceSU3& element)
{
    return element.OppositeC();
}

template<INT N, INT NoE> __device__ __inline__ deviceSUN<N, NoE> _oppoC(const deviceSUN<N, NoE>& element)
{
    return element.OppositeC();
}

template<INT N> __device__ __inline__ deviceZN<N> _oppoC(const deviceZN<N>& element)
{
    return element.OppositeC();
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _oppoC<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& element)
{
    return element.OppositeC();
}

template<typename T> __device__ __inline__ T _rcpC(const T& element) = delete;

template<> __device__ __inline__ Real _rcpC<Real>(const Real& element)
{
    return F(1.0) / element;
}

template<> __device__ __inline__ CLGComplex _rcpC<CLGComplex>(const CLGComplex& element)
{
    return _cuCdivf(_onec, element);
}

template<> __device__ __inline__ deviceSU2 _rcpC<deviceSU2>(const deviceSU2& element)
{
    return element.Inverse();
}

template<> __device__ __inline__ deviceSU3 _rcpC<deviceSU3>(const deviceSU3& element)
{
    return element.Inverse();
}

template<INT N, INT NoE> __device__ __inline__ deviceSUN<N, NoE> _rcpC(const deviceSUN<N, NoE>& element)
{
    return element.Inverse();
}

template<INT N, INT NoE> __device__ __inline__ deviceSLNC<N, NoE> _rcpC(const deviceSLNC<N, NoE>& element)
{
    return element.InverseC();
}


template<typename TLeft, typename TRight> __device__ __inline__ TLeft _addC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _add(TLeft& left, const TRight& right)
{
    left = _addC(left, right);
}


template<> __device__ __inline__ Real _addC<Real, Real>(const Real& left, const Real& right)
{
    return left + right;
}
template<> __device__ __inline__ CLGComplex _addC<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __device__ __inline__ void _add<CLGComplex, Real>(CLGComplex& left, const Real& right)
{
    left.x = left.x + right;
}
template<> __device__ __inline__ CLGComplex _addC<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCaddf(left, right);
}
template<> __device__ __inline__ void _add<CLGComplex, CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left.x = left.x + right.x;
    left.y = left.y + right.y;
}

template<> __device__ __inline__ cuComplex _addC<cuComplex, cuDoubleComplex>(const cuComplex& left, const cuDoubleComplex& right)
{
    return make_cuComplex(left.x + static_cast<FLOAT>(right.x), left.y + static_cast<FLOAT>(right.y));
}
template<> __device__ __inline__ cuDoubleComplex _addC<cuDoubleComplex, cuComplex>(const cuDoubleComplex& left, const cuComplex& right)
{
    return make_cuDoubleComplex(left.x + static_cast<DOUBLE>(right.x), left.y + static_cast<DOUBLE>(right.y));
}

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ FLOAT _addC<FLOAT, Real>(const FLOAT& left, const Real& right)
{
    return left + static_cast<FLOAT>(right);
}
template<> __device__ __inline__ FLOAT _addC<FLOAT, FLOAT>(const FLOAT& left, const FLOAT& right)
{
    return left + right;
}
template<> __device__ __inline__ CLGComplex _addC<CLGComplex, FLOAT>(const CLGComplex& left, const FLOAT& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + static_cast<DOUBLE>(right);
    return ret;
}
template<> __device__ __inline__ void _add<CLGComplex, FLOAT>(CLGComplex& left, const FLOAT& right)
{
    left.x = left.x + static_cast<DOUBLE>(right);
}
template<> __device__ __inline__ cuComplex _addC<cuComplex, FLOAT>(const cuComplex& left, const FLOAT& right)
{
    cuComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __device__ __inline__ void _add<cuComplex, FLOAT>(cuComplex& left, const FLOAT& right)
{
    left.x = left.x + right;
}
template<> __device__ __inline__ cuComplex _addC<cuComplex, cuComplex>(const cuComplex& left, const cuComplex& right)
{
    return make_cuComplex(left.x + right.x, left.y + right.y);
}
#else
template<> __device__ __inline__ DOUBLE _addC<DOUBLE, Real>(const DOUBLE& left, const Real& right)
{
    return left + static_cast<DOUBLE>(right);
}
template<> __device__ __inline__ DOUBLE _addC<DOUBLE, DOUBLE>(const DOUBLE& left, const DOUBLE& right)
{
    return left + right;
}
template<> __device__ __inline__ CLGComplex _addC<CLGComplex, DOUBLE>(const CLGComplex& left, const DOUBLE& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + static_cast<Real>(right);
    return ret;
}
template<> __device__ __inline__ void _add<CLGComplex, DOUBLE>(CLGComplex& left, const DOUBLE& right)
{
    left.x = left.x + static_cast<Real>(right);
}
template<> __device__ __inline__ cuDoubleComplex _addC<cuDoubleComplex, DOUBLE>(const cuDoubleComplex& left, const DOUBLE& right)
{
    cuDoubleComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __device__ __inline__ void _add<cuDoubleComplex, DOUBLE>(cuDoubleComplex& left, const DOUBLE& right)
{
    left.x = left.x + right;
}
template<> __device__ __inline__ cuDoubleComplex _addC<cuDoubleComplex, cuDoubleComplex>(const cuDoubleComplex& left, const cuDoubleComplex& right)
{
    return make_cuDoubleComplex(left.x + right.x, left.y + right.y);
}
#endif

#if _CLG_DOUBLEFLOAT
#define __DEFINE_TWO_ELEMENT_Func(TYPENAME, FUNC1, FUNC2) \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, Real>(const TYPENAME& left, const Real& right) \
{ \
    return left.FUNC2##RealC(right); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, FLOAT>(const TYPENAME& left, const FLOAT& right) \
{ \
    return left.FUNC2##RealC(static_cast<Real>(right)); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, CLGComplex>(const TYPENAME& left, const CLGComplex& right) \
{ \
    return left.FUNC2##CompC(right); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, TYPENAME>(const TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2##C(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, Real>(TYPENAME& left, const Real& right) \
{ \
    return left.FUNC2##Real(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, FLOAT>(TYPENAME& left, const FLOAT& right) \
{ \
    return left.FUNC2##Real(static_cast<Real>(right)); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, CLGComplex>(TYPENAME& left, const CLGComplex& right) \
{ \
    return left.FUNC2##Comp(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, TYPENAME>(TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2(right); \
} 
#else
#define __DEFINE_TWO_ELEMENT_Func(TYPENAME, FUNC1, FUNC2) \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, Real>(const TYPENAME& left, const Real& right) \
{ \
    return left.FUNC2##RealC(right); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, DOUBLE>(const TYPENAME& left, const DOUBLE& right) \
{ \
    return left.FUNC2##RealC(static_cast<Real>(right)); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, CLGComplex>(const TYPENAME& left, const CLGComplex& right) \
{ \
    return left.FUNC2##CompC(right); \
} \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, TYPENAME>(const TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2##C(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, Real>(TYPENAME& left, const Real& right) \
{ \
    return left.FUNC2##Real(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, DOUBLE>(TYPENAME& left, const DOUBLE& right) \
{ \
    return left.FUNC2##Real(static_cast<Real>(right)); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, CLGComplex>(TYPENAME& left, const CLGComplex& right) \
{ \
    return left.FUNC2##Comp(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, TYPENAME>(TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2(right); \
} 
#endif

__DEFINE_TWO_ELEMENT_Func(deviceSU2Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU3Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _add, Add)

template<INT N> __device__ __inline__ void _add(deviceZN<N>& left, const deviceZN<N>& right) { left.Add(right); }
template<INT N> __device__ __inline__ void _add(deviceZN<N>& left, const Real& right) { left.AddReal(right); }
template<INT N> __device__ __inline__ deviceZN<N> _addC(const deviceZN<N>& left, const deviceZN<N>& right) { return left.AddC(right); }
template<INT N> __device__ __inline__ deviceZN<N> _addC(const deviceZN<N>& left, const Real& right) { return left.AddRealC(right); }

template<> __device__ __inline__ void _add<deviceWilsonVectorSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
{
    left.Add(right);
}

template<typename TLeft, typename TRight> __host__ __inline__ TLeft _addCHost(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __host__ __inline__ void _addHost(TLeft& left, const TRight& right)
{
    left = _addCHost(left, right);
}


template<> __host__ __inline__ Real _addCHost<Real, Real>(const Real& left, const Real& right)
{
    return left + right;
}
template<> __host__ __inline__ CLGComplex _addCHost<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __host__ __inline__ void _addHost<CLGComplex, Real>(CLGComplex& left, const Real& right)
{
    left.x = left.x + right;
}
template<> __host__ __inline__ CLGComplex _addCHost<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCaddf(left, right);
}
template<> __host__ __inline__ void _addHost<CLGComplex, CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left.x = left.x + right.x;
    left.y = left.y + right.y;
}

template<> __host__ __inline__ cuComplex _addCHost<cuComplex, cuDoubleComplex>(const cuComplex& left, const cuDoubleComplex& right)
{
    return make_cuComplex(left.x + static_cast<FLOAT>(right.x), left.y + static_cast<FLOAT>(right.y));
}
template<> __host__ __inline__ cuDoubleComplex _addCHost<cuDoubleComplex, cuComplex>(const cuDoubleComplex& left, const cuComplex& right)
{
    return make_cuDoubleComplex(left.x + static_cast<DOUBLE>(right.x), left.y + static_cast<DOUBLE>(right.y));
}

#if _CLG_DOUBLEFLOAT
template<> __host__ __inline__ FLOAT _addCHost<FLOAT, Real>(const FLOAT& left, const Real& right)
{
    return left + static_cast<FLOAT>(right);
}
template<> __host__ __inline__ FLOAT _addCHost<FLOAT, FLOAT>(const FLOAT& left, const FLOAT& right)
{
    return left + right;
}
template<> __host__ __inline__ CLGComplex _addCHost<CLGComplex, FLOAT>(const CLGComplex& left, const FLOAT& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + static_cast<DOUBLE>(right);
    return ret;
}
template<> __host__ __inline__ void _addHost<CLGComplex, FLOAT>(CLGComplex& left, const FLOAT& right)
{
    left.x = left.x + static_cast<DOUBLE>(right);
}
template<> __host__ __inline__ cuComplex _addCHost<cuComplex, FLOAT>(const cuComplex& left, const FLOAT& right)
{
    cuComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __host__ __inline__ void _addHost<cuComplex, FLOAT>(cuComplex& left, const FLOAT& right)
{
    left.x = left.x + right;
}
template<> __host__ __inline__ cuComplex _addCHost<cuComplex, cuComplex>(const cuComplex& left, const cuComplex& right)
{
    return make_cuComplex(left.x + right.x, left.y + right.y);
}
#else
template<> __host__ __inline__ DOUBLE _addCHost<DOUBLE, Real>(const DOUBLE& left, const Real& right)
{
    return left + static_cast<DOUBLE>(right);
}
template<> __host__ __inline__ DOUBLE _addCHost<DOUBLE, DOUBLE>(const DOUBLE& left, const DOUBLE& right)
{
    return left + right;
}
template<> __host__ __inline__ CLGComplex _addCHost<CLGComplex, DOUBLE>(const CLGComplex& left, const DOUBLE& right)
{
    CLGComplex ret = left;
    ret.x = ret.x + static_cast<Real>(right);
    return ret;
}
template<> __host__ __inline__ void _addHost<CLGComplex, DOUBLE>(CLGComplex& left, const DOUBLE& right)
{
    left.x = left.x + static_cast<Real>(right);
}
template<> __host__ __inline__ cuDoubleComplex _addCHost<cuDoubleComplex, DOUBLE>(const cuDoubleComplex& left, const DOUBLE& right)
{
    cuDoubleComplex ret = left;
    ret.x = ret.x + right;
    return ret;
}
template<> __host__ __inline__ void _addHost<cuDoubleComplex, DOUBLE>(cuDoubleComplex& left, const DOUBLE& right)
{
    left.x = left.x + right;
}
template<> __host__ __inline__ cuDoubleComplex _addCHost<cuDoubleComplex, cuDoubleComplex>(const cuDoubleComplex& left, const cuDoubleComplex& right)
{
    return make_cuDoubleComplex(left.x + right.x, left.y + right.y);
}
#endif


template<typename TLeft, typename TRight> __device__ __inline__ TLeft _subC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _sub(TLeft& left, const TRight& right) = delete;


template<> __device__ __inline__ Real _subC<Real, Real>(const Real& left, const Real& right)
{
    return left - right;
}
template<> __device__ __inline__ void _sub<Real, Real>(Real& left, const Real& right)
{
    left = left - right;
}
template<> __device__ __inline__ CLGComplex _subC<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    CLGComplex ret = left;
    ret.x = ret.x - right;
    return ret;
}
template<> __device__ __inline__ void _sub<CLGComplex, Real>(CLGComplex& left, const Real& right)
{
    left.x = left.x - right;
}
template<> __device__ __inline__ CLGComplex _subC<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCsubf(left, right);
}
template<> __device__ __inline__ void _sub<CLGComplex, CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left.x = left.x - right.x;
    left.y = left.y - right.y;
}

__DEFINE_TWO_ELEMENT_Func(deviceSU2Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU3Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _sub, Sub)

template<INT N> __device__ __inline__ void _sub(deviceZN<N>& left, const deviceZN<N>& right) { left.Sub(right); }
template<INT N> __device__ __inline__ void _sub(deviceZN<N>& left, const Real& right) { left.SubReal(right); }
template<INT N> __device__ __inline__ deviceZN<N> _subC(const deviceZN<N>& left, const deviceZN<N>& right) { return left.SubC(right); }
template<INT N> __device__ __inline__ deviceZN<N> _subC(const deviceZN<N>& left, const Real& right) { return left.SubRealC(right); }

template<> __device__ __inline__ void _sub<deviceWilsonVectorSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
{
    left.Sub(right);
}

template<typename TLeft, typename TRight> __device__ __inline__ TLeft _divC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _div(TLeft& left, const TRight& right)
{
    left = _divC(left, right);
}
template<typename TLeft, typename TRight> __host__ __inline__ TLeft _divCHost(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __host__ __inline__ void _divHost(TLeft& left, const TRight& right)
{
    left = _divCHost(left, right);
}

#define __DEFINE_TWO_ELEMENT_Func_Div_Slash(t1, t2) \
template<> __device__ __inline__ t1 _divC<t1, t2>(const t1& left, const t2& right) \
{ \
    return static_cast<t1>(left / right); \
} \
template<> __host__ __inline__ t1 _divCHost<t1, t2>(const t1& left, const t2& right) \
{ \
    return static_cast<t1>(left / right); \
}

__DEFINE_TWO_ELEMENT_Func_Div_Slash(Real, Real)
#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_Slash(Real, FLOAT)
#else
__DEFINE_TWO_ELEMENT_Func_Div_Slash(Real, DOUBLE)
#endif
__DEFINE_TWO_ELEMENT_Func_Div_Slash(Real, INT)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(Real, UINT)
#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_Slash(FLOAT, Real)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(FLOAT, FLOAT)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(FLOAT, INT)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(FLOAT, UINT)
#else
__DEFINE_TWO_ELEMENT_Func_Div_Slash(DOUBLE, Real)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(DOUBLE, DOUBLE)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(DOUBLE, INT)
__DEFINE_TWO_ELEMENT_Func_Div_Slash(DOUBLE, UINT)
#endif

__DEFINE_TWO_ELEMENT_Func_Div_Slash(INT, Real)
#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_Slash(INT, FLOAT)
#else
__DEFINE_TWO_ELEMENT_Func_Div_Slash(INT, DOUBLE)
#endif
__DEFINE_TWO_ELEMENT_Func_Div_Slash(INT, INT)

#define __DEFINE_TWO_ELEMENT_Func_Div_FC(t1, t2, method, castint1, castint2, castout) \
template<> __device__ __inline__ t1 _divC<t1, t2>(const t1& left, const t2& right) \
{ \
    return castout(method(castint1(left), castint2(right))); \
} \
template<> __host__ __inline__ t1 _divCHost<t1, t2>(const t1& left, const t2& right) \
{ \
    return castout(method(castint1(left), castint2(right))); \
}

#define __DEFINE_TWO_ELEMENT_Func_Div_FC_DH(t1, t2, method, methodhost, castint1, castint2, castout) \
template<> __device__ __inline__ t1 _divC<t1, t2>(const t1& left, const t2& right) \
{ \
    return castout(method(castint1(left), castint2(right))); \
} \
template<> __host__ __inline__ t1 _divCHost<t1, t2>(const t1& left, const t2& right) \
{ \
    return castout(methodhost(castint1(left), castint2(right))); \
}

__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(CLGComplex, Real, cuCdivf_cr, cuCdivf_cr_host, , , )
__DEFINE_TWO_ELEMENT_Func_Div_FC(CLGComplex, CLGComplex, _cuCdivf, , , )
#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(CLGComplex, FLOAT, cuCdivf_cr, cuCdivf_cr_host, , static_cast<Real>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC(CLGComplex, cuComplex, _cuCdivf, , _cToRealC, )
#else
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(CLGComplex, DOUBLE, cuCdivf_cr, cuCdivf_cr_host, , static_cast<Real>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC(CLGComplex, cuDoubleComplex, _cuCdivf, , _cToRealC, )
#endif
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(CLGComplex, INT, cuCdivf_cr, cuCdivf_cr_host, , static_cast<Real>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(CLGComplex, UINT, cuCdivf_cr, cuCdivf_cr_host, , static_cast<Real>, )

#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuComplex, Real, cuCdivf_cr, cuCdivf_cr_host, _cToRealC, , _cToFloat)
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuComplex, FLOAT, cuCdivf_cr, cuCdivf_cr_host, _cToRealC, static_cast<Real>, _cToFloat)
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuComplex, INT, cuCdivf_cr, cuCdivf_cr_host, _cToRealC, static_cast<Real>, _cToFloat)
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuComplex, UINT, cuCdivf_cr, cuCdivf_cr_host, _cToRealC, static_cast<Real>, _cToFloat)
__DEFINE_TWO_ELEMENT_Func_Div_FC(cuComplex, CLGComplex, cuCdivf, , _cToFloat, )
__DEFINE_TWO_ELEMENT_Func_Div_FC(cuComplex, cuComplex, cuCdivf, , , )
#else
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuDoubleComplex, Real, cuCdivf_cd, cuCdivf_cd_host, , static_cast<DOUBLE>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuDoubleComplex, INT, cuCdivf_cd, cuCdivf_cd_host, , static_cast<DOUBLE>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuDoubleComplex, UINT, cuCdivf_cd, cuCdivf_cd_host, , static_cast<DOUBLE>, )
__DEFINE_TWO_ELEMENT_Func_Div_FC_DH(cuDoubleComplex, DOUBLE, cuCdivf_cd, cuCdivf_cd_host, , , )
__DEFINE_TWO_ELEMENT_Func_Div_FC(cuDoubleComplex, CLGComplex, cuCdiv, , _cToDouble, )
__DEFINE_TWO_ELEMENT_Func_Div_FC(cuDoubleComplex, cuDoubleComplex, cuCdiv, , , )
#endif



#define __DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(t1, t2, methodname, castname) \
template<INT e, INT me> __device__ __inline__ t1<e, me> _divC(const t1<e, me>& left, const t2& right) \
{ \
    return left.Div##methodname##C(castname(right)); \
} \
template<INT e, INT me> __device__ __inline__ void _div(t1<e, me>& left, const t2& right) \
{ \
    left.Div##methodname(castname(right)); \
}

__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, Real, Real, )
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, INT, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, UINT, Real, static_cast<Real>)

__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, Real, Real, )
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, INT, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, UINT, Real, static_cast<Real>)

__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, CLGComplex, Comp, _cToRealC)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, CLGComplex, Comp, _cToRealC)

#if _CLG_DOUBLEFLOAT
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, FLOAT, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, FLOAT, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, cuComplex, Comp, _cToRealC)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, cuComplex, Comp, _cToRealC)
#else
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, DOUBLE, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, DOUBLE, Real, static_cast<Real>)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUN, cuDoubleComplex, Comp, _cToRealC)
__DEFINE_TWO_ELEMENT_Func_Div_t_on_fc(deviceSUNVector, cuDoubleComplex, Comp, _cToRealC)
#endif

#if defined(__cplusplus) && defined(__CUDACC__)
template<typename T> __device__ __inline__ void _atomicAdd(T* left, const T& right) = delete;

#define __DEFINE_SINLE_ATOMIC_ADD(t) \
template<> __device__ __inline__ void _atomicAdd<t>(t* left, const t& right) { atomicAdd(left, right); } \

#define __DEFINE_SINLE_ATOMIC_ADDC(t) \
template<> __device__ __inline__ void _atomicAdd<t>(t* left, const t& right) { atomicAdd(&(left->x), right.x); atomicAdd(&(left->y), right.y); } \

__DEFINE_SINLE_ATOMIC_ADD(Real)
#if _CLG_DOUBLEFLOAT
__DEFINE_SINLE_ATOMIC_ADD(FLOAT)
#else
__DEFINE_SINLE_ATOMIC_ADD(DOUBLE)
#endif
__DEFINE_SINLE_ATOMIC_ADD(INT)
__DEFINE_SINLE_ATOMIC_ADD(UINT)

#if _CLG_DOUBLEFLOAT
__DEFINE_SINLE_ATOMIC_ADDC(cuComplex)
__DEFINE_SINLE_ATOMIC_ADDC(CLGComplex)
#else
__DEFINE_SINLE_ATOMIC_ADDC(cuDoubleComplex)
__DEFINE_SINLE_ATOMIC_ADDC(CLGComplex)
#endif

template<INT e, INT me> __device__ __inline__ void _atomicAdd(deviceSUN<e, me>* left, const deviceSUN<e, me>& right) 
{ 
    for (INT i = 0; i < e; ++i)
    {
        _atomicAdd(left->m_me + i, right.m_me[i]);
    }
}

template<INT e, INT me> __device__ __inline__ void _atomicAdd(deviceSUNVector<e, me>* left, const deviceSUNVector<e, me>& right)
{
    for (INT i = 0; i < e; ++i)
    {
        _atomicAdd(left->m_ve + i, right.m_ve[i]);
    }
}
#endif

template<typename TLeft, typename TRight> __device__ __inline__ TLeft _mulC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _mul(TLeft& left, const TRight& right)
{
    left = _mulC(left, right);
}

template<> __device__ __inline__ FLOAT _mulC<FLOAT, FLOAT>(const FLOAT& left, const FLOAT& right)
{
    return left * right;
}
template<> __device__ __inline__ DOUBLE _mulC<DOUBLE, DOUBLE>(const DOUBLE& left, const DOUBLE& right)
{
    return left * right;
}
template<> __device__ __inline__ FLOAT _mulC<FLOAT, DOUBLE>(const FLOAT& left, const DOUBLE& right)
{
    return left * static_cast<FLOAT>(right);
}
template<> __device__ __inline__ DOUBLE _mulC<DOUBLE, FLOAT>(const DOUBLE& left, const FLOAT& right)
{
    return left * static_cast<DOUBLE>(right);
}

template<> __device__ __inline__ Real _mulC<Real, CLGComplex>(const Real& left, const CLGComplex& right)
{
    return left * _cuCabsf(right);
}
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    return cuCmulf_cr(left, right);
}


#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, FLOAT>(const CLGComplex& left, const FLOAT& right)
{
    return cuCmulf_cr(left, static_cast<Real>(right));
}
template<> __device__ __inline__ cuComplex _mulC<cuComplex, FLOAT>(const cuComplex& left, const FLOAT& right)
{
    return _cToFloat(cuCmulf_cr(_cToRealC(left), static_cast<Real>(right)));
}
template<> __device__ __inline__ cuComplex _mulC<cuComplex, DOUBLE>(const cuComplex& left, const DOUBLE& right)
{
    return _cToFloat(cuCmulf_cr(_cToRealC(left), static_cast<Real>(right)));
}
#else
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, DOUBLE>(const CLGComplex& left, const DOUBLE& right)
{
    return cuCmulf_cr(left, static_cast<Real>(right));
}
template<> __device__ __inline__ cuDoubleComplex _mulC<cuDoubleComplex, FLOAT>(const cuDoubleComplex& left, const FLOAT& right)
{
    return cuCmulf_cd(left, static_cast<DOUBLE>(right));
}
template<> __device__ __inline__ cuDoubleComplex _mulC<cuDoubleComplex, DOUBLE>(const cuDoubleComplex& left, const DOUBLE& right)
{
    return cuCmulf_cd(left, right);
}
#endif

template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(left, right);
}
#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, cuComplex>(const CLGComplex& left, const cuComplex& right)
{
    return _cuCmulf(left, _cToRealC(right));
}
#else
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, cuDoubleComplex>(const CLGComplex& left, const cuDoubleComplex& right)
{
    return _cuCmulf(left, _cToRealC(right));
}
#endif

#if _CLG_DOUBLEFLOAT
template<> __device__ __inline__ cuComplex _mulC<cuComplex, CLGComplex>(const cuComplex& left, const CLGComplex& right)
{
    return cuCmulf(left, _cToFloat(right));
}
template<> __device__ __inline__ cuComplex _mulC<cuComplex, cuComplex>(const cuComplex& left, const cuComplex& right)
{
    return cuCmulf(left, right);
}
//template<> __device__ __inline__ cuComplex _mulC<cuComplex, cuDoubleComplex>(const cuComplex& left, const cuDoubleComplex& right)
//{
//    return cuCmulf(left, _cToFloat(right));
//}
#else
template<> __device__ __inline__ cuDoubleComplex _mulC<cuDoubleComplex, CLGComplex>(const cuDoubleComplex& left, const CLGComplex& right)
{
    return cuCmul(left, _cToDouble(right));
}
//template<> __device__ __inline__ cuDoubleComplex _mulC<cuDoubleComplex, cuComplex>(const cuDoubleComplex& left, const cuComplex& right)
//{
//    return cuCmul(left, _cToDouble(right));
//}
template<> __device__ __inline__ cuDoubleComplex _mulC<cuDoubleComplex, cuDoubleComplex>(const cuDoubleComplex& left, const cuDoubleComplex& right)
{
    return cuCmul(left, right);
}
#endif

__DEFINE_TWO_ELEMENT_Func(deviceSU2Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU3Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _mul, Mul)

template<INT N> __device__ __inline__ void _mul(deviceZN<N>& left, const deviceZN<N>& right) { left.Mul(right); }
template<INT N> __device__ __inline__ void _mul(deviceZN<N>& left, const Real& right) { left.MulReal(right); }
template<INT N> __device__ __inline__ void _mul(deviceZN<N>& left, const CLGComplex& right) { left.MulComp(right); }
template<INT N> __device__ __inline__ deviceZN<N> _mulC(const deviceZN<N>& left, const deviceZN<N>& right) { return left.MulC(right); }
template<INT N> __device__ __inline__ deviceZN<N> _mulC(const deviceZN<N>& left, const Real& right) { return left.MulRealC(right); }
template<INT N> __device__ __inline__ deviceZN<N> _mulC(const deviceZN<N>& left, const CLGComplex& right) { return left.MulCompC(right); }

template<> __device__ __inline__ void _mul<deviceWilsonVectorSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
{
    left.Mul(right);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _mulC<deviceWilsonVectorSU3, deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
{
    return left.MulC(right);
}

template<> __device__ __inline__ void _mul<deviceWilsonVectorSU3, CLGComplex>(deviceWilsonVectorSU3& left, const CLGComplex& right)
{
    left.MulComp(right);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _mulC<deviceWilsonVectorSU3, CLGComplex>(const deviceWilsonVectorSU3& left, const CLGComplex& right)
{
    return left.MulCompC(right);
}

template<> __device__ __inline__ void _mul<deviceWilsonVectorSU3, Real>(deviceWilsonVectorSU3& left, const Real& right)
{
    left.MulReal(right);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _mulC<deviceWilsonVectorSU3, Real>(const deviceWilsonVectorSU3& left, const Real& right)
{
    return left.MulRealC(right);
}

template<typename TLeft, typename TRight> __host__ __inline__ TLeft _mulCHost(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __host__ __inline__ void _mulHost(TLeft& left, const TRight& right)
{
    left = _mulCHost(left, right);
}

template<> __host__ __inline__ FLOAT _mulCHost<FLOAT, FLOAT>(const FLOAT& left, const FLOAT& right)
{
    return left * right;
}
template<> __host__ __inline__ DOUBLE _mulCHost<DOUBLE, DOUBLE>(const DOUBLE& left, const DOUBLE& right)
{
    return left * right;
}
template<> __host__ __inline__ FLOAT _mulCHost<FLOAT, DOUBLE>(const FLOAT& left, const DOUBLE& right)
{
    return left * static_cast<FLOAT>(right);
}
template<> __host__ __inline__ DOUBLE _mulCHost<DOUBLE, FLOAT>(const DOUBLE& left, const FLOAT& right)
{
    return left * static_cast<DOUBLE>(right);
}

template<> __host__ __inline__ Real _mulCHost<Real, CLGComplex>(const Real& left, const CLGComplex& right)
{
    return left * _cuCabsf(right);
}
template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    return cuCmulf_cr(left, right);
}


#if _CLG_DOUBLEFLOAT
template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, FLOAT>(const CLGComplex& left, const FLOAT& right)
{
    return cuCmulf_cr(left, static_cast<Real>(right));
}
template<> __host__ __inline__ cuComplex _mulCHost<cuComplex, FLOAT>(const cuComplex& left, const FLOAT& right)
{
    return _cToFloat(cuCmulf_cr(_cToRealC(left), static_cast<Real>(right)));
}
template<> __host__ __inline__ cuComplex _mulCHost<cuComplex, DOUBLE>(const cuComplex& left, const DOUBLE& right)
{
    return _cToFloat(cuCmulf_cr(_cToRealC(left), static_cast<Real>(right)));
}
#else
template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, DOUBLE>(const CLGComplex& left, const DOUBLE& right)
{
    return cuCmulf_cr(left, static_cast<Real>(right));
}
template<> __host__ __inline__ cuDoubleComplex _mulCHost<cuDoubleComplex, FLOAT>(const cuDoubleComplex& left, const FLOAT& right)
{
    return cuCmulf_cd(left, static_cast<DOUBLE>(right));
}
template<> __host__ __inline__ cuDoubleComplex _mulCHost<cuDoubleComplex, DOUBLE>(const cuDoubleComplex& left, const DOUBLE& right)
{
    return cuCmulf_cd(left, right);
}
#endif

template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(left, right);
}
#if _CLG_DOUBLEFLOAT
template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, cuComplex>(const CLGComplex& left, const cuComplex& right)
{
    return _cuCmulf(left, _cToRealC(right));
}
#else
template<> __host__ __inline__ CLGComplex _mulCHost<CLGComplex, cuDoubleComplex>(const CLGComplex& left, const cuDoubleComplex& right)
{
    return _cuCmulf(left, _cToRealC(right));
}
#endif

#if _CLG_DOUBLEFLOAT
template<> __host__ __inline__ cuComplex _mulCHost<cuComplex, CLGComplex>(const cuComplex& left, const CLGComplex& right)
{
    return cuCmulf(left, _cToFloat(right));
}
template<> __host__ __inline__ cuComplex _mulCHost<cuComplex, cuComplex>(const cuComplex& left, const cuComplex& right)
{
    return cuCmulf(left, right);
}
#else
template<> __host__ __inline__ cuDoubleComplex _mulCHost<cuDoubleComplex, CLGComplex>(const cuDoubleComplex& left, const CLGComplex& right)
{
    return cuCmul(left, _cToDouble(right));
}
template<> __host__ __inline__ cuDoubleComplex _mulCHost<cuDoubleComplex, cuDoubleComplex>(const cuDoubleComplex& left, const cuDoubleComplex& right)
{
    return cuCmul(left, right);
}
#endif

template<typename T> __device__ __inline__ T _dagmulC(const T& left, const T& right) = delete;
template<typename T> __device__ __inline__ void _dagmul(T& left, const T& right) = delete;
template<> __device__ __inline__ CLGComplex _dagmulC<CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(_cuConjf(left), right);
}
template<> __device__ __inline__ void _dagmul<CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left = _cuCmulf(_cuConjf(left), right);
}

#define __DEFINE_TWO_ELEMENT_Func2(TYPENAME, FUNC1, FUNC2) \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME>(const TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2##C(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME>(TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2(right); \
} \


__DEFINE_TWO_ELEMENT_Func2(deviceSU2Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU2, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3, _dagmul, DaggerMul)

template<INT N> __device__ __inline__ void _dagmul(deviceZN<N>& left, const deviceZN<N>& right) { left.DaggerMul(right); }
template<INT N> __device__ __inline__ deviceZN<N> _dagmulC(const deviceZN<N>& left, const deviceZN<N>& right) { return left.DaggerMulC(right); }

template<typename T> __device__ __inline__ T _muldagC(const T& left, const T& right) = delete;
template<typename T> __device__ __inline__ void _muldag(T& left, const T& right) = delete;
template<> __device__ __inline__ Real _muldagC<Real>(const Real& left, const Real& right)
{
    return left * right;
}
template<> __device__ __inline__ void _muldag<Real>(Real& left, const Real& right)
{
    left = left * right;
}
template<> __device__ __inline__ CLGComplex _muldagC<CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(left, _cuConjf(right));
}
template<> __device__ __inline__ void _muldag<CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left = _cuCmulf(left, _cuConjf(right));
}

__DEFINE_TWO_ELEMENT_Func2(deviceSU2Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceWilsonVectorSU3, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU2, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3, _muldag, MulDagger)

template<INT N> __device__ __inline__ void _muldag(deviceZN<N>& left, const deviceZN<N>& right) { left.MulDagger(right); }
template<INT N> __device__ __inline__ deviceZN<N> _muldagC(const deviceZN<N>& left, const deviceZN<N>& right) { return left.MulDaggerC(right); }

template<typename T> __device__ __inline__ void _mul(T& left, const T& right, CLGComplex* buffer) = delete;
template<typename T> __device__ __inline__ void _muldag(T& left, const T& right, CLGComplex* buffer) = delete;

#define _DEFINE_BUFFER_MUL(TYPENAME) \
template<> __device__ __inline__ void _mul<TYPENAME>(TYPENAME& left, const TYPENAME& right, CLGComplex* buffer) \
{ \
    left.Mul(right, buffer); \
} \
template<> __device__ __inline__ void _muldag<TYPENAME>(TYPENAME& left, const TYPENAME& right, CLGComplex* buffer) \
{ \
    left.MulDagger(right, buffer); \
}

template<> __device__ __inline__ void _mul<CLGComplex>(CLGComplex& left, const CLGComplex& right, CLGComplex* buffer)
{ 
    _mul(left, right); 
} 
template<> __device__ __inline__ void _muldag<CLGComplex>(CLGComplex& left, const CLGComplex& right, CLGComplex* buffer)
{ 
    _muldag(left, right);
}
_DEFINE_BUFFER_MUL(deviceSU2);
_DEFINE_BUFFER_MUL(deviceSU3);

template<typename TMatrix, typename TVector> __device__ __inline__ TVector _mulVec(const TMatrix& matrix, const TVector& vector) = delete;

template<> __device__ __inline__ Real _mulVec<Real, Real>(const Real& matrix, const Real& vector)
{
    return matrix * vector;
}

template<> __device__ __inline__ CLGComplex _mulVec<CLGComplex, CLGComplex>(const CLGComplex& matrix, const CLGComplex& vector)
{
    return _cuCmulf(matrix, vector);
}

template<> __device__ __inline__ deviceSU2Vector _mulVec<deviceSU2, deviceSU2Vector>(const deviceSU2& matrix, const deviceSU2Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceSU3Vector _mulVec<deviceSU3, deviceSU3Vector>(const deviceSU3& matrix, const deviceSU3Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _mulVec<deviceSU3, deviceWilsonVectorSU3>(const deviceSU3& matrix, const deviceWilsonVectorSU3& vector)
{
    return matrix.MulWilsonVector(vector);
}

template<typename TMatrix, typename TVector> __device__ __inline__ TVector _dagmulVec(const TMatrix& matrix, const TVector& vector) = delete;

template<> __device__ __inline__ Real _dagmulVec<Real, Real>(const Real& matrix, const Real& vector)
{
    return matrix * vector;
}

template<> __device__ __inline__ CLGComplex _dagmulVec<CLGComplex, CLGComplex>(const CLGComplex& matrix, const CLGComplex& vector)
{
    return _cuCmulf(_cuConjf(matrix), vector);
}

template<> __device__ __inline__ deviceSU2Vector _dagmulVec<deviceSU2, deviceSU2Vector>(const deviceSU2& matrix, const deviceSU2Vector& vector)
{
    return matrix.DagMulVector(vector);
}

template<> __device__ __inline__ deviceSU3Vector _dagmulVec<deviceSU3, deviceSU3Vector>(const deviceSU3& matrix, const deviceSU3Vector& vector)
{
    return matrix.DagMulVector(vector);
}

template<> __device__ __inline__ deviceWilsonVectorSU3 _dagmulVec<deviceSU3, deviceWilsonVectorSU3>(const deviceSU3& matrix, const deviceWilsonVectorSU3& vector)
{
    return matrix.DaggerC().MulWilsonVector(vector);
}

template<typename TMatrix, typename TVector> __device__ __inline__ TVector _transposemulVec(const TMatrix& matrix, const TVector& vector) = delete;

template<typename TMatrix> __device__ __inline__ void _ta(TMatrix& matrix) = delete;

template<> __device__ __inline__ void _ta<Real>(Real& matrix)
{
    //do nothing
}

template<> __device__ __inline__ void _ta<CLGComplex>(CLGComplex& matrix)
{
    matrix.x = F(0.0);
    //matrix = _make_cuComplex(F(0.0), __cuCargf(matrix));
}

template<> __device__ __inline__ void _ta<deviceSU2>(deviceSU2& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU3>(deviceSU3& matrix)
{
    matrix.Ta();
}

template<INT N, INT NoE> __device__ __inline__ void _ta(deviceSUN<N, NoE>& matrix)
{
    matrix.Ta();
}

template<INT N> __device__ __inline__ void _ta(deviceZN<N>& matrix)
{
    matrix.Ta();
}

template<INT N, INT NoE> __device__ __inline__ void _ta(deviceSLNC<N, NoE>& matrix)
{
    matrix.Ta();
}

template<INT N, INT NoE> __device__ __inline__ void _ta(deviceUN<N, NoE>& matrix)
{
    matrix.Ta();
}

template<INT N, INT NoE> __device__ __inline__ void _ta(deviceON<N, NoE>& matrix)
{
    matrix.Ta();
}

template<INT N, INT NoE> __device__ __inline__ void _ta(deviceSON<N, NoE>& matrix)
{
    matrix.Ta();
}

template<typename TMatrix> __device__ __inline__ void _traceless(TMatrix& matrix) = delete;

template<INT N, INT NoE> __device__ __inline__ void _traceless(deviceSLNC<N, NoE>& matrix)
{
    matrix.Traceless();
}

template<typename TMatrix> __device__ __inline__ void _th(TMatrix& matrix) = delete;

template<> __device__ __inline__ void _th<Real>(Real& matrix)
{
    //do nothing
}

template<> __device__ __inline__ void _th<CLGComplex>(CLGComplex& matrix)
{
    matrix.y = F(0.0);
    //matrix = _make_cuComplex(F(0.0), __cuCargf(matrix));
}

template<> __device__ __inline__ void _th<deviceSU2>(deviceSU2& matrix)
{
    matrix.Th();
}

template<> __device__ __inline__ void _th<deviceSU3>(deviceSU3& matrix)
{
    matrix.Th();
}

template<INT N, INT NoE> __device__ __inline__ void _th(deviceSUN<N, NoE>& matrix)
{
    matrix.Th();
}

template<INT N> __device__ __inline__ void _th(deviceZN<N>& matrix)
{
    matrix.Th();
}


#define _impl_sub_add_mul(n) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n##Vector, _add, Add) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n, _add, Add) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n##Vector, _sub, Sub) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n, _sub, Sub) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n##Vector, _mul, Mul) \
__DEFINE_TWO_ELEMENT_Func(deviceSU##n, _mul, Mul) \
__DEFINE_TWO_ELEMENT_Func2(deviceSU##n##Vector, _dagmul, DaggerMul) \
__DEFINE_TWO_ELEMENT_Func2(deviceSU##n, _dagmul, DaggerMul) \
__DEFINE_TWO_ELEMENT_Func2(deviceSU##n##Vector, _muldag, MulDagger) \
__DEFINE_TWO_ELEMENT_Func2(deviceSU##n, _muldag, MulDagger) \
_DEFINE_BUFFER_MUL(deviceSU##n); \
template<> __device__ __inline__ deviceSU##n##Vector _mulVec<deviceSU##n, deviceSU##n##Vector>(const deviceSU##n& matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.MulVector(vector); \
} \
template<> __device__ __inline__ deviceSU##n##Vector _dagmulVec<deviceSU##n, deviceSU##n##Vector>(const deviceSU##n& matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.DagMulVector(vector); \
} 

#define _make_all_sub_add_mul(n) _DEF_F_N(n, _impl_sub_add_mul)
_make_all_sub_add_mul(_MAX_SUN)

#define _impl_slnc_sub_add_mul(n) \
__DEFINE_TWO_ELEMENT_Func(deviceSL##n##C, _add, Add) \
__DEFINE_TWO_ELEMENT_Func(deviceSL##n##C, _sub, Sub) \
__DEFINE_TWO_ELEMENT_Func(deviceSL##n##C, _mul, Mul) \
__DEFINE_TWO_ELEMENT_Func2(deviceSL##n##C, _dagmul, DaggerMul) \
__DEFINE_TWO_ELEMENT_Func2(deviceSL##n##C, _muldag, MulDagger) \
_DEFINE_BUFFER_MUL(deviceSL##n##C); \
template<> __device__ __inline__ deviceSU##n##Vector _mulVec<deviceSL##n##C, deviceSU##n##Vector>(const deviceSL##n##C & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.MulVector(vector); \
} \
template<> __device__ __inline__ deviceSU##n##Vector _dagmulVec<deviceSL##n##C, deviceSU##n##Vector>(const deviceSL##n##C & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.DagMulVector(vector); \
} 

#define _make_all_sub_add_mul_slnc(n) _DEF_FSLNC_N(n, _impl_slnc_sub_add_mul)
_make_all_sub_add_mul_slnc(_MAX_SLNC)

#define _impl_un_sub_add_mul(n) \
__DEFINE_TWO_ELEMENT_Func(deviceU##n, _add, Add) \
__DEFINE_TWO_ELEMENT_Func(deviceU##n, _sub, Sub) \
__DEFINE_TWO_ELEMENT_Func(deviceU##n, _mul, Mul) \
__DEFINE_TWO_ELEMENT_Func2(deviceU##n, _dagmul, DaggerMul) \
__DEFINE_TWO_ELEMENT_Func2(deviceU##n, _muldag, MulDagger) \
_DEFINE_BUFFER_MUL(deviceU##n); \
template<> __device__ __inline__ deviceSU##n##Vector _mulVec<deviceU##n, deviceSU##n##Vector>(const deviceU##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.MulVector(vector); \
} \
template<> __device__ __inline__ deviceSU##n##Vector _dagmulVec<deviceU##n, deviceSU##n##Vector>(const deviceU##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.DagMulVector(vector); \
}

#define _make_all_sub_add_mul_un(n) _DEF_FUN_N(n, _impl_un_sub_add_mul)
_make_all_sub_add_mul_un(_MAX_UN)

#define _impl_on_sub_add_mul(n) \
template<> __device__ __inline__ void _add<deviceO##n, deviceO##n>(deviceO##n& left, const deviceO##n& right) { left.Add(right); } \
template<> __device__ __inline__ void _add<deviceO##n, Real>(deviceO##n& left, const Real& right) { left.AddReal(right); } \
template<> __device__ __inline__ deviceO##n _addC<deviceO##n, deviceO##n>(const deviceO##n& left, const deviceO##n& right) { return left.AddC(right); } \
template<> __device__ __inline__ deviceO##n _addC<deviceO##n, Real>(const deviceO##n& left, const Real& right) { return left.AddRealC(right); } \
template<> __device__ __inline__ void _sub<deviceO##n, deviceO##n>(deviceO##n& left, const deviceO##n& right) { left.Sub(right); } \
template<> __device__ __inline__ void _sub<deviceO##n, Real>(deviceO##n& left, const Real& right) { left.SubReal(right); } \
template<> __device__ __inline__ deviceO##n _subC<deviceO##n, deviceO##n>(const deviceO##n& left, const deviceO##n& right) { return left.SubC(right); } \
template<> __device__ __inline__ deviceO##n _subC<deviceO##n, Real>(const deviceO##n& left, const Real& right) { return left.SubRealC(right); } \
template<> __device__ __inline__ void _mul<deviceO##n, deviceO##n>(deviceO##n& left, const deviceO##n& right) { left.Mul(right); } \
template<> __device__ __inline__ void _mul<deviceO##n, Real>(deviceO##n& left, const Real& right) { left.MulReal(right); } \
template<> __device__ __inline__ deviceO##n _mulC<deviceO##n, deviceO##n>(const deviceO##n& left, const deviceO##n& right) { return left.MulC(right); } \
template<> __device__ __inline__ deviceO##n _mulC<deviceO##n, Real>(const deviceO##n& left, const Real& right) { return left.MulRealC(right); } \
template<> __device__ __inline__ deviceSU##n##Vector _mulVec<deviceO##n, deviceSU##n##Vector>(const deviceO##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.MulVector(vector); \
} \
template<> __device__ __inline__ deviceSU##n##Vector _transposemulVec<deviceO##n, deviceSU##n##Vector>(const deviceO##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.TransposeMulVector(vector); \
} \
template<> __device__ __inline__ void _mul<deviceO##n, CLGComplex>(deviceO##n& left, const CLGComplex& right) { left.MulReal(right.x); } \
template<> __device__ __inline__ deviceO##n _mulC<deviceO##n, CLGComplex>(const deviceO##n& left, const CLGComplex& right) { return left.MulRealC(right.x); } \
template<> __device__ __inline__ void _add<deviceO##n, CLGComplex>(deviceO##n& left, const CLGComplex& right) { left.AddReal(right.x); } \
template<> __device__ __inline__ deviceO##n _addC<deviceO##n, CLGComplex>(const deviceO##n& left, const CLGComplex& right) { return left.AddRealC(right.x); } \
template<> __device__ __inline__ void _sub<deviceO##n, CLGComplex>(deviceO##n& left, const CLGComplex& right) { left.SubReal(right.x); } \
template<> __device__ __inline__ deviceO##n _subC<deviceO##n, CLGComplex>(const deviceO##n& left, const CLGComplex& right) { return left.SubRealC(right.x); } \
template<> __device__ __inline__ void _dagmul<deviceO##n>(deviceO##n& left, const deviceO##n& right) { left.TransposeMul(right); } \
template<> __device__ __inline__ deviceO##n _dagmulC<deviceO##n>(const deviceO##n& left, const deviceO##n& right) { return left.TransposeMulC(right); } \
template<> __device__ __inline__ void _muldag<deviceO##n>(deviceO##n& left, const deviceO##n& right) { left.MulTranspose(right); } \
template<> __device__ __inline__ deviceO##n _muldagC<deviceO##n>(const deviceO##n& left, const deviceO##n& right) { return left.MulTransposeC(right); }

#define _make_all_sub_add_mul_on(n) _DEF_FON_N(n, _impl_on_sub_add_mul)
_make_all_sub_add_mul_on(_MAX_ON)

#define _impl_son_sub_add_mul(n) \
template<> __device__ __inline__ void _add<deviceSO##n, deviceSO##n>(deviceSO##n& left, const deviceSO##n& right) { left.Add(right); } \
template<> __device__ __inline__ void _add<deviceSO##n, Real>(deviceSO##n& left, const Real& right) { left.AddReal(right); } \
template<> __device__ __inline__ deviceSO##n _addC<deviceSO##n, deviceSO##n>(const deviceSO##n& left, const deviceSO##n& right) { return left.AddC(right); } \
template<> __device__ __inline__ deviceSO##n _addC<deviceSO##n, Real>(const deviceSO##n& left, const Real& right) { return left.AddRealC(right); } \
template<> __device__ __inline__ void _sub<deviceSO##n, deviceSO##n>(deviceSO##n& left, const deviceSO##n& right) { left.Sub(right); } \
template<> __device__ __inline__ void _sub<deviceSO##n, Real>(deviceSO##n& left, const Real& right) { left.SubReal(right); } \
template<> __device__ __inline__ deviceSO##n _subC<deviceSO##n, deviceSO##n>(const deviceSO##n& left, const deviceSO##n& right) { return left.SubC(right); } \
template<> __device__ __inline__ deviceSO##n _subC<deviceSO##n, Real>(const deviceSO##n& left, const Real& right) { return left.SubRealC(right); } \
template<> __device__ __inline__ void _mul<deviceSO##n, deviceSO##n>(deviceSO##n& left, const deviceSO##n& right) { left.Mul(right); } \
template<> __device__ __inline__ void _mul<deviceSO##n, Real>(deviceSO##n& left, const Real& right) { left.MulReal(right); } \
template<> __device__ __inline__ deviceSO##n _mulC<deviceSO##n, deviceSO##n>(const deviceSO##n& left, const deviceSO##n& right) { return left.MulC(right); } \
template<> __device__ __inline__ deviceSO##n _mulC<deviceSO##n, Real>(const deviceSO##n& left, const Real& right) { return left.MulRealC(right); } \
template<> __device__ __inline__ deviceSU##n##Vector _mulVec<deviceSO##n, deviceSU##n##Vector>(const deviceSO##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.MulVector(vector); \
} \
template<> __device__ __inline__ deviceSU##n##Vector _transposemulVec<deviceSO##n, deviceSU##n##Vector>(const deviceSO##n & matrix, const deviceSU##n##Vector& vector) \
{ \
    return matrix.TransposeMulVector(vector); \
} \
template<> __device__ __inline__ void _mul<deviceSO##n, CLGComplex>(deviceSO##n& left, const CLGComplex& right) { left.MulReal(right.x); } \
template<> __device__ __inline__ deviceSO##n _mulC<deviceSO##n, CLGComplex>(const deviceSO##n& left, const CLGComplex& right) { return left.MulRealC(right.x); } \
template<> __device__ __inline__ void _add<deviceSO##n, CLGComplex>(deviceSO##n& left, const CLGComplex& right) { left.AddReal(right.x); } \
template<> __device__ __inline__ deviceSO##n _addC<deviceSO##n, CLGComplex>(const deviceSO##n& left, const CLGComplex& right) { return left.AddRealC(right.x); } \
template<> __device__ __inline__ void _sub<deviceSO##n, CLGComplex>(deviceSO##n& left, const CLGComplex& right) { left.SubReal(right.x); } \
template<> __device__ __inline__ deviceSO##n _subC<deviceSO##n, CLGComplex>(const deviceSO##n& left, const CLGComplex& right) { return left.SubRealC(right.x); } \
template<> __device__ __inline__ void _dagmul<deviceSO##n>(deviceSO##n& left, const deviceSO##n& right) { left.TransposeMul(right); } \
template<> __device__ __inline__ deviceSO##n _dagmulC<deviceSO##n>(const deviceSO##n& left, const deviceSO##n& right) { return left.TransposeMulC(right); } \
template<> __device__ __inline__ void _muldag<deviceSO##n>(deviceSO##n& left, const deviceSO##n& right) { left.MulTranspose(right); } \
template<> __device__ __inline__ deviceSO##n _muldagC<deviceSO##n>(const deviceSO##n& left, const deviceSO##n& right) { return left.MulTransposeC(right); }

#define _make_all_sub_add_mul_son(n) _DEF_FSON_N(n, _impl_son_sub_add_mul)
_make_all_sub_add_mul_son(_MAX_SON)

template<typename T> __device__ __inline__ CLGComplex _dot(const T& x, const T& y) = delete;

template<> __device__ __inline__ CLGComplex _dot<Real>(const Real& x, const Real& y)
{
    return _make_cuComplex(x * y, F(0.0));
}
template<> __device__ __inline__ CLGComplex _dot<CLGComplex>(const CLGComplex& x, const CLGComplex& y)
{
    return _cuCmulf(_cuConjf(x), y);
}
template<> __device__ __inline__ CLGComplex _dot<deviceSU2Vector>(const deviceSU2Vector& x, const deviceSU2Vector& y)
{
    return x.ConjugateDotC(y);
}
template<> __device__ __inline__ CLGComplex _dot<deviceSU3Vector>(const deviceSU3Vector& x, const deviceSU3Vector& y)
{
    return x.ConjugateDotC(y);
}
template<> __device__ __inline__ CLGComplex _dot<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& x, const deviceWilsonVectorSU3& y)
{
    return x.ConjugateDotC(y);
}
template<INT N, INT NoE> 
__device__ __inline__ CLGComplex _dot(const deviceSUNVector<N, NoE>& x, const deviceSUNVector<N, NoE>& y)
{
    return x.ConjugateDotC(y);
}

template<> __device__ __inline__ CLGComplex _dot<deviceSU2>(const deviceSU2& x, const deviceSU2& y)
{
    return x.DaggerMulC(y).Tr();
}

template<> __device__ __inline__ CLGComplex _dot<deviceSU3>(const deviceSU3& x, const deviceSU3& y)
{
    return x.DaggerMulC(y).Tr();
}

template<INT N, INT NoE>
__device__ __inline__ CLGComplex _dot(const deviceSUN<N, NoE>& x, const deviceSUN<N, NoE>& y)
{
    return x.DaggerMulC(y).Tr();
}

template<INT N, INT NoE>
__device__ __inline__ CLGComplex _dot(const deviceSLNC<N, NoE>& x, const deviceSLNC<N, NoE>& y)
{
    return x.DaggerMulC(y).Tr();
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _dot(const deviceUN<N, NoE>& x, const deviceUN<N, NoE>& y)
{
    return x.DaggerMulC(y).Tr();
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _dot(const deviceON<N, NoE>& x, const deviceON<N, NoE>& y)
{
    return _make_cuComplex(x.TransposeMulC(y).Tr(), F(0.0));
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _dot(const deviceSON<N, NoE>& x, const deviceSON<N, NoE>& y)
{
    return _make_cuComplex(x.TransposeMulC(y).Tr(), F(0.0));
}

template<INT N>
__device__ __inline__ CLGComplex _dot(const deviceZN<N>& x, const deviceZN<N>& y)
{
    return x.DaggerMulC(y).Tr();
}

template<typename T> __device__ __inline__ Real _lensq(const T& x) = delete;

template<> __device__ __inline__ Real _lensq<Real>(const Real& x)
{
    return x * x;
}
template<> __device__ __inline__ Real _lensq<CLGComplex>(const CLGComplex& x)
{
    return __cuCabsSqf(x);
}
template<> __device__ __inline__ Real _lensq<deviceSU2Vector>(const deviceSU2Vector& x)
{
    return x.ConjugateDotC(x).x;
}
template<> __device__ __inline__ Real _lensq<deviceSU3Vector>(const deviceSU3Vector& x)
{
    return x.ConjugateDotC(x).x;
}
template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceSUNVector<N, NoE>& x)
{
    return x.ConjugateDotC(x).x;
}

template<> __device__ __inline__ Real _lensq<deviceSU2>(const deviceSU2& x)
{
    return x.DaggerMulC(x).ReTr();
}

template<> __device__ __inline__ Real _lensq<deviceSU3>(const deviceSU3& x)
{
    return x.DaggerMulC(x).ReTr();
}

template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceSUN<N, NoE>& x)
{
    return x.DaggerMulC(x).ReTr();
}

template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceSLNC<N, NoE>& x)
{
    return x.DaggerMulC(x).ReTr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceUN<N, NoE>& x)
{
    return x.DaggerMulC(x).ReTr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceON<N, NoE>& x)
{
    return x.TransposeMulC(x).Tr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _lensq(const deviceSON<N, NoE>& x)
{
    return x.TransposeMulC(x).Tr();
}

template<INT N>
__device__ __inline__ Real _lensq(const deviceZN<N>& x)
{
    return x.DaggerMulC(x).ReTr();
}

template<> __device__ __inline__ Real _lensq<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& x)
{
    return x.ConjugateDotC(x).x;
}

template<typename T> __device__ __inline__  Real _retr(const T& x) = delete;

template<> __device__ __inline__ Real _retr<Real>(const Real& x)
{
    return x;
}

template<> __device__ __inline__ Real _retr<CLGComplex>(const CLGComplex& x)
{
    return x.x;
}

template<> __device__ __inline__ Real _retr<deviceSU2>(const deviceSU2& x)
{
    return x.ReTr();
}

template<> __device__ __inline__ Real _retr<deviceSU3>(const deviceSU3& x)
{
    return x.ReTr();
}

template<INT N, INT NoE>
__device__ __inline__ Real _retr(const deviceSUN<N, NoE>& x)
{
    return x.ReTr();
}

template<INT N, INT NoE>
__device__ __inline__ Real _retr(const deviceSLNC<N, NoE>& x)
{
    return x.ReTr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _retr(const deviceUN<N, NoE>& x)
{
    return x.ReTr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _retr(const deviceON<N, NoE>& x)
{
    return x.Tr();
}
template<INT N, INT NoE>
__device__ __inline__ Real _retr(const deviceSON<N, NoE>& x)
{
    return x.Tr();
}

template<INT N>
__device__ __inline__ Real _retr(const deviceZN<N>& x)
{
    return x.ReTr();
}

template<typename T> __device__ __inline__  CLGComplex _tr(const T& x) = delete;

template<> __device__ __inline__ CLGComplex _tr<CLGComplex>(const CLGComplex& x)
{
    return x;
}

template<> __device__ __inline__ CLGComplex _tr<deviceSU2>(const deviceSU2& x)
{
    return x.Tr();
}

template<> __device__ __inline__ CLGComplex _tr<deviceSU3>(const deviceSU3& x)
{
    return x.Tr();
}

template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceSUN<N, NoE>& x)
{
    return x.Tr();
}

template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceSLNC<N, NoE>& x)
{
    return x.Tr();
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceUN<N, NoE>& x)
{
    return x.Tr();
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceON<N, NoE>& x)
{
    return _make_cuComplex(x.Tr(), F(0.0));
}
template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceSON<N, NoE>& x)
{
    return _make_cuComplex(x.Tr(), F(0.0));
}

template<INT N>
__device__ __inline__ CLGComplex _tr(const deviceZN<N>& x)
{
    return x.Tr();
}

template<> __device__ __inline__ CLGComplex _tr<deviceSU2Vector>(const deviceSU2Vector& x)
{
    return x.Sum();
}

template<> __device__ __inline__ CLGComplex _tr<deviceSU3Vector>(const deviceSU3Vector& x)
{
    return x.Sum();
}

template<INT N, INT NoE>
__device__ __inline__ CLGComplex _tr(const deviceSUNVector<N, NoE>& x)
{
    return x.Sum();
}

template<> __device__ __inline__ CLGComplex _tr<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& x)
{
    return x.Sum();
}

template<typename T> __device__ __inline__ void _re(T& v) = delete;
template<> __device__ __inline__ void _re<CLGComplex>(CLGComplex& v) { v.y = F(0.0); }
template<> __device__ __inline__ void _re<deviceSU2Vector>(deviceSU2Vector& v) { v.Re(); }
template<> __device__ __inline__ void _re<deviceSU3Vector>(deviceSU3Vector& v) { v.Re(); }
template<INT N, INT NoVE> __device__ __inline__ void _re(deviceSUNVector<N, NoVE>& v) { v.Re(); }
template<INT N, INT NoE> __device__ __inline__ void _re(deviceSLNC<N, NoE>& v) { v.Re(); }
template<INT N, INT NoE> __device__ __inline__ void _re(deviceUN<N, NoE>& v) { v.Re(); }

template<typename T> __device__ __inline__ void _re2(T& v) = delete;
template<> __device__ __inline__ void _re2<CLGComplex>(CLGComplex& v) { v.y = F(0.0); v.x = F(2.0) * v.x; }
template<> __device__ __inline__ void _re2<deviceSU2>(deviceSU2& v) { v.Re2(); }
template<> __device__ __inline__ void _re2<deviceSU3>(deviceSU3& v) { v.Re2(); }
template<INT N, INT NoVE> __device__ __inline__ void _re2(deviceSUN<N, NoVE>& v) { v.Re2(); }
template<INT N> __device__ __inline__ void _re2(deviceZN<N>& v) { v.Re2(); }

template<typename T> __device__ __inline__ void _iim2(T& v) = delete;
template<> __device__ __inline__ void _iim2<CLGComplex>(CLGComplex& v) { v.y = F(2.0) * v.y;  v.x = F(0.0); }
template<> __device__ __inline__ void _iim2<deviceSU2>(deviceSU2& v) { v.iIm2(); }
template<> __device__ __inline__ void _iim2<deviceSU3>(deviceSU3& v) { v.iIm2(); }
template<INT N, INT NoVE> __device__ __inline__ void _iim2(deviceSUN<N, NoVE>& v) { v.iIm2(); }
template<INT N> __device__ __inline__ void _iim2(deviceZN<N>& v) { v.iIm2(); }

template<typename T> __device__ __inline__ Real _trim(const T& left, const T& right) = delete;
template<> __device__ __inline__ Real _trim<CLGComplex>(const CLGComplex& left, const CLGComplex& right) { return left.y * right.y; }
template<> __device__ __inline__ Real _trim<deviceSU2>(const deviceSU2& left, const deviceSU2& right) { return deviceSU2::TrIm(left, right); }
template<> __device__ __inline__ Real _trim<deviceSU3>(const deviceSU3& left, const deviceSU3& right) { return deviceSU3::TrIm(left, right); }
template<INT N, INT NoE>  __device__ __inline__ Real _trim(const deviceSUN<N, NoE>& left, const deviceSUN<N, NoE>& right) { return deviceSUN<N, NoE>::TrIm(left, right); }
template<INT N>  __device__ __inline__ Real _trim(const deviceZN<N>& left, const deviceZN<N>& right) { return deviceZN<N>::TrIm(left, right); }

template<typename T> __device__ __host__ __inline__  BYTE _dim() = delete;

template<> __device__ __host__ __inline__  BYTE _dim<Real>() { return 1; }
template<> __device__ __host__ __inline__  BYTE _dim<CLGComplex>() { return 1; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU2Vector>() { return 2; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU3Vector>() { return 3; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU2>() { return 2; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU3>() { return 3; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceWilsonVectorSU3>() { return 12; }


template<typename T> __device__ __host__ __inline__ UINT _generator_count_helper(T*) = delete;
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSUN<N, NoE>*) { return static_cast<UINT>(N * N - 1); }
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceUN<N, NoE>*) { return static_cast<UINT>(N * N); }
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceON<N, NoE>*) { return static_cast<UINT>(N * (N - 1) / 2); }
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSON<N, NoE>*) { return static_cast<UINT>(N * (N - 1) / 2); }
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSLNC<N, NoE>*) { return static_cast<UINT>(2 * (N * N - 1)); }
template<INT N, INT NoE> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSUNVector<N, NoE>*) { return 0; }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSU2*) { return static_cast<UINT>(2 * 2 - 1); }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSU3*) { return static_cast<UINT>(3 * 3 - 1); }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSU2Vector*) { return 0; }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(deviceSU3Vector*) { return 0; }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(deviceWilsonVectorSU3*) { return 0; }

template<> __device__ __host__ __inline__ UINT _generator_count_helper(CLGComplex*) { return 1; }
template<> __device__ __host__ __inline__ UINT _generator_count_helper(Real*) { return 1; }
template<INT N> __device__ __host__ __inline__ UINT _generator_count_helper(deviceZN<N>*) { return 1; }
template<INT N> __device__ __host__ __inline__ UINT _generator_count_helper(deviceDN<N>*) { return 2; }

template<typename T> __device__ __host__ __inline__ UINT _generator_count() { return _generator_count_helper((T*)NULL); }

template<typename T> __device__ __inline__  T _expreal(const T& x, Real a) = delete;

template<> __device__ __inline__  Real _expreal<Real>(const Real& x, Real a) { return _exp(x * a); }
template<> __device__ __inline__  CLGComplex _expreal<CLGComplex>(const CLGComplex& x, Real a) 
{ 
    const Real fAngle = x.y * a;
    return _make_cuComplex(_cos(fAngle), _sin(fAngle));
    //return __cuCexpf(cuCmulf_cr(x, a)); 
}
template<> __device__ __inline__  deviceSU2 _expreal<deviceSU2>(const deviceSU2& x, Real a) { return x.QuickExp(a); }

template<> __device__ __inline__  deviceSU3 _expreal<deviceSU3>(const deviceSU3& x, Real a)
{ 
    if (0 == _DC_ExpPrecision)
    {
        return x.QuickExp(a);
    }
    else if (1 == _DC_ExpPrecision)
    {
        return x.StrictExpTA(a);
    }
    return x.ExpReal(a, static_cast<BYTE>(_DC_ExpPrecision));
}

template<INT N, INT NoE>
__device__ __inline__  deviceSUN<N, NoE> _expreal(const deviceSUN<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}

template<INT N, INT NoE>
__device__ __inline__  deviceSLNC<N, NoE> _expreal(const deviceSLNC<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}
template<INT N, INT NoE>
__device__ __inline__  deviceUN<N, NoE> _expreal(const deviceUN<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}
template<INT N, INT NoE>
__device__ __inline__  deviceON<N, NoE> _expreal(const deviceON<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}
template<INT N, INT NoE>
__device__ __inline__  deviceSON<N, NoE> _expreal(const deviceSON<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}

template<INT N>
__device__ __inline__  deviceZN<N> _expreal(const deviceZN<N>& x, Real a)
{
    return x.ExpReal(a, N + 1);
}

template<typename T> __device__ __inline__  T _strictexp(const T& x) = delete;

template<> __device__ __inline__  CLGComplex _strictexp<CLGComplex>(const CLGComplex& x)
{
    return __cuCexpf(x);
}

template<> __device__ __inline__  deviceSU2 _strictexp<deviceSU2>(const deviceSU2& x)
{
    return x.StrictExp();
}

template<> __device__ __inline__  deviceSU3 _strictexp<deviceSU3>(const deviceSU3& x)
{
    return x.StrictExp();
}

template<INT N, INT NoE>
__device__ __inline__  deviceSUN<N, NoE> _strictexp(const deviceSUN<N, NoE>& x)
{
    return x.StrictExp();
}

template<INT N, INT NoE>
__device__ __inline__  deviceSLNC<N, NoE> _strictexp(const deviceSLNC<N, NoE>& x)
{
    return x.StrictExp();
}
template<INT N, INT NoE>
__device__ __inline__  deviceUN<N, NoE> _strictexp(const deviceUN<N, NoE>& x)
{
    return x.StrictExp();
}

template<INT N>
__device__ __inline__  deviceZN<N> _strictexp(const deviceZN<N>& x)
{
    return x.StrictExp();
}

template<INT N>
__device__ __inline__  deviceDN<N> _strictexp(const deviceDN<N>& x)
{
    return x.StrictExp();
}

template<INT N, INT NoE>
__device__ __inline__ deviceON<N, NoE> _strictexp(const deviceON<N, NoE>& x)
{
    return x.ExpReal(F(1.0));
}

template<INT N, INT NoE>
__device__ __inline__ deviceSON<N, NoE> _strictexp(const deviceSON<N, NoE>& x)
{
    return x.ExpReal(F(1.0));
}

template<typename T> __device__ __inline__  T _strictlog(const T& x) = delete;

template<> __device__ __inline__  CLGComplex _strictlog<CLGComplex>(const CLGComplex& x)
{
    return __cuClogf(x);
}

template<> __device__ __inline__  deviceSU2 _strictlog<deviceSU2>(const deviceSU2& x)
{
    return x.Log();
}

template<> __device__ __inline__  deviceSU3 _strictlog<deviceSU3>(const deviceSU3& x)
{
    return x.Log();
}

template<INT N, INT NoE> 
__device__ __inline__  deviceSUN<N, NoE> _strictlog(const deviceSUN<N, NoE>& x)
{
    return x.Log();
}

template<INT N, INT NoE>
__device__ __inline__  deviceSLNC<N, NoE> _strictlog(const deviceSLNC<N, NoE>& x)
{
    return x.Log();
}
template<INT N, INT NoE>
__device__ __inline__  deviceUN<N, NoE> _strictlog(const deviceUN<N, NoE>& x)
{
    return x.Log();
}
template<INT N, INT NoE>
__device__ __inline__  deviceON<N, NoE> _strictlog(const deviceON<N, NoE>& x)
{
    return x.Log();
}
template<INT N, INT NoE>
__device__ __inline__  deviceSON<N, NoE> _strictlog(const deviceSON<N, NoE>& x)
{
    return x.Log();
}

template<INT N>
__device__ __inline__  deviceZN<N> _strictlog(const deviceZN<N>& x)
{
    return x.Log();
}

template<INT N>
__device__ __inline__  deviceDN<N> _strictlog(const deviceDN<N>& x)
{
    return x.Log();
}

template<typename T> __device__ __inline__  void _norm(T& x) = delete;
template<> __device__ __inline__  void _norm<Real>(Real& x)
{
    //meaningless
    x = F(1.0);
}
template<> __device__ __inline__  void _norm<CLGComplex>(CLGComplex& x) 
{ 
    const Real fArg = __cuCargf(x);
    x = _make_cuComplex(_cos(fArg), _sin(fArg));
}
template<> __device__ __inline__  void _norm<deviceSU2>(deviceSU2& x) { x.Norm(); }
template<> __device__ __inline__  void _norm<deviceSU3>(deviceSU3& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceSUN<N, NoE>& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceSLNC<N, NoE>& x) { x.NormalizeDet(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceUN<N, NoE>& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceON<N, NoE>& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceSON<N, NoE>& x) { x.Norm(); }
template<INT N> __device__ __inline__  void _norm(deviceZN<N>& x) { x.Norm(); }
template<> __device__ __inline__  void _norm<deviceSU2Vector>(deviceSU2Vector& x) { x.Norm(); }
template<> __device__ __inline__  void _norm<deviceSU3Vector>(deviceSU3Vector& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceSUNVector<N, NoE>& x) { x.Norm(); }
template<> __device__ __inline__  void _norm<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& x) { x.Norm(); }

template<typename T> __device__ __host__ __inline__  Real _element(const T& x, INT idx) = delete;

template<> __device__ __host__ __inline__ Real _element<Real>(const Real& x, INT idx)
{
    return x;
}

template<> __device__ __host__ __inline__ Real _element<CLGComplex>(const CLGComplex& x, INT idx)
{
    if (0 == idx)
    {
        return x.x;
    }
    return x.y;
}

template<> __device__ __host__ __inline__ Real _element<deviceSU2Vector>(const deviceSU2Vector& x, INT idx)
{
    if (0 == idx)
    {
        return x.m_ve[0].x;
    }
    else if (1 == idx)
    {
        return x.m_ve[0].y;
    }
    else if (2 == idx)
    {
        return x.m_ve[1].x;
    }
    else if (3 == idx)
    {
        return x.m_ve[1].y;
    }

    return F(0.0);
}

template<> __device__ __host__ __inline__ Real _element<deviceSU3Vector>(const deviceSU3Vector& x, INT idx)
{
    if (0 == idx)
    {
        return x.m_ve[0].x;
    }
    else if (1 == idx)
    {
        return x.m_ve[0].y;
    }
    else if (2 == idx)
    {
        return x.m_ve[1].x;
    }
    else if (3 == idx)
    {
        return x.m_ve[1].y;
    }
    else if (4 == idx)
    {
        return x.m_ve[2].x;
    }
    else if (5 == idx)
    {
        return x.m_ve[2].y;
    }

    return F(0.0);
}

template<INT N, INT NoE> 
__device__ __host__ __inline__ Real _element(const deviceSUNVector<N, NoE>& x, INT idx)
{
    if (idx < 2 * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_ve[idxc].y;
        }
        return x.m_ve[idxc].x;
    }
    return F(0.0);
}

template<> __device__ __host__ __inline__ Real _element<deviceSU2>(const deviceSU2& x, INT idx)
{
    if (idx < 8)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}

template<> __device__ __host__ __inline__ Real _element<deviceSU3>(const deviceSU3& x, INT idx)
{
    if (idx < 18)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}

template<> __device__ __host__ __inline__ Real _element<deviceSU3_12>(const deviceSU3_12& x, INT idx)
{
    if (idx < 12)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}

template<> __device__ __host__ __inline__ Real _element<deviceWilsonVectorSU3>(const deviceWilsonVectorSU3& x, INT idx)
{
    if (idx < 24)
    {
        return x.m_rme[idx];
    }
    return F(0.0);
}

template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceSUN<N, NoE>& x, INT idx)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}

template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceSLNC<N, NoE>& x, INT idx)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}
template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceUN<N, NoE>& x, INT idx)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_me[idxc].y;
        }
        return x.m_me[idxc].x;
    }
    return F(0.0);
}
template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceON<N, NoE>& x, INT idx)
{
    if (idx < N * N)
    {
        return x.m_me[idx];
    }
    return F(0.0);
}
template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceSON<N, NoE>& x, INT idx)
{
    if (idx < N * N)
    {
        return x.m_me[idx];
    }
    return F(0.0);
}

template<INT N>
__device__ __host__ __inline__ Real _element(const deviceZN<N>& x, INT idx)
{
    if (0 == idx) { return x.m_me.x; }
    if (1 == idx) { return x.m_me.y; }
    return F(0.0);
}

template<typename T> __device__ __host__ __inline__  void _setelement(T& x, INT idx, Real v) = delete;

template<> __device__ __host__ __inline__ void _setelement<Real>(Real& x, INT idx, Real v)
{
    x = v;
}

template<> __device__ __host__ __inline__ void _setelement<CLGComplex>(CLGComplex& x, INT idx, Real v)
{
    if (0 == idx)
    {
        x.x = v;
        return;
    }
    x.y = v;
}

template<> __device__ __host__ __inline__ void _setelement<deviceSU2Vector>(deviceSU2Vector& x, INT idx, Real v)
{
    if (0 == idx)
    {
        x.m_ve[0].x = v;
    }
    else if (1 == idx)
    {
        x.m_ve[0].y = v;
    }
    else if (2 == idx)
    {
        x.m_ve[1].x = v;
    }
    else if (3 == idx)
    {
        x.m_ve[1].y = v;
    }
}

template<> __device__ __host__ __inline__ void _setelement<deviceSU3Vector>(deviceSU3Vector& x, INT idx, Real v)
{
    if (0 == idx)
    {
        x.m_ve[0].x = v;
    }
    else if (1 == idx)
    {
        x.m_ve[0].y = v;
    }
    else if (2 == idx)
    {
        x.m_ve[1].x = v;
    }
    else if (3 == idx)
    {
        x.m_ve[1].y = v;
    }
    else if (4 == idx)
    {
        x.m_ve[2].x = v;
    }
    else if (5 == idx)
    {
        x.m_ve[2].y = v;
    }
}

template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceSUNVector<N, NoE>& x, INT idx, Real v)
{
    if (idx < 2 * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_ve[idxc].y = v;
            return;
        }
        x.m_ve[idxc].x = v;
    }
}

template<> __device__ __host__ __inline__ void _setelement<deviceSU2>(deviceSU2& x, INT idx, Real v)
{
    if (idx < 8)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}

template<> __device__ __host__ __inline__ void _setelement<deviceSU3>(deviceSU3& x, INT idx, Real v)
{
    if (idx < 18)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}

template<> __device__ __host__ __inline__ void _setelement<deviceSU3_12>(deviceSU3_12& x, INT idx, Real v)
{
    if (idx < 12)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}

template<> __device__ __host__ __inline__ void _setelement<deviceWilsonVectorSU3>(deviceWilsonVectorSU3& x, INT idx, Real v)
{
    if (idx < 24)
    {
        x.m_rme[idx] = v;
    }
}

template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceSUN<N, NoE>& x, INT idx, Real v)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}

template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceSLNC<N, NoE>& x, INT idx, Real v)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}
template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceUN<N, NoE>& x, INT idx, Real v)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            x.m_me[idxc].y = v;
            return;
        }
        x.m_me[idxc].x = v;
    }
}
template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceON<N, NoE>& x, INT idx, Real v)
{
    if (idx < N * N)
    {
        x.m_me[idx] = v;
    }
}
template<INT N, INT NoE>
__device__ __host__ __inline__ void _setelement(deviceSON<N, NoE>& x, INT idx, Real v)
{
    if (idx < N * N)
    {
        x.m_me[idx] = v;
    }
}

template<INT N>
__device__ __host__ __inline__ void _setelement(deviceZN<N>& x, INT idx, Real v)
{
    if (0 == idx) { x.m_me.x = v; return; }
    if (1 == idx) { x.m_me.y = v; }
}

template<typename T> __device__ __host__ __inline__  WORD _elementdim() = delete;

template<> __device__ __host__ __inline__  WORD _elementdim<Real>() { return 1; }
template<> __device__ __host__ __inline__  WORD _elementdim<CLGComplex>() { return 2; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU2Vector>() { return 4; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU3Vector>() { return 6; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU2>() { return 8; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU3>() { return 18; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU3_12>() { return 12; }
template<> __device__ __host__ __inline__  WORD _elementdim<deviceWilsonVectorSU3>() { return 24; }


#define _impl_dim_func(n) \
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU##n##Vector>() { return n; } \
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU##n>() { return n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU##n##Vector>() { return 2 * n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSU##n>() { return 2 * n * n; } \


#define _make_all_dim_func(n) _DEF_F_N(n, _impl_dim_func)
_make_all_dim_func(_MAX_SUN)

#define _impl_dim_func_slnc(n) \
template<> __device__ __host__ __inline__  BYTE _dim<deviceSL##n##C>() { return n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSL##n##C>() { return 2 * n * n; } \

#define _make_all_dim_func_slnc(n) _DEF_FSLNC_N(n, _impl_dim_func_slnc)
_make_all_dim_func_slnc(_MAX_SLNC)

#define _impl_dim_func_un(n) \
template<> __device__ __host__ __inline__  BYTE _dim<deviceU##n>() { return n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceU##n>() { return 2 * n * n; } \

#define _make_all_dim_func_un(n) _DEF_FUN_N(n, _impl_dim_func_un)
_make_all_dim_func_un(_MAX_UN)

#define _impl_dim_func_on(n) \
template<> __device__ __host__ __inline__  BYTE _dim<deviceO##n>() { return n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceO##n>() { return n * n; } \

#define _make_all_dim_func_on(n) _DEF_FON_N(n, _impl_dim_func_on)
_make_all_dim_func_on(_MAX_ON)

#define _impl_dim_func_son(n) \
template<> __device__ __host__ __inline__  BYTE _dim<deviceSO##n>() { return n; } \
template<> __device__ __host__ __inline__  WORD _elementdim<deviceSO##n>() { return n * n; } \

#define _make_all_dim_func_son(n) _DEF_FSON_N(n, _impl_dim_func_son)
_make_all_dim_func_son(_MAX_SON)

template<typename T> __device__ __host__ __inline__ CLGComplex _vn(const T& v, INT idx) = delete;
template<> __device__ __host__ __inline__ CLGComplex _vn<Real>(const Real& v, INT idx)
{
    return _make_cuComplex(v, F(0.0));
}
template<> __device__ __host__ __inline__ CLGComplex _vn<CLGComplex>(const CLGComplex& v, INT idx)
{
    return v;
}
template<> __device__ __host__ __inline__ CLGComplex _vn<deviceSU2Vector>(const deviceSU2Vector& v, INT idx)
{
    return v.m_ve[idx];
}
template<> __device__ __host__ __inline__ CLGComplex _vn<deviceSU3Vector>(const deviceSU3Vector& v, INT idx)
{
    return v.m_ve[idx];
}
template<INT N, INT NoVE> __device__ __host__ __inline__ CLGComplex _vn(const deviceSUNVector<N, NoVE>& v, INT idx)
{
    return v.m_ve[idx];
}

template<INT N>
inline CCString appToString(const deviceZN<N>& v)
{
    CCString ret;
    ret.Format(_T("{%f %s %f I}"),
        v.m_me.x, v.m_me.y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me.y)
    );
    return ret;
}

template<>
inline CCString appToString<deviceSU2Vector>(const deviceSU2Vector& v)
{
    CCString ret;
    ret.Format(_T("{%f %s %f I, %f %s %f I}"), 
        v.m_ve[0].x, v.m_ve[0].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[0].y),
        v.m_ve[1].x, v.m_ve[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[1].y)
        );
    return ret;
}

template<>
inline CCString appToString<deviceSU3Vector>(const deviceSU3Vector& v)
{
    CCString ret;
    ret.Format(_T("{%f %s %f I, %f %s %f I, %f %s %f I}"),
        v.m_ve[0].x, v.m_ve[0].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[0].y),
        v.m_ve[1].x, v.m_ve[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[1].y),
        v.m_ve[2].x, v.m_ve[2].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[2].y)
    );
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceSUNVector<N, NoE>& v)
{
    CCString ret = _T("{");
    for (INT i = 0; i < N; ++i)
    {
        CCString stoadd;
        stoadd.Format(_T("%f %s %f I"), v.m_ve[i].x, v.m_ve[i].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_ve[i].y));
        ret = ret + stoadd;
        if (i != (N - 1))
        {
            ret = ret + _T(", ");
        }
        else
        {
            ret = ret + _T("}");
        }
    }
    return ret;
}

template<>
inline CCString appToString<deviceSU2>(const deviceSU2& v)
{
    CCString ret;
    ret.Format(_T("{{%f %s %f I, %f %s %f I}, {%f %s %f I, %f %s %f I}}"),
        v.m_me[0].x, v.m_me[0].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[0].y),
        v.m_me[1].x, v.m_me[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[1].y),
        v.m_me[2].x, v.m_me[2].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[2].y),
        v.m_me[3].x, v.m_me[3].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[3].y)
    );
    return ret;
}

template<INT N>
inline CCString appToString(const deviceDN<N>& v)
{
    CCString ret;
    ret.Format(_T("{{%f %s %f I, %f %s %f I}, {%f %s %f I, %f %s %f I}}"),
        v.m_me[0].x, v.m_me[0].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[0].y),
        v.m_me[1].x, v.m_me[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[1].y),
        v.m_me[2].x, v.m_me[2].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[2].y),
        v.m_me[3].x, v.m_me[3].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[3].y)
    );
    return ret;
}

template<>
inline CCString appToString<deviceSU3>(const deviceSU3& v)
{
    CCString ret;
    ret.Format(_T("{{%f %s %f I, %f %s %f I, %f %s %f I},\n {%f %s %f I, %f %s %f I, %f %s %f I},\n {%f %s %f I, %f %s %f I, %f %s %f I}}"),
        v.m_me[0].x, v.m_me[0].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[0].y),
        v.m_me[1].x, v.m_me[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[1].y),
        v.m_me[2].x, v.m_me[2].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[2].y),
        v.m_me[3].x, v.m_me[3].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[3].y),
        v.m_me[4].x, v.m_me[4].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[4].y),
        v.m_me[5].x, v.m_me[5].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[5].y),
        v.m_me[6].x, v.m_me[6].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[6].y),
        v.m_me[7].x, v.m_me[7].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[7].y),
        v.m_me[8].x, v.m_me[8].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[8].y)
    );
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceSUN<N, NoE>& v)
{
    CCString ret = _T("{{");
    for (INT y = 0; y < N; ++y)
    {
        for (INT x = 0; x < N; ++x)
        {
            CCString stoAdd;
            stoAdd.Format(_T("%f %s %f I"), v.m_me[y * N + x].x, v.m_me[y * N + x].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[y * N + x].y));
            ret = ret + stoAdd;
            if (x != (N - 1))
            {
                ret = ret + _T(", ");
            }
            else
            {
                if (y != (N - 1))
                {
                    ret = ret + _T("},\n {");
                }
            }
        }
    }
    ret = ret + _T("}}");
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceUN<N, NoE>& v)
{
    CCString ret = _T("{{");
    for (INT y = 0; y < N; ++y)
    {
        for (INT x = 0; x < N; ++x)
        {
            CCString stoAdd;
            stoAdd.Format(_T("%f %s %f I"), v.m_me[y * N + x].x, v.m_me[y * N + x].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[y * N + x].y));
            ret = ret + stoAdd;
            if (x != (N - 1))
            {
                ret = ret + _T(", ");
            }
            else
            {
                if (y != (N - 1))
                {
                    ret = ret + _T("},\n {");
                }
            }
        }
    }
    ret = ret + _T("}}");
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceSLNC<N, NoE>& v)
{
    CCString ret = _T("{{");
    for (INT y = 0; y < N; ++y)
    {
        for (INT x = 0; x < N; ++x)
        {
            CCString stoAdd;
            stoAdd.Format(_T("%f %s %f I"), v.m_me[y * N + x].x, v.m_me[y * N + x].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[y * N + x].y));
            ret = ret + stoAdd;
            if (x != (N - 1))
            {
                ret = ret + _T(", ");
            }
            else
            {
                if (y != (N - 1))
                {
                    ret = ret + _T("},\n {");
                }
            }
        }
    }
    ret = ret + _T("}}");
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceON<N, NoE>& v)
{
    CCString ret = _T("{{");
    for (INT y = 0; y < N; ++y)
    {
        for (INT x = 0; x < N; ++x)
        {
            CCString stoAdd;
            stoAdd.Format(_T("%f"), v.m_me[y * N + x]);
            ret = ret + stoAdd;
            if (x != (N - 1))
            {
                ret = ret + _T(", ");
            }
            else
            {
                if (y != (N - 1))
                {
                    ret = ret + _T("},\n {");
                }
            }
        }
    }
    ret = ret + _T("}}");
    return ret;
}

template<INT N, INT NoE>
inline CCString appToString(const deviceSON<N, NoE>& v)
{
    CCString ret = _T("{{");
    for (INT y = 0; y < N; ++y)
    {
        for (INT x = 0; x < N; ++x)
        {
            CCString stoAdd;
            stoAdd.Format(_T("%f"), v.m_me[y * N + x]);
            ret = ret + stoAdd;
            if (x != (N - 1))
            {
                ret = ret + _T(", ");
            }
            else
            {
                if (y != (N - 1))
                {
                    ret = ret + _T("},\n {");
                }
            }
        }
    }
    ret = ret + _T("}}");
    return ret;
}

template<typename T> __device__ __inline__  void _print(const T& x) = delete;
template<typename T> __device__ __inline__  void _print(const T& x, const char* head) = delete;

template<> __device__ __inline__  void _print<Real>(const Real& x)
{
    printf("%f\n", x);
}
template<> __device__ __inline__  void _print<CLGComplex>(const CLGComplex& x)
{
    printf("%f + %f I\n", x.x, x.y);
}
template<> __device__ __inline__  void _print<deviceSU2>(const deviceSU2& x)
{
    x.DebugPrint();
}
template<> __device__ __inline__  void _print<deviceSU2>(const deviceSU2& x, const char* head)
{
    x.DebugPrint(head);
}
template<> __device__ __inline__  void _print<deviceSU3>(const deviceSU3& x)
{
    x.DebugPrint();
}
template<> __device__ __inline__  void _print<deviceSU3>(const deviceSU3& x, const char* head)
{
    x.DebugPrint(head);
}
template<> __device__ __inline__  void _print<deviceSU2Vector>(const deviceSU2Vector& x)
{
    x.DebugPrint();
}
template<> __device__ __inline__  void _print<deviceSU3Vector>(const deviceSU3Vector& x)
{
    x.DebugPrint();
}
template<INT N, INT NoE> __device__ __inline__  void _print(const deviceSUN<N, NoE>& x)
{
    x.DebugPrint();
}
template<INT N, INT NoE> __device__ __inline__  void _print(const deviceSUN<N, NoE>& x, const char* head)
{
    x.DebugPrint(head);
}
template<INT N, INT NoVE> __device__ __inline__  void _print(const deviceSUNVector<N, NoVE>& x)
{
    x.DebugPrint();
}
template<INT N> __device__ __inline__  void _print(const deviceZN<N>& x)
{
    x.DebugPrint();
}
template<INT N> __device__ __inline__  void _print(const deviceZN<N>& x, const char* head)
{
    x.DebugPrint(head);
}

template<typename T> __device__ __inline__  CLGComplex _detv(const T& x) = delete;

template<> __device__ __inline__  CLGComplex _detv<Real>(const Real& x)
{
    return _make_cuComplex(x, F(0.0));
}
template<> __device__ __inline__  CLGComplex _detv<CLGComplex>(const CLGComplex& x)
{
    return x;
}
template<> __device__ __inline__  CLGComplex _detv<deviceSU2>(const deviceSU2& x)
{
    return x.Det();
}
template<> __device__ __inline__  CLGComplex _detv<deviceSU3>(const deviceSU3& x)
{
    return x.Det();
}
template<> __device__ __inline__  CLGComplex _detv<deviceSU2Vector>(const deviceSU2Vector& x)
{
    return _make_cuComplex(x.Abs(), F(0.0));
}
template<> __device__ __inline__  CLGComplex _detv<deviceSU3Vector>(const deviceSU3Vector& x)
{
    return _make_cuComplex(x.Abs(), F(0.0));
}
template<INT N, INT NoE> __device__ __inline__  CLGComplex _detv(const deviceSUN<N, NoE>& x)
{
    return x.Determinent();
}
template<INT N, INT NoVE> __device__ __inline__  CLGComplex _detv(const deviceSUNVector<N, NoVE>& x)
{
    return _make_cuComplex(x.Abs(), F(0.0));
}
template<INT N> __device__ __inline__  CLGComplex _detv(const deviceZN<N>& x)
{
    return x.Tr();
}

template<typename T> __device__ __inline__  Real _absv(const T& x) = delete;

template<> __device__ __inline__  Real _absv<Real>(const Real& x)
{
    return abs(x);
}
template<> __device__ __inline__  Real _absv<CLGComplex>(const CLGComplex& x)
{
    return _cuCabsf(x);
}
template<> __device__ __inline__  Real _absv<deviceSU2>(const deviceSU2& x)
{
    return _cuCabsf(x.Det());
}
template<> __device__ __inline__  Real _absv<deviceSU3>(const deviceSU3& x)
{
    return _cuCabsf(x.Det());
}
template<> __device__ __inline__  Real _absv<deviceSU2Vector>(const deviceSU2Vector& x)
{
    return x.Abs();
}
template<> __device__ __inline__  Real _absv<deviceSU3Vector>(const deviceSU3Vector& x)
{
    return x.Abs();
}
template<INT N, INT NoE> __device__ __inline__  Real _absv(const deviceSUN<N, NoE>& x)
{
    return _cuCabsf(x.Determinent());
}
template<INT N, INT NoVE> __device__ __inline__  Real _absv(const deviceSUNVector<N, NoVE>& x)
{
    return x.Abs();
}
template<INT N> __device__ __inline__  Real _absv(const deviceZN<N>& x)
{
    return _cuCabsf(x.Tr());
}

#pragma endregion

// ==============================
// deviceZN<N> specializations
// ==============================

#define _make_impl_zn(n) \
template<> __device__ __inline__ deviceZN<n> _makeId<deviceZN<n>>() { return deviceZN<n>::makeZNId(); } \
template<> __device__ __inline__ deviceZN<n> _makeZero<deviceZN<n>>() { return deviceZN<n>::makeZNZero(); } \
template<> __device__ __inline__ deviceZN<n> _makeAsK<deviceZN<n>>(UINT k) { return deviceZN<n>::makeAsK(k); } \
template<> __device__ __inline__ deviceZN<n> _makeRandom<deviceZN<n>>(UINT fatIdx) { return deviceZN<n>::makeZNRandom(fatIdx); } \
template<> __device__ __inline__ deviceZN<n> _makeGaussian<deviceZN<n>>(UINT fatIdx) { return deviceZN<n>::makeZNId(); } \
template<> __device__ __inline__ deviceZN<n> _makeSumGenerator<deviceZN<n>>(Real factor) { return deviceZN<n>::makeZNId(); }

#define _impl_dim_zn(n) \
template<> __device__ __host__ __inline__ BYTE _dim<deviceZN<n>>() { return 1; } \
template<> __device__ __host__ __inline__ WORD _elementdim<deviceZN<n>>() { return 2; }

#define _DEF_ZN_SPECIALIZATION(N, unuse) \
_make_impl_zn(N) \
_impl_dim_zn(N) 

#define _make_all_imp_zn(n) _DEF_F_ZN(n, _DEF_ZN_SPECIALIZATION)
_make_all_imp_zn(_MAX_ZN)


// ==============================
// deviceDN<N> specializations
// ==============================

#define _DEF_DN_MAKE(N) \
template<> __device__ __inline__ deviceDN<N> _makeId<deviceDN<N>>() { return deviceDN<N>::makeDNId(); } \
template<> __device__ __inline__ deviceDN<N> _makeZero<deviceDN<N>>() { return deviceDN<N>::makeDNZero(); } \
template<> __device__ __inline__ deviceDN<N> _makeRandom<deviceDN<N>>(UINT fatIdx) { return deviceDN<N>::makeDNRandom(fatIdx); } \
template<> __device__ __inline__ deviceDN<N> _makeAsK<deviceDN<N>>(UINT k) { return deviceDN<N>::makeAsK(k); } \
template<> __device__ __inline__ deviceDN<N> _makeGaussian<deviceDN<N>>(UINT fatIdx) { return deviceDN<N>::makeDNId(); } \
template<> __device__ __inline__ deviceDN<N> _makeSumGenerator<deviceDN<N>>(Real factor) { return deviceDN<N>::makeDNId(); }

#define _DEF_DN_DIM(N) \
template<> __device__ __host__ __inline__ BYTE _dim<deviceDN<N>>() { return 2; } \
template<> __device__ __host__ __inline__ WORD _elementdim<deviceDN<N>>() { return 8; }



#define _DEF_DN_SPECIALIZATION(N, unuse) \
_DEF_DN_MAKE(N) \
_DEF_DN_DIM(N) 

#define _make_all_imp_dn(n) _DEF_F_DN(n, _DEF_DN_SPECIALIZATION)
_make_all_imp_dn(_MAX_DN)

// --- creation ---

template<INT N> __device__ __inline__ void _Id(deviceDN<N>& v) { v.Id(); }
template<INT N> __device__ __inline__ void _Zero(deviceDN<N>& v) { v.Zero(); }

// --- dagger ---

template<INT N> __device__ __inline__ void _dagger(deviceDN<N>& element) { element.Dagger(); }
template<INT N> __device__ __inline__ deviceDN<N> _daggerC(const deviceDN<N>& element) { return element.DaggerC(); }

// --- opposite ---

template<INT N> __device__ __inline__ void _oppo(deviceDN<N>& element) { element.Opposite(); }
template<INT N> __device__ __inline__ deviceDN<N> _oppoC(const deviceDN<N>& element) { return element.OppositeC(); }

// --- add ---

template<INT N> __device__ __inline__ void _add(deviceDN<N>& left, const deviceDN<N>& right) { left.Add(right); }
template<INT N> __device__ __inline__ void _add(deviceDN<N>& left, const Real& right) { left.AddReal(right); }
template<INT N> __device__ __inline__ deviceDN<N> _addC(const deviceDN<N>& left, const deviceDN<N>& right) { return left.AddC(right); }
template<INT N> __device__ __inline__ deviceDN<N> _addC(const deviceDN<N>& left, const Real& right) { return left.AddRealC(right); }

// --- sub ---

template<INT N> __device__ __inline__ void _sub(deviceDN<N>& left, const deviceDN<N>& right) { left.Sub(right); }
template<INT N> __device__ __inline__ void _sub(deviceDN<N>& left, const Real& right) { left.SubReal(right); }
template<INT N> __device__ __inline__ deviceDN<N> _subC(const deviceDN<N>& left, const deviceDN<N>& right) { return left.SubC(right); }
template<INT N> __device__ __inline__ deviceDN<N> _subC(const deviceDN<N>& left, const Real& right) { return left.SubRealC(right); }

// --- mul ---

template<INT N> __device__ __inline__ void _mul(deviceDN<N>& left, const deviceDN<N>& right) { left.Mul(right); }
template<INT N> __device__ __inline__ void _mul(deviceDN<N>& left, const Real& right) { left.MulReal(right); }
template<INT N> __device__ __inline__ void _mul(deviceDN<N>& left, const CLGComplex& right) { left.MulComp(right); }
template<INT N> __device__ __inline__ deviceDN<N> _mulC(const deviceDN<N>& left, const deviceDN<N>& right) { return left.MulC(right); }
template<INT N> __device__ __inline__ deviceDN<N> _mulC(const deviceDN<N>& left, const Real& right) { return left.MulRealC(right); }
template<INT N> __device__ __inline__ deviceDN<N> _mulC(const deviceDN<N>& left, const CLGComplex& right) { return left.MulCompC(right); }

// --- dagmul / muldag ---

template<INT N> __device__ __inline__ void _dagmul(deviceDN<N>& left, const deviceDN<N>& right) { left.DaggerMul(right); }
template<INT N> __device__ __inline__ deviceDN<N> _dagmulC(const deviceDN<N>& left, const deviceDN<N>& right) { return left.DaggerMulC(right); }
template<INT N> __device__ __inline__ void _muldag(deviceDN<N>& left, const deviceDN<N>& right) { left.MulDagger(right); }
template<INT N> __device__ __inline__ deviceDN<N> _muldagC(const deviceDN<N>& left, const deviceDN<N>& right) { return left.MulDaggerC(right); }

// --- ta / th ---

template<INT N> __device__ __inline__ void _ta(deviceDN<N>& matrix) { matrix.Ta(); }
template<INT N> __device__ __inline__ void _th(deviceDN<N>& matrix) { matrix.Th(); }

// --- reduction ---

template<INT N> __device__ __inline__ CLGComplex _dot(const deviceDN<N>& x, const deviceDN<N>& y)
{
    return _cuCaddf(
        _cuCaddf(_cuCmulf(_cuConjf(x.m_me[0]), y.m_me[0]), _cuCmulf(_cuConjf(x.m_me[1]), y.m_me[1])),
        _cuCaddf(_cuCmulf(_cuConjf(x.m_me[2]), y.m_me[2]), _cuCmulf(_cuConjf(x.m_me[3]), y.m_me[3]))
    );
}

template<INT N> __device__ __inline__ Real _lensq(const deviceDN<N>& x)
{
    return _cuCabsf(x.m_me[0]) * _cuCabsf(x.m_me[0])
         + _cuCabsf(x.m_me[1]) * _cuCabsf(x.m_me[1])
         + _cuCabsf(x.m_me[2]) * _cuCabsf(x.m_me[2])
         + _cuCabsf(x.m_me[3]) * _cuCabsf(x.m_me[3]);
}

template<INT N> __device__ __inline__ Real _retr(const deviceDN<N>& x) { return x.ReTr(); }
template<INT N> __device__ __inline__ CLGComplex _tr(const deviceDN<N>& x) { return x.Tr(); }

// --- re2 / iim2 / trim ---

template<INT N> __device__ __inline__ void _re2(deviceDN<N>& v) { v.Re2(); }
template<INT N> __device__ __inline__ void _iim2(deviceDN<N>& v) { v.iIm2(); }
template<INT N> __device__ __inline__ Real _trim(const deviceDN<N>& left, const deviceDN<N>& right) { return deviceDN<N>::TrIm(left, right); }

// --- exponential / norm ---

template<INT N> __device__ __inline__ deviceDN<N> _expreal(const deviceDN<N>& x, Real a) { return x.ExpReal(a); }
template<INT N> __device__ __inline__ void _norm(deviceDN<N>& x) { x.Norm(); }

// --- element access ---

template<INT N> __device__ __host__ __inline__ Real _element(const deviceDN<N>& x, INT idx)
{
    switch (idx)
    {
    case 0: return x.m_me[0].x;
    case 1: return x.m_me[0].y;
    case 2: return x.m_me[1].x;
    case 3: return x.m_me[1].y;
    case 4: return x.m_me[2].x;
    case 5: return x.m_me[2].y;
    case 6: return x.m_me[3].x;
    case 7: return x.m_me[3].y;
    default: return F(0.0);
    }
}

template<INT N> __device__ __host__ __inline__ void _setelement(deviceDN<N>& x, INT idx, Real v)
{
    switch (idx)
    {
    case 0: x.m_me[0].x = v; break;
    case 1: x.m_me[0].y = v; break;
    case 2: x.m_me[1].x = v; break;
    case 3: x.m_me[1].y = v; break;
    case 4: x.m_me[2].x = v; break;
    case 5: x.m_me[2].y = v; break;
    case 6: x.m_me[3].x = v; break;
    case 7: x.m_me[3].y = v; break;
    default: break;
    }
}

// --- print ---

template<INT N> __device__ __inline__ void _print(const deviceDN<N>& x) { x.DebugPrint(); }
template<INT N> __device__ __inline__ void _print(const deviceDN<N>& x, const char* head) { x.DebugPrint(head); }

// --- detv / absv ---

template<INT N> __device__ __inline__ CLGComplex _detv(const deviceDN<N>& x) { return x.Det(); }
template<INT N> __device__ __inline__ Real _absv(const deviceDN<N>& x) { return _cuCabsf(x.Det()); }

__END_NAMESPACE

#include "DeviceTemplates/DeviceInlineGauge.h"
#include "DeviceTemplates/DeviceInlineGaugeChair.h"
#include "DeviceTemplates/DeviceInlineStaggeredGamma.h"
#include "DeviceTemplates/DeviceInlineStaggeredRotation.h"


#endif //#ifndef _DEVICEINLINETEMPLATE_H_

//=============================================================================
// END OF FILE
//=============================================================================
