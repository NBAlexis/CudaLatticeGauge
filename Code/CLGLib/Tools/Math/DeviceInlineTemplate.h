//=============================================================================
// FILENAME : DeviceInlineTemplate.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

#ifndef _DEVICEINLINETEMPLATE_H_
#define _DEVICEINLINETEMPLATE_H_

__BEGIN_NAMESPACE

#pragma region element-wise

template<typename T> __device__ __inline__ T _makeId() = delete;

template<> __device__ __inline__ Real _makeId<Real>()
{
    return F(1.0);
}

template<> __device__ __inline__ CLGComplex _makeId<CLGComplex>()
{
    return _onec;
}

template<> __device__ __inline__ deviceSU2 _makeId<deviceSU2>()
{
    return deviceSU2::makeSU2Id();
}

template<> __device__ __inline__ deviceSU3 _makeId<deviceSU3>()
{
    return deviceSU3::makeSU3Id();
}

template<> __device__ __inline__ deviceSU4 _makeId<deviceSU4>()
{
    return deviceSU4::makeSUNId();
}

template<> __device__ __inline__ deviceSU5 _makeId<deviceSU5>()
{
    return deviceSU5::makeSUNId();
}

template<> __device__ __inline__ deviceSU6 _makeId<deviceSU6>()
{
    return deviceSU6::makeSUNId();
}

template<> __device__ __inline__ deviceSU7 _makeId<deviceSU7>()
{
    return deviceSU7::makeSUNId();
}

template<> __device__ __inline__ deviceSU8 _makeId<deviceSU8>()
{
    return deviceSU8::makeSUNId();
}

template<> __device__ __inline__ deviceSU2Vector _makeId<deviceSU2Vector>()
{
    return deviceSU2Vector::makeOneSU2Vector();
}

template<> __device__ __inline__ deviceSU3Vector _makeId<deviceSU3Vector>()
{
    return deviceSU3Vector::makeOneSU3Vector();
}

template<> __device__ __inline__ deviceSU4Vector _makeId<deviceSU4Vector>()
{
    return deviceSU4Vector::makeOneSUNVector();
}

template<> __device__ __inline__ deviceSU5Vector _makeId<deviceSU5Vector>()
{
    return deviceSU5Vector::makeOneSUNVector();
}

template<> __device__ __inline__ deviceSU6Vector _makeId<deviceSU6Vector>()
{
    return deviceSU6Vector::makeOneSUNVector();
}

template<> __device__ __inline__ deviceSU7Vector _makeId<deviceSU7Vector>()
{
    return deviceSU7Vector::makeOneSUNVector();
}

template<> __device__ __inline__ deviceSU8Vector _makeId<deviceSU8Vector>()
{
    return deviceSU8Vector::makeOneSUNVector();
}



template<typename T> __device__ __inline__ T _makeZero() = delete;

template<> __device__ __inline__ Real _makeZero<Real>()
{
    return F(0.0);
}

template<> __device__ __inline__ CLGComplex _makeZero<CLGComplex>()
{
    return _zeroc;
}

template<> __device__ __inline__ deviceSU2 _makeZero<deviceSU2>()
{
    return deviceSU2::makeSU2Zero();
}

template<> __device__ __inline__ deviceSU3 _makeZero<deviceSU3>()
{
    return deviceSU3::makeSU3Zero();
}

template<> __device__ __inline__ deviceSU4 _makeZero<deviceSU4>()
{
    return deviceSU4::makeSUNZero();
}

template<> __device__ __inline__ deviceSU5 _makeZero<deviceSU5>()
{
    return deviceSU5::makeSUNZero();
}

template<> __device__ __inline__ deviceSU6 _makeZero<deviceSU6>()
{
    return deviceSU6::makeSUNZero();
}

template<> __device__ __inline__ deviceSU7 _makeZero<deviceSU7>()
{
    return deviceSU7::makeSUNZero();
}

template<> __device__ __inline__ deviceSU8 _makeZero<deviceSU8>()
{
    return deviceSU8::makeSUNZero();
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

template<> __device__ __inline__ deviceSU4Vector _makeZero<deviceSU4Vector>()
{
    return deviceSU4Vector::makeZeroSUNVector();
}

template<> __device__ __inline__ deviceSU5Vector _makeZero<deviceSU5Vector>()
{
    return deviceSU5Vector::makeZeroSUNVector();
}

template<> __device__ __inline__ deviceSU6Vector _makeZero<deviceSU6Vector>()
{
    return deviceSU6Vector::makeZeroSUNVector();
}

template<> __device__ __inline__ deviceSU7Vector _makeZero<deviceSU7Vector>()
{
    return deviceSU7Vector::makeZeroSUNVector();
}

template<> __device__ __inline__ deviceSU8Vector _makeZero<deviceSU8Vector>()
{
    return deviceSU8Vector::makeZeroSUNVector();
}

template<typename TMatrix, typename TVector> __device__ __inline__ TMatrix _makeContract(const TVector& left, const TVector& right) = delete;

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

template<> __device__ __inline__ deviceSU4 _makeContract<deviceSU4, deviceSU4Vector>(const deviceSU4Vector& left, const deviceSU4Vector& right)
{
    return deviceSU4::makeSUNContractV(left, right);
}

template<> __device__ __inline__ deviceSU5 _makeContract<deviceSU5, deviceSU5Vector>(const deviceSU5Vector& left, const deviceSU5Vector& right)
{
    return deviceSU5::makeSUNContractV(left, right);
}

template<> __device__ __inline__ deviceSU6 _makeContract<deviceSU6, deviceSU6Vector>(const deviceSU6Vector& left, const deviceSU6Vector& right)
{
    return deviceSU6::makeSUNContractV(left, right);
}

template<> __device__ __inline__ deviceSU7 _makeContract<deviceSU7, deviceSU7Vector>(const deviceSU7Vector& left, const deviceSU7Vector& right)
{
    return deviceSU7::makeSUNContractV(left, right);
}

template<> __device__ __inline__ deviceSU8 _makeContract<deviceSU8, deviceSU8Vector>(const deviceSU8Vector& left, const deviceSU8Vector& right)
{
    return deviceSU8::makeSUNContractV(left, right);
}

//This is to make white noise, so for gauge, it is random generator
template<typename T> __device__ __inline__ T _makeGaussian(UINT fatIdx) = delete;

template<> __device__ __inline__ Real _makeGaussian<Real>(UINT fatIdx)
{
    return _deviceRandomGaussFSqrt2(fatIdx);
}

template<> __device__ __inline__ CLGComplex _makeGaussian<CLGComplex>(UINT fatIdx)
{
    return _deviceRandomGaussC(fatIdx);
}

template<> __device__ __inline__ deviceSU2Vector _makeGaussian<deviceSU2Vector>(UINT fatIdx)
{
    return deviceSU2Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU3Vector _makeGaussian<deviceSU3Vector>(UINT fatIdx)
{
    return deviceSU3Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU4Vector _makeGaussian<deviceSU4Vector>(UINT fatIdx)
{
    return deviceSU4Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU5Vector _makeGaussian<deviceSU5Vector>(UINT fatIdx)
{
    return deviceSU5Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU6Vector _makeGaussian<deviceSU6Vector>(UINT fatIdx)
{
    return deviceSU6Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU7Vector _makeGaussian<deviceSU7Vector>(UINT fatIdx)
{
    return deviceSU7Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU8Vector _makeGaussian<deviceSU8Vector>(UINT fatIdx)
{
    return deviceSU8Vector::makeRandomGaussian(fatIdx);
}

template<> __device__ __inline__ deviceSU2 _makeGaussian<deviceSU2>(UINT fatIdx)
{
    return deviceSU2::makeSU2RandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU3 _makeGaussian<deviceSU3>(UINT fatIdx)
{
    return deviceSU3::makeSU3RandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU4 _makeGaussian<deviceSU4>(UINT fatIdx)
{
    return deviceSU4::makeSUNRandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU5 _makeGaussian<deviceSU5>(UINT fatIdx)
{
    return deviceSU5::makeSUNRandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU6 _makeGaussian<deviceSU6>(UINT fatIdx)
{
    return deviceSU6::makeSUNRandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU7 _makeGaussian<deviceSU7>(UINT fatIdx)
{
    return deviceSU7::makeSUNRandomGenerator(fatIdx);
}

template<> __device__ __inline__ deviceSU8 _makeGaussian<deviceSU8>(UINT fatIdx)
{
    return deviceSU8::makeSUNRandomGenerator(fatIdx);
}

template<typename T> __device__ __inline__ T _makeZ4(UINT fatIdx) = delete;

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

template<> __device__ __inline__ deviceSU4Vector _makeZ4<deviceSU4Vector>(UINT fatIdx)
{
    return deviceSU4Vector::makeRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU5Vector _makeZ4<deviceSU5Vector>(UINT fatIdx)
{
    return deviceSU5Vector::makeRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU6Vector _makeZ4<deviceSU6Vector>(UINT fatIdx)
{
    return deviceSU6Vector::makeRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU7Vector _makeZ4<deviceSU7Vector>(UINT fatIdx)
{
    return deviceSU7Vector::makeRandomZ4(fatIdx);
}

template<> __device__ __inline__ deviceSU8Vector _makeZ4<deviceSU8Vector>(UINT fatIdx)
{
    return deviceSU8Vector::makeRandomZ4(fatIdx);
}


template<typename T> __device__ __inline__ T _makeRandom(UINT fatIdx) = delete;

template<> __device__ __inline__ Real _makeRandom<Real>(UINT fatIdx)
{
    return _deviceRandomF(fatIdx);
}

template<> __device__ __inline__ CLGComplex _makeRandom<CLGComplex>(UINT fatIdx)
{
    return _deviceRandomGaussC(fatIdx);
}

template<> __device__ __inline__ deviceSU2 _makeRandom<deviceSU2>(UINT fatIdx)
{
    return deviceSU2::makeSU2Random(fatIdx);
}

template<> __device__ __inline__ deviceSU3 _makeRandom<deviceSU3>(UINT fatIdx)
{
    return deviceSU3::makeSU3Random(fatIdx);
}

template<> __device__ __inline__ deviceSU4 _makeRandom<deviceSU4>(UINT fatIdx)
{
    return deviceSU4::makeSUNRandom(fatIdx);
}

template<> __device__ __inline__ deviceSU5 _makeRandom<deviceSU5>(UINT fatIdx)
{
    return deviceSU5::makeSUNRandom(fatIdx);
}

template<> __device__ __inline__ deviceSU6 _makeRandom<deviceSU6>(UINT fatIdx)
{
    return deviceSU6::makeSUNRandom(fatIdx);
}

template<> __device__ __inline__ deviceSU7 _makeRandom<deviceSU7>(UINT fatIdx)
{
    return deviceSU7::makeSUNRandom(fatIdx);
}

template<> __device__ __inline__ deviceSU8 _makeRandom<deviceSU8>(UINT fatIdx)
{
    return deviceSU8::makeSUNRandom(fatIdx);
}

template<typename T> __device__ __inline__ void _dagger(T& element) = delete;
template<typename T> __device__ __inline__ T _daggerC(const T& element) = delete;

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

template<typename TLeft, typename TRight> __device__ __inline__ TLeft _addC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _add(TLeft& left, const TRight& right) = delete;


template<> __device__ __inline__ Real _addC<Real, Real>(const Real& left, const Real& right)
{
    return left + right;
}
template<> __device__ __inline__ void _add<Real, Real>(Real& left, const Real& right)
{
    left = left + right;
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

#define __DEFINE_TWO_ELEMENT_Func(TYPENAME, FUNC1, FUNC2) \
template<> __device__ __inline__ TYPENAME FUNC1##C<TYPENAME, Real>(const TYPENAME& left, const Real& right) \
{ \
    return left.FUNC2##RealC(right); \
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
template<> __device__ __inline__ void FUNC1<TYPENAME, CLGComplex>(TYPENAME& left, const CLGComplex& right) \
{ \
    return left.FUNC2##Comp(right); \
} \
template<> __device__ __inline__ void FUNC1<TYPENAME, TYPENAME>(TYPENAME& left, const TYPENAME& right) \
{ \
    return left.FUNC2(right); \
} \


__DEFINE_TWO_ELEMENT_Func(deviceSU2Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU3Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU4Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU5Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU6Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU7Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU8Vector, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU4, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU5, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU6, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU7, _add, Add)
__DEFINE_TWO_ELEMENT_Func(deviceSU8, _add, Add)


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
__DEFINE_TWO_ELEMENT_Func(deviceSU4Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU5Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU6Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU7Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU8Vector, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU4, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU5, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU6, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU7, _sub, Sub)
__DEFINE_TWO_ELEMENT_Func(deviceSU8, _sub, Sub)


template<typename TLeft, typename TRight> __device__ __inline__ TLeft _mulC(const TLeft& left, const TRight& right) = delete;
template<typename TLeft, typename TRight> __device__ __inline__ void _mul(TLeft& left, const TRight& right) = delete;


template<> __device__ __inline__ Real _mulC<Real, Real>(const Real& left, const Real& right)
{
    return left * right;
}
template<> __device__ __inline__ void _mul<Real, Real>(Real& left, const Real& right)
{
    left = left * right;
}
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, Real>(const CLGComplex& left, const Real& right)
{
    return cuCmulf_cr(left, right);
}
template<> __device__ __inline__ void _mul<CLGComplex, Real>(CLGComplex& left, const Real& right)
{
    left = cuCmulf_cr(left, right);
}
template<> __device__ __inline__ CLGComplex _mulC<CLGComplex, CLGComplex>(const CLGComplex& left, const CLGComplex& right)
{
    return _cuCmulf(left, right);
}
template<> __device__ __inline__ void _mul<CLGComplex, CLGComplex>(CLGComplex& left, const CLGComplex& right)
{
    left = _cuCmulf(left, right);
}


__DEFINE_TWO_ELEMENT_Func(deviceSU2Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU3Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU4Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU5Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU6Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU7Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU8Vector, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU2, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU3, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU4, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU5, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU6, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU7, _mul, Mul)
__DEFINE_TWO_ELEMENT_Func(deviceSU8, _mul, Mul)

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
__DEFINE_TWO_ELEMENT_Func2(deviceSU4Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU5Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU6Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU7Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU8Vector, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU2, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU4, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU5, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU6, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU7, _dagmul, DaggerMul)
__DEFINE_TWO_ELEMENT_Func2(deviceSU8, _dagmul, DaggerMul)


template<typename T> __device__ __inline__ T _muldagC(const T& left, const T& right) = delete;
template<typename T> __device__ __inline__ void _muldag(T& left, const T& right) = delete;
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
__DEFINE_TWO_ELEMENT_Func2(deviceSU4Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU5Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU6Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU7Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU8Vector, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU2, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU3, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU4, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU5, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU6, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU7, _muldag, MulDagger)
__DEFINE_TWO_ELEMENT_Func2(deviceSU8, _muldag, MulDagger)

template<typename TMatrix, typename TVector> __device__ __inline__ TVector _mulVec(const TMatrix& matrix, const TVector& vector) = delete;

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

template<> __device__ __inline__ deviceSU4Vector _mulVec<deviceSU4, deviceSU4Vector>(const deviceSU4& matrix, const deviceSU4Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceSU5Vector _mulVec<deviceSU5, deviceSU5Vector>(const deviceSU5& matrix, const deviceSU5Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceSU6Vector _mulVec<deviceSU6, deviceSU6Vector>(const deviceSU6& matrix, const deviceSU6Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceSU7Vector _mulVec<deviceSU7, deviceSU7Vector>(const deviceSU7& matrix, const deviceSU7Vector& vector)
{
    return matrix.MulVector(vector);
}

template<> __device__ __inline__ deviceSU8Vector _mulVec<deviceSU8, deviceSU8Vector>(const deviceSU8& matrix, const deviceSU8Vector& vector)
{
    return matrix.MulVector(vector);
}

template<typename TMatrix> __device__ __inline__ void _ta(TMatrix& matrix) = delete;

template<> __device__ __inline__ void _ta<CLGComplex>(CLGComplex& matrix)
{
    matrix.x = F(0.0);
}

template<> __device__ __inline__ void _ta<deviceSU2>(deviceSU2& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU3>(deviceSU3& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU4>(deviceSU4& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU5>(deviceSU5& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU6>(deviceSU6& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU7>(deviceSU7& matrix)
{
    matrix.Ta();
}

template<> __device__ __inline__ void _ta<deviceSU8>(deviceSU8& matrix)
{
    matrix.Ta();
}

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

template<typename T> __device__ __host__ __inline__  BYTE _dim() = delete;

template<> __device__ __host__ __inline__  BYTE _dim<Real>() { return 1; }
template<> __device__ __host__ __inline__  BYTE _dim<CLGComplex>() { return 1; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU2Vector>() { return 2; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU3Vector>() { return 3; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU4Vector>() { return 4; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU5Vector>() { return 5; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU6Vector>() { return 6; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU7Vector>() { return 7; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU8Vector>() { return 8; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU2>() { return 2; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU3>() { return 3; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU4>() { return 4; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU5>() { return 5; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU6>() { return 6; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU7>() { return 7; }
template<> __device__ __host__ __inline__  BYTE _dim<deviceSU8>() { return 8; }

template<typename T> __device__ __inline__  T _expreal(const T& x, Real a) = delete;

template<> __device__ __inline__  Real _expreal<Real>(const Real& x, Real a) { return _exp(x * a); }
template<> __device__ __inline__  CLGComplex _expreal<CLGComplex>(const CLGComplex& x, Real a) { return __cuCexpf(cuCmulf_cr(x, a)); }
template<> __device__ __inline__  deviceSU2 _expreal<deviceSU2>(const deviceSU2& x, Real a) { return x.QuickExp(a); }

template<> __device__ __inline__  deviceSU3 _expreal<deviceSU3>(const deviceSU3& x, Real a)
{ 
    return (0 == _DC_ExpPrecision) ? x.QuickExp(a) : x.ExpReal(a, static_cast<BYTE>(_DC_ExpPrecision));
}

template<INT N, INT NoE> 
__device__ __inline__  deviceSUN<N, NoE> _expreal(const deviceSUN<N, NoE>& x, Real a)
{
    return x.ExpReal(a, _DC_ExpPrecision > N ? _DC_ExpPrecision : (N + 1));
}

template<typename T> __device__ __inline__  void _norm(T& x) = delete;
template<> __device__ __inline__  void _norm<deviceSU2>(deviceSU2& x) { x.Norm(); }
template<> __device__ __inline__  void _norm<deviceSU3>(deviceSU3& x) { x.Norm(); }
template<INT N, INT NoE> __device__ __inline__  void _norm(deviceSUN<N, NoE>& x) { x.Norm(); }

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

template<INT N, INT NoE>
__device__ __host__ __inline__ Real _element(const deviceSUN<N, NoE>& x, INT idx)
{
    if (idx < 2 * N * N)
    {
        const UINT idxc = (idx >> 1);
        if (idx & 1)
        {
            return x.m_ve[idxc].y;
        }
        return x.m_ve[idxc].x;
    }
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
        }
        x.m_me[idxc].x = v;
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
            x.m_ve[idxc].y = v;
        }
        x.m_ve[idxc].x = v;
    }
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
        v.m_me[1].x, v.m_me[1].y > F(0.0) ? _T("+") : _T("-"), appAbs(v.m_me[1].y)
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

#pragma endregion

#pragma region Gauge

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
template<typename deviceGauge>
static __device__ __inline__ const deviceGauge& _deviceGetGaugeBCT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryOne<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
template<typename deviceGauge>
static __device__ __inline__ const deviceGauge& _deviceGetGaugeBCDirT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        ((CFieldBoundaryOne<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ]
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirOneT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        _makeId<deviceGauge>()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirZeroT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        _makeZero<deviceGauge>()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    deviceGauge ret = idx.IsDirichlet() ?
        ((CFieldBoundaryOne<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        _dagger(ret);
    }
    return ret;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirOneSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return _makeId<deviceGauge>();
    }
    if (idx.NeedToDagger())
    {
        return _daggerC(pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)]);
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * Note that, when get zero instead of one, it is minus not dagger
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirZeroSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return _makeZero<deviceGauge>();
    }
    if (idx.NeedToDagger())
    {
        return _mulC(pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)], F(-1.0));
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * calculate D_mu A _nu = Delta _mu + [A_mu, A _nu]
 * Use U now to calculate A pure
 * me will be changed, so, if me is A phys, copy me first
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceDPureMuT(
    const deviceGauge* __restrict__ piA,
    const deviceGauge* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));

    deviceGauge res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceGauge res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    _mul(res2, res); //A _nu Apure _mu
    _mul(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    _sub(res, res2); //[Apure, A]
    _add(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu));
    _sub(res, _deviceGetGaugeBCDirZeroSIndexT(piA,__idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
    return res;
}

/**
 * test using (A(N+mu)-A(N-mu))/2
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceDPureMu2T(
    const deviceGauge* __restrict__ piA,
    const deviceGauge* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n+mu) - A_nu (n-mu))/2 + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));
    const UINT uiSiteBig_p_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, byMu + 1));

    deviceGauge res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceGauge res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    _mul(res2, res); //A _nu Apure _mu
    _mul(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    _sub(res, res2); //[Apure, A]
    _add(res, _mulC(_deviceGetGaugeBCDirZeroSIndexT(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]), F(0.5)));
    _sub(res, _mulC(_deviceGetGaugeBCDirZeroSIndexT(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]), F(0.5)));
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _devicePlaqutteT(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceGauge toAdd(_deviceGetGaugeBCDirOneSIndexT(pDeviceData, first));
    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceGauge toMul(_deviceGetGaugeBCDirOneSIndexT(pDeviceData, first));
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
        }
    }
    return toAdd;
}

/**
 * pDir[] is dirs of path, the dir is:
 *  x,y,z,t : 1,2,3,4
 *  -x,-y,-z,-t: -1,-2,-3,-4
 *
 * NOTE: This function assumes the boundary is always unity
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkT(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _dagger(sRet);
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

/**
 * After every move, it maps to inside the lattice
 * Do NOT use it in projective plane boundary condition
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkLongT(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];
        sStartSite = __deviceSiteIndexToInt4(newLink.m_uiSiteIndex);

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _dagger(sRet);
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

/**
* big index is the index of walking table.
* The plaqutte index may not be cached because n may out of boundary, so we calculate every one
* n, n+mu, n+nu, n
*
*   <----- ^
*   |      |
*   |      |
*   V      |
* O ------->
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermPPT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    //For any boundary condition it is always, site->mu, site_p_mu->nu, site_p_nu->mu+, site->nu+
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite4, byMu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite4, byNu + 1);
    const UINT uiB4 = uiBigIdx * _DC_Dir;
    const SIndex& s_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byMu];
    const SIndex& s_p_mu_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_mu) * _DC_Dir + byNu];
    const SIndex& s_p_nu_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_nu) * _DC_Dir + byMu];
    const SIndex& s_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byNu];

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_mu, byFieldId);
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_p_mu_nu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_p_nu_mu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_nu, byFieldId));

    return u;
}

/**
* U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
*
*    ------->
*    ^      |
*    |      |
*    |      V
*    <------- O
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermMPT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__mu, byFieldId);
    _dagmul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__nu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s__nu, byFieldId));

    return u;
}

/**
* U(mu,-nu) = U(N) U^+(N+mu-nu) U^+(N-nu) U(N-nu)
*
* O  ------->
*    ^      |
*    |      |
*    |      V
*    <-------
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermPMT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s__mu, byFieldId);
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__mu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}

/**
* U(-mu,-nu) = U^+(N-mu) U^+(N-mu-nu) U(N-mu-nu) U(N-nu)
*
* <----- ^ O
* |      |
* |      |
* V      |
* ------->
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermMMT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_m_mu = _deviceSmallInt4OffsetC(n_m_nu, -static_cast<INT>(byMu) - 1);
    const UINT uin_m_nu_m_mub4 = __idx->_deviceGetBigIndex(n_m_nu_m_mu) * _DC_Dir;

    const SIndex& s_m_nu_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byNu];
    const SIndex& s_m_nu_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byMu];
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir + byNu];

    //u1^+ u2^+ u3 u4
    //= (u2 u1)^+ u3 u4
    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_m_mu__nu, byFieldId);
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__mu, byFieldId));
    _dagmul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}

/**
 * U_{mu,nu}(n)+U^+_{-mu,nu}(n)+U^+_{mu,-nu}(n)+U_{-mu,-nu}(n)
 * or
 * U_{mu,nu}(n)+U_{nu,-mu}(n)+U_{-nu,mu}(n)+U_{-mu,-nu}(n) <--- we are using this one
 * or
 * U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu)
 *
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceCloverT(const deviceGauge* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceGauge ret = _device1PlaqutteTermPPT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId);
    _add(ret, _device1PlaqutteTermMMT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    _add(ret, _device1PlaqutteTermPMT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    _add(ret, _device1PlaqutteTermMPT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
template<typename deviceGauge>
static __device__ __inline__ Real _deviceCloverRetrT(const deviceGauge* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _retr(_device1PlaqutteTermPPT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermMMT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermPMT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermMPT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
}

#pragma endregion


__END_NAMESPACE

#endif //#ifndef _DEVICEINLINETEMPLATE_H_

//=============================================================================
// END OF FILE
//=============================================================================
