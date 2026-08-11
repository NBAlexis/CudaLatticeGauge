//=============================================================================
// FILENAME : DeviceInlineTemplate.h
// 
// DESCRIPTION:
// This is a test, to see whether we can speed up by reduce the usage of registers
// template<typename TMatrix, typename TVector> __device__ __inline__ void _mulVecP(TVector* res, const TMatrix* __restrict__ matrix, const TVector* __restrict__ vector) = delete;
// This doesn't work
// Why?
//
// REVISION:
//  [mm/dd/yy]
//  [01/27/2025 nbale]
//=============================================================================

#ifndef _DEVICEINLINETEMPLATE2_H_
#define _DEVICEINLINETEMPLATE2_H_



__BEGIN_NAMESPACE

//It turns out these functions cannot reduce the usage of registers
//Maybe ptx has already optimized these functions

#if _CLG_SAVE_REGISTOR_STRATEGY

#pragma region add

#pragma endregion

#pragma region sub

#pragma endregion

#pragma region multiply

template<typename T1, typename T2> __device__ __inline__ void _mulP(T1* a, const T2* b) = delete;
template<typename T1, typename T2, typename T3> __device__ __inline__ void _mulP(T1* a, const T2* b, const T3* c) = delete;

//=================== mm ================================= it seems we can directly optimize the Mul function

//this complex is useless, just for complex-valued field template
//template<> __device__ __inline__ void _mulP<CLGComplex, CLGComplex>(CLGComplex* a, const CLGComplex* c)
//{
//    *a = _cuCmulf(*a, *c);
//}
//
//template<> __device__ __inline__ void _mulP<CLGComplex, CLGComplex, CLGComplex>(CLGComplex* a, const CLGComplex* b, const CLGComplex* c)
//{
//    *a = _cuCmulf(*b, *c);
//}
//
//template<> __device__ __inline__ void _mulP<deviceSU2, deviceSU2>(deviceSU2* a, const deviceSU2* c)
//{
//    CLGComplex tmp1 = _cuCmulf(a->m_me[0], c->m_me[0]);
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[1], c->m_me[2]));
//
//    CLGComplex tmp2 = _cuCmulf(a->m_me[0], c->m_me[1]);
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[1], c->m_me[3]));
//
//    a->m_me[0] = tmp1;
//    a->m_me[1] = tmp2;
//
//    tmp1 = _cuCmulf(a->m_me[2], c->m_me[0]);
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[3], c->m_me[2]));
//
//    tmp2 = _cuCmulf(a->m_me[2], c->m_me[1]);
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[3], c->m_me[3]));
//
//    a->m_me[2] = tmp1;
//    a->m_me[3] = tmp2;
//}
//
//template<> __device__ __inline__ void _mulP<deviceSU3, deviceSU3>(deviceSU3* a, const deviceSU3* c)
//{
//    CLGComplex tmp1 = _cuCmulf(a->m_me[0], c->m_me[0]);
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[1], c->m_me[3]));
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[2], c->m_me[6]));
//
//    CLGComplex tmp2 = _cuCmulf(a->m_me[0], c->m_me[1]);
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[1], c->m_me[4]));
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[2], c->m_me[7]));
//
//    CLGComplex tmp3 = _cuCmulf(a->m_me[0], c->m_me[2]);
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[1], c->m_me[5]));
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[2], c->m_me[7]));
//
//    a->m_me[0] = tmp1;
//    a->m_me[1] = tmp2;
//    a->m_me[2] = tmp3;
//
//    tmp1 = _cuCmulf(a->m_me[3], c->m_me[0]);
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[4], c->m_me[3]));
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[5], c->m_me[6]));
//
//    tmp2 = _cuCmulf(a->m_me[3], c->m_me[1]);
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[4], c->m_me[4]));
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[5], c->m_me[7]));
//
//    tmp3 = _cuCmulf(a->m_me[3], c->m_me[2]);
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[4], c->m_me[5]));
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[5], c->m_me[7]));
//
//    a->m_me[3] = tmp1;
//    a->m_me[4] = tmp2;
//    a->m_me[5] = tmp3;
//
//    tmp1 = _cuCmulf(a->m_me[6], c->m_me[0]);
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[7], c->m_me[3]));
//    tmp1 = _cuCaddf(tmp1, _cuCmulf(a->m_me[8], c->m_me[6]));
//
//    tmp2 = _cuCmulf(a->m_me[6], c->m_me[1]);
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[7], c->m_me[4]));
//    tmp2 = _cuCaddf(tmp2, _cuCmulf(a->m_me[8], c->m_me[7]));
//
//    tmp3 = _cuCmulf(a->m_me[6], c->m_me[2]);
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[7], c->m_me[5]));
//    tmp3 = _cuCaddf(tmp3, _cuCmulf(a->m_me[8], c->m_me[7]));
//
//    a->m_me[6] = tmp1;
//    a->m_me[7] = tmp2;
//    a->m_me[8] = tmp3;
//}

//=================== mv =================================

template<> __device__ __inline__ void _mulP<deviceSU2Vector, deviceSU2, deviceSU2Vector>(deviceSU2Vector* res, const deviceSU2* matrix, const deviceSU2Vector* vector)
{
    res->m_ve[0] = _cuCmulf(matrix->m_me[0], vector->m_ve[0]);
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));

    res->m_ve[1] = _cuCmulf(matrix->m_me[2], vector->m_ve[0]);
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[3], vector->m_ve[1]));
}

template<> __device__ __inline__ void _mulP<deviceSU3Vector, deviceSU3, deviceSU3Vector>(deviceSU3Vector* res, const deviceSU3* matrix, const deviceSU3Vector* vector)
{
    res->m_ve[0] = _cuCmulf(matrix->m_me[0], vector->m_ve[0]);
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[2], vector->m_ve[2]));

    res->m_ve[1] = _cuCmulf(matrix->m_me[3], vector->m_ve[0]);
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[4], vector->m_ve[1]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[5], vector->m_ve[2]));

    res->m_ve[2] = _cuCmulf(matrix->m_me[6], vector->m_ve[0]);
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(matrix->m_me[7], vector->m_ve[1]));
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(matrix->m_me[8], vector->m_ve[2]));
}

template<> __device__ __inline__ void _mulP<deviceWilsonVectorSU3, deviceSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3* res, const deviceSU3* matrix, const deviceWilsonVectorSU3* vector)
{
    _mulP(res->m_d, matrix, vector->m_d);
    _mulP(res->m_d + 1, matrix, vector->m_d + 1);
    _mulP(res->m_d + 2, matrix, vector->m_d + 2);
    _mulP(res->m_d + 3, matrix, vector->m_d + 3);
}

template<INT n, INT noeM, INT noeV> __device__ __inline__ void _mulP(
    deviceSUNVector<n, noeV>* res, 
    const deviceSUN<n, noeM>* __restrict__ matrix, 
    const deviceSUNVector<n, noeV>* __restrict__ vector)
{
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (BYTE x = 0U; x < n; ++x)
    {
        res->m_ve[x] = _cuCmulf(matrix->m_me[x * n], vector->m_ve[0]);
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (BYTE y = 1U; y < n; ++y)
        {
            res->m_ve[x] = _cuCaddf(res->m_ve[x], _cuCmulf(matrix->m_me[x * n + y], vector->m_ve[y]));
        }
    }
}


#pragma endregion

#pragma region dagger-multiply



#pragma endregion

#pragma region add-multiply

template<typename T1, typename T2, typename T3> __device__ __inline__ void _addmulP(T1* a, const T2* b, const T3* c) = delete;

template<> __device__ __inline__ void _addmulP<Real, Real, Real>(Real* res, const Real* matrix, const Real* vector)
{
    *res += (*matrix) * (*vector);
}

template<> __device__ __inline__ void _addmulP<CLGComplex, CLGComplex, CLGComplex>(CLGComplex* res, const CLGComplex* matrix, const CLGComplex* vector)
{
    *res = _cuCaddf(*res, _cuCmulf(*matrix, *vector));
}

template<> __device__ __inline__ void _addmulP<deviceSU2Vector, deviceSU2, deviceSU2Vector>(deviceSU2Vector* res, const deviceSU2* matrix, const deviceSU2Vector* vector)
{
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[0], vector->m_ve[0]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));

    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[2], vector->m_ve[0]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[3], vector->m_ve[1]));
}

template<> __device__ __inline__ void _addmulP<deviceSU3Vector, deviceSU3, deviceSU3Vector>(deviceSU3Vector* res, const deviceSU3* matrix, const deviceSU3Vector* vector)
{
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[0], vector->m_ve[0]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(matrix->m_me[2], vector->m_ve[2]));

    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[3], vector->m_ve[0]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[4], vector->m_ve[1]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(matrix->m_me[5], vector->m_ve[2]));

    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(matrix->m_me[6], vector->m_ve[0]));
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(matrix->m_me[7], vector->m_ve[1]));
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(matrix->m_me[8], vector->m_ve[2]));
}

template<> __device__ __inline__ void _addmulP<deviceWilsonVectorSU3, deviceSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3* res, const deviceSU3* matrix, const deviceWilsonVectorSU3* vector)
{
    _addmulP(res->m_d, matrix, vector->m_d);
    _addmulP(res->m_d + 1, matrix, vector->m_d + 1);
    _addmulP(res->m_d + 2, matrix, vector->m_d + 2);
    _addmulP(res->m_d + 3, matrix, vector->m_d + 3);
}

template<INT n, INT noeM, INT noeV> __device__ __inline__ void _addmulP(
    deviceSUNVector<n, noeV>* res,
    const deviceSUN<n, noeM>* __restrict__ matrix,
    const deviceSUNVector<n, noeV>* __restrict__ vector)
{
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (BYTE x = 0U; x < n; ++x)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (BYTE y = 0U; y < n; ++y)
        {
            res->m_ve[x] = _cuCaddf(res->m_ve[x], _cuCmulf(matrix->m_me[x * n + y], vector->m_ve[y]));
        }
    }
}

#pragma endregion

#pragma region sub-multiply

template<typename T1, typename T2, typename T3> __device__ __inline__ void _submulP(T1* a, const T2* b, const T3* c) = delete;

template<> __device__ __inline__ void _submulP<Real, Real, Real>(Real* res, const Real* matrix, const Real* vector)
{
    *res -= (*matrix) * (*vector);
}

template<> __device__ __inline__ void _submulP<CLGComplex, CLGComplex, CLGComplex>(CLGComplex* res, const CLGComplex* matrix, const CLGComplex* vector)
{
    *res = _cuCsubf(*res, _cuCmulf(*matrix, *vector));
}

template<> __device__ __inline__ void _submulP<deviceSU2Vector, deviceSU2, deviceSU2Vector>(deviceSU2Vector* res, const deviceSU2* matrix, const deviceSU2Vector* vector)
{
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(matrix->m_me[0], vector->m_ve[0]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));

    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(matrix->m_me[2], vector->m_ve[0]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(matrix->m_me[3], vector->m_ve[1]));
}

template<> __device__ __inline__ void _submulP<deviceSU3Vector, deviceSU3, deviceSU3Vector>(deviceSU3Vector* res, const deviceSU3* matrix, const deviceSU3Vector* vector)
{
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(matrix->m_me[0], vector->m_ve[0]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(matrix->m_me[1], vector->m_ve[1]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(matrix->m_me[2], vector->m_ve[2]));

    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(matrix->m_me[3], vector->m_ve[0]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(matrix->m_me[4], vector->m_ve[1]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(matrix->m_me[5], vector->m_ve[2]));

    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(matrix->m_me[6], vector->m_ve[0]));
    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(matrix->m_me[7], vector->m_ve[1]));
    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(matrix->m_me[8], vector->m_ve[2]));
}

template<> __device__ __inline__ void _submulP<deviceWilsonVectorSU3, deviceSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3* res, const deviceSU3* matrix, const deviceWilsonVectorSU3* vector)
{
    _submulP(res->m_d, matrix, vector->m_d);
    _submulP(res->m_d + 1, matrix, vector->m_d + 1);
    _submulP(res->m_d + 2, matrix, vector->m_d + 2);
    _submulP(res->m_d + 3, matrix, vector->m_d + 3);
}

template<INT n, INT noeM, INT noeV> __device__ __inline__ void _submulP(
    deviceSUNVector<n, noeV>* res,
    const deviceSUN<n, noeM>* __restrict__ matrix,
    const deviceSUNVector<n, noeV>* __restrict__ vector)
{
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (BYTE x = 0U; x < n; ++x)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (BYTE y = 0U; y < n; ++y)
        {
            res->m_ve[x] = _cuCsubf(res->m_ve[x], _cuCmulf(matrix->m_me[x * n + y], vector->m_ve[y]));
        }
    }
}

#pragma endregion

#pragma region add-dagger-multiply

template<typename T1, typename T2, typename T3> __device__ __inline__ void _adddagmulP(T1* a, const T2* b, const T3* c) = delete;

template<> __device__ __inline__ void _adddagmulP<Real, Real, Real>(Real* res, const Real* matrix, const Real* vector)
{
    *res += (*matrix) * (*vector);
}

template<> __device__ __inline__ void _adddagmulP<CLGComplex, CLGComplex, CLGComplex>(CLGComplex* res, const CLGComplex* matrix, const CLGComplex* vector)
{
    *res = _cuCaddf(_cuConjf(*res), _cuCmulf(*matrix, *vector));
}

template<> __device__ __inline__ void _adddagmulP<deviceSU2Vector, deviceSU2, deviceSU2Vector>(deviceSU2Vector* res, const deviceSU2* matrix, const deviceSU2Vector* vector)
{
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[0]), vector->m_ve[0]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[2]), vector->m_ve[1]));

    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[1]), vector->m_ve[0]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[3]), vector->m_ve[1]));
}

template<> __device__ __inline__ void _adddagmulP<deviceSU3Vector, deviceSU3, deviceSU3Vector>(deviceSU3Vector* res, const deviceSU3* matrix, const deviceSU3Vector* vector)
{
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[0]), vector->m_ve[0]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[3]), vector->m_ve[1]));
    res->m_ve[0] = _cuCaddf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[6]), vector->m_ve[2]));

    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[1]), vector->m_ve[0]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[4]), vector->m_ve[1]));
    res->m_ve[1] = _cuCaddf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[7]), vector->m_ve[2]));

    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[2]), vector->m_ve[0]));
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[5]), vector->m_ve[1]));
    res->m_ve[2] = _cuCaddf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[8]), vector->m_ve[2]));
}

template<> __device__ __inline__ void _adddagmulP<deviceWilsonVectorSU3, deviceSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3* res, const deviceSU3* matrix, const deviceWilsonVectorSU3* vector)
{
    _adddagmulP(res->m_d, matrix, vector->m_d);
    _adddagmulP(res->m_d + 1, matrix, vector->m_d + 1);
    _adddagmulP(res->m_d + 2, matrix, vector->m_d + 2);
    _adddagmulP(res->m_d + 3, matrix, vector->m_d + 3);
}

template<INT n, INT noeM, INT noeV> __device__ __inline__ void _adddagmulP(
    deviceSUNVector<n, noeV>* res,
    const deviceSUN<n, noeM>* __restrict__ matrix,
    const deviceSUNVector<n, noeV>* __restrict__ vector)
{
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (BYTE x = 0U; x < n; ++x)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (BYTE y = 0U; y < n; ++y)
        {
            res->m_ve[x] = _cuCaddf(res->m_ve[x], _cuCmulf(_cuConjf(matrix->m_me[y * n + x]), vector->m_ve[y]));
        }
    }
}

#pragma endregion

#pragma region sub-sub-multiply

template<typename T1, typename T2, typename T3> __device__ __inline__ void _subdagmulP(T1* a, const T2* b, const T3* c) = delete;

template<> __device__ __inline__ void _subdagmulP<Real, Real, Real>(Real* res, const Real* matrix, const Real* vector)
{
    *res -= (*matrix) * (*vector);
}

template<> __device__ __inline__ void _subdagmulP<CLGComplex, CLGComplex, CLGComplex>(CLGComplex* res, const CLGComplex* matrix, const CLGComplex* vector)
{
    *res = _cuCsubf(_cuConjf(*res), _cuCmulf(*matrix, *vector));
}

template<> __device__ __inline__ void _subdagmulP<deviceSU2Vector, deviceSU2, deviceSU2Vector>(deviceSU2Vector* res, const deviceSU2* matrix, const deviceSU2Vector* vector)
{
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[0]), vector->m_ve[0]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[2]), vector->m_ve[1]));

    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[1]), vector->m_ve[0]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[3]), vector->m_ve[1]));
}

template<> __device__ __inline__ void _subdagmulP<deviceSU3Vector, deviceSU3, deviceSU3Vector>(deviceSU3Vector* res, const deviceSU3* matrix, const deviceSU3Vector* vector)
{
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[0]), vector->m_ve[0]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[3]), vector->m_ve[1]));
    res->m_ve[0] = _cuCsubf(res->m_ve[0], _cuCmulf(_cuConjf(matrix->m_me[6]), vector->m_ve[2]));

    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[1]), vector->m_ve[0]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[4]), vector->m_ve[1]));
    res->m_ve[1] = _cuCsubf(res->m_ve[1], _cuCmulf(_cuConjf(matrix->m_me[7]), vector->m_ve[2]));

    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[2]), vector->m_ve[0]));
    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[5]), vector->m_ve[1]));
    res->m_ve[2] = _cuCsubf(res->m_ve[2], _cuCmulf(_cuConjf(matrix->m_me[8]), vector->m_ve[2]));
}

template<> __device__ __inline__ void _subdagmulP<deviceWilsonVectorSU3, deviceSU3, deviceWilsonVectorSU3>(deviceWilsonVectorSU3* res, const deviceSU3* matrix, const deviceWilsonVectorSU3* vector)
{
    _subdagmulP(res->m_d, matrix, vector->m_d);
    _subdagmulP(res->m_d + 1, matrix, vector->m_d + 1);
    _subdagmulP(res->m_d + 2, matrix, vector->m_d + 2);
    _subdagmulP(res->m_d + 3, matrix, vector->m_d + 3);
}

template<INT n, INT noeM, INT noeV> __device__ __inline__ void _subdagmulP(
    deviceSUNVector<n, noeV>* res,
    const deviceSUN<n, noeM>* __restrict__ matrix,
    const deviceSUNVector<n, noeV>* __restrict__ vector)
{
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
    for (BYTE x = 0U; x < n; ++x)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (BYTE y = 0U; y < n; ++y)
        {
            res->m_ve[x] = _cuCsubf(res->m_ve[x], _cuCmulf(_cuConjf(matrix->m_me[y * n + x]), vector->m_ve[y]));
        }
    }
}

#pragma endregion

#endif

__END_NAMESPACE
//#endif //__CUDACC__

#endif //#ifndef _DEVICEINLINETEMPLATE2_H_

//=============================================================================
// END OF FILE
//=============================================================================
