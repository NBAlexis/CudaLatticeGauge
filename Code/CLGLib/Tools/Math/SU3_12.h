//=============================================================================
// FILENAME : SU3_12.h
//
// DESCRIPTION:
// Compact SU(3) representation storing only the first 2 columns (6 complex = 12 reals).
// Third column is reconstructed via cross product: c3 = conj(c1 x c2).
// This assumes det=1, so it is only valid for SU(3) matrices.
//
// Operations that preserve SU(3) work directly on 6 elements.
// Operations that don't preserve SU(3) but are const return deviceSU3 (full 18 reals).
// Mutating operations that don't preserve SU(3) are not supported.
//
// REVISION:
//  [05/10/2026 nbale]
//=============================================================================

#ifndef _SU3_12_H_
#define _SU3_12_H_

#include "SU3.h"

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

struct deviceSU3_12
{
    // Storage: first 2 columns of a 3x3 complex matrix
    // col0: m_me[0], m_me[1], m_me[2]
    // col1: m_me[3], m_me[4], m_me[5]
    // col2: reconstructed via cross product
    CLGComplex m_me[6];

#pragma region helper

    // Reconstruct third column via cross product: c3 = conj(c1 x c2)
    // c1 = (m_me[0], m_me[1], m_me[2]), c2 = (m_me[3], m_me[4], m_me[5])
    // c1 x c2 = (c1[1]*c2[2] - c1[2]*c2[1], c1[2]*c2[0] - c1[0]*c2[2], c1[0]*c2[1] - c1[1]*c2[0])
    // c3 = conj(c1 x c2)
    __device__ __inline__ void _reconstructCol2(CLGComplex& c6, CLGComplex& c7, CLGComplex& c8) const
    {
        c6 = _cuConjf(_cuCsubf(_cuCmulf(m_me[1], m_me[5]), _cuCmulf(m_me[2], m_me[4])));
        c7 = _cuConjf(_cuCsubf(_cuCmulf(m_me[2], m_me[3]), _cuCmulf(m_me[0], m_me[5])));
        c8 = _cuConjf(_cuCsubf(_cuCmulf(m_me[0], m_me[4]), _cuCmulf(m_me[1], m_me[3])));
    }

#pragma endregion helper

#pragma region create

    __device__ deviceSU3_12()
    {
    }

    __device__ deviceSU3_12(const deviceSU3_12& other)
    {
        m_me[0] = other.m_me[0];
        m_me[1] = other.m_me[1];
        m_me[2] = other.m_me[2];
        m_me[3] = other.m_me[3];
        m_me[4] = other.m_me[4];
        m_me[5] = other.m_me[5];
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12Zero()
    {
        deviceSU3_12 ret;
        ret.m_me[0] = _zeroc;
        ret.m_me[1] = _zeroc;
        ret.m_me[2] = _zeroc;
        ret.m_me[3] = _zeroc;
        ret.m_me[4] = _zeroc;
        ret.m_me[5] = _zeroc;
        return ret;
    }

    __device__ __inline__ void Zero()
    {
        m_me[0] = _zeroc;
        m_me[1] = _zeroc;
        m_me[2] = _zeroc;
        m_me[3] = _zeroc;
        m_me[4] = _zeroc;
        m_me[5] = _zeroc;
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12Id()
    {
        deviceSU3_12 ret;
        ret.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
        ret.m_me[1] = _zeroc;
        ret.m_me[2] = _zeroc;
        ret.m_me[3] = _zeroc;
        ret.m_me[4] = _make_cuComplex(F(1.0), F(0.0));
        ret.m_me[5] = _zeroc;
        return ret;
    }

    __device__ __inline__ void Id()
    {
        m_me[0] = _make_cuComplex(F(1.0), F(0.0));
        m_me[1] = _zeroc;
        m_me[2] = _zeroc;
        m_me[3] = _zeroc;
        m_me[4] = _make_cuComplex(F(1.0), F(0.0));
        m_me[5] = _zeroc;
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12Random(UINT fatIndex)
    {
        deviceSU3 full = deviceSU3::makeSU3Random(fatIndex);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12RandomGenerator(UINT fatIndex)
    {
        deviceSU3 full = deviceSU3::makeSU3RandomGenerator(fatIndex);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12SumGenerator(Real fDivide)
    {
        deviceSU3 full = deviceSU3::makeSU3SumGenerator(fDivide);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12Generator(UINT uiGenerator)
    {
        deviceSU3 full = deviceSU3::makeSU3Generator(uiGenerator);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12TA(
        const CLGComplex& m12, const CLGComplex& m13, const CLGComplex& m23,
        Real m11, Real m22)
    {
        deviceSU3 full = deviceSU3::makeSU3TA(m12, m13, m23, m11, m22);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12Contract(
        const deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
    {
        deviceSU3 full = deviceSU3::makeSU3Contract(left, right);
        return deviceSU3_12(full);
    }

    __device__ __inline__ static deviceSU3_12 makeSU3_12ContractV(
        const deviceSU3Vector& left, const deviceSU3Vector& right)
    {
        deviceSU3 full = deviceSU3::makeSU3ContractV(left, right);
        return deviceSU3_12(full);
    }

    __device__ __inline__ void DebugPrint(const char* header = NULL) const
    {
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        printf("%s%s{{%f%s%fi, %f%s%fi, %f%s%fi},\n {%f%s%fi, %f%s%fi, %f%s%fi},\n {%f%s%fi, %f%s%fi, %f%s%fi}};\n",
            NULL == header ? "" : header,
            NULL == header ? "" : "=",
            m_me[0].x, m_me[0].y < 0 ? "" : "+", m_me[0].y,
            m_me[1].x, m_me[1].y < 0 ? "" : "+", m_me[1].y,
            m_me[2].x, m_me[2].y < 0 ? "" : "+", m_me[2].y,
            m_me[3].x, m_me[3].y < 0 ? "" : "+", m_me[3].y,
            m_me[4].x, m_me[4].y < 0 ? "" : "+", m_me[4].y,
            m_me[5].x, m_me[5].y < 0 ? "" : "+", m_me[5].y,
            c6.x, c6.y < 0 ? "" : "+", c6.y,
            c7.x, c7.y < 0 ? "" : "+", c7.y,
            c8.x, c8.y < 0 ? "" : "+", c8.y);
    }

#pragma endregion create

#pragma region conversion

    // Compress: deviceSU3 -> deviceSU3_12 (drop row2)
    // Both deviceSU3 and deviceSU3_12 use row-major: m_me[0..2]=row0, m_me[3..5]=row1
    // Row 2 is reconstructed via cross product: row2 = conj(row0 x row1)
    __device__ __inline__ explicit deviceSU3_12(const deviceSU3& full)
    {
        m_me[0] = full.m_me[0];
        m_me[1] = full.m_me[1];
        m_me[2] = full.m_me[2];
        m_me[3] = full.m_me[3];
        m_me[4] = full.m_me[4];
        m_me[5] = full.m_me[5];
    }

    // Expand: deviceSU3_12 -> deviceSU3 (reconstruct row2)
    __device__ __inline__ deviceSU3 toSU3() const
    {
        deviceSU3 ret;
        ret.m_me[0] = m_me[0];
        ret.m_me[1] = m_me[1];
        ret.m_me[2] = m_me[2];
        ret.m_me[3] = m_me[3];
        ret.m_me[4] = m_me[4];
        ret.m_me[5] = m_me[5];
        _reconstructCol2(ret.m_me[6], ret.m_me[7], ret.m_me[8]);
        return ret;
    }

    // Compress helper (free function style)
    __device__ __inline__ static deviceSU3_12 compress(const deviceSU3& full)
    {
        return deviceSU3_12(full);
    }

    // Expand helper (free function style)
    __device__ __inline__ static deviceSU3 expand(const deviceSU3_12& compact)
    {
        return compact.toSU3();
    }

#pragma endregion conversion

#pragma region operators_su3_preserving

    //--- Group multiply (SU(3) x SU(3) -> SU(3)) ---

    // this = this * right
    // SU3 is row-major: m_me[0..2]=row0, m_me[3..5]=row1, c6/c7/c8=row2
    // (A*B)(i,j) = sum_k A(i,k)*B(k,j)
    __device__ __inline__ void Mul(const deviceSU3_12& right)
    {
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        // Row 0: A(0,k) = {m_me[0], m_me[1], m_me[2]}
        CLGComplex res0 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[0]), _cuCmulf(m_me[1], right.m_me[3])), _cuCmulf(m_me[2], c6r));
        CLGComplex res1 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[1]), _cuCmulf(m_me[1], right.m_me[4])), _cuCmulf(m_me[2], c7r));
        CLGComplex res2 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[2]), _cuCmulf(m_me[1], right.m_me[5])), _cuCmulf(m_me[2], c8r));
        // Row 1: A(1,k) = {m_me[3], m_me[4], m_me[5]}
        CLGComplex res3 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[0]), _cuCmulf(m_me[4], right.m_me[3])), _cuCmulf(m_me[5], c6r));
        CLGComplex res4 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[1]), _cuCmulf(m_me[4], right.m_me[4])), _cuCmulf(m_me[5], c7r));
        CLGComplex res5 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[2]), _cuCmulf(m_me[4], right.m_me[5])), _cuCmulf(m_me[5], c8r));

        m_me[0] = res0; m_me[1] = res1; m_me[2] = res2;
        m_me[3] = res3; m_me[4] = res4; m_me[5] = res5;
    }

    // return this * right
    __device__ __inline__ deviceSU3_12 MulC(const deviceSU3_12& right) const
    {
        deviceSU3_12 ret;
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);

        // Row 0: A(0,k) = {m_me[0], m_me[1], m_me[2]}
        ret.m_me[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[0]), _cuCmulf(m_me[1], right.m_me[3])), _cuCmulf(m_me[2], c6r));
        ret.m_me[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[1]), _cuCmulf(m_me[1], right.m_me[4])), _cuCmulf(m_me[2], c7r));
        ret.m_me[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], right.m_me[2]), _cuCmulf(m_me[1], right.m_me[5])), _cuCmulf(m_me[2], c8r));
        // Row 1: A(1,k) = {m_me[3], m_me[4], m_me[5]}
        ret.m_me[3] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[0]), _cuCmulf(m_me[4], right.m_me[3])), _cuCmulf(m_me[5], c6r));
        ret.m_me[4] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[1]), _cuCmulf(m_me[4], right.m_me[4])), _cuCmulf(m_me[5], c7r));
        ret.m_me[5] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], right.m_me[2]), _cuCmulf(m_me[4], right.m_me[5])), _cuCmulf(m_me[5], c8r));
        return ret;
    }

    // this = this * right^dagger
    // (A*B^†)(i,j) = sum_k A(i,k)*conj(B(j,k))
    __device__ __inline__ void MulDagger(const deviceSU3_12& right)
    {
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);

        // Row 0: A(0,k) = {m_me[0], m_me[1], m_me[2]}
        CLGComplex res0 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(right.m_me[0])), _cuCmulf(m_me[1], _cuConjf(right.m_me[1]))), _cuCmulf(m_me[2], _cuConjf(right.m_me[2])));
        CLGComplex res1 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(right.m_me[3])), _cuCmulf(m_me[1], _cuConjf(right.m_me[4]))), _cuCmulf(m_me[2], _cuConjf(right.m_me[5])));
        CLGComplex res2 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(c6r)), _cuCmulf(m_me[1], _cuConjf(c7r))), _cuCmulf(m_me[2], _cuConjf(c8r)));
        // Row 1: A(1,k) = {m_me[3], m_me[4], m_me[5]}
        CLGComplex res3 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(right.m_me[0])), _cuCmulf(m_me[4], _cuConjf(right.m_me[1]))), _cuCmulf(m_me[5], _cuConjf(right.m_me[2])));
        CLGComplex res4 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(right.m_me[3])), _cuCmulf(m_me[4], _cuConjf(right.m_me[4]))), _cuCmulf(m_me[5], _cuConjf(right.m_me[5])));
        CLGComplex res5 = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(c6r)), _cuCmulf(m_me[4], _cuConjf(c7r))), _cuCmulf(m_me[5], _cuConjf(c8r)));

        m_me[0] = res0; m_me[1] = res1; m_me[2] = res2;
        m_me[3] = res3; m_me[4] = res4; m_me[5] = res5;
    }

    // return this * right^dagger
    __device__ __inline__ deviceSU3_12 MulDaggerC(const deviceSU3_12& right) const
    {
        deviceSU3_12 ret;
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);

        // Row 0: A(0,k) = {m_me[0], m_me[1], m_me[2]}
        ret.m_me[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(right.m_me[0])), _cuCmulf(m_me[1], _cuConjf(right.m_me[1]))), _cuCmulf(m_me[2], _cuConjf(right.m_me[2])));
        ret.m_me[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(right.m_me[3])), _cuCmulf(m_me[1], _cuConjf(right.m_me[4]))), _cuCmulf(m_me[2], _cuConjf(right.m_me[5])));
        ret.m_me[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], _cuConjf(c6r)), _cuCmulf(m_me[1], _cuConjf(c7r))), _cuCmulf(m_me[2], _cuConjf(c8r)));
        // Row 1: A(1,k) = {m_me[3], m_me[4], m_me[5]}
        ret.m_me[3] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(right.m_me[0])), _cuCmulf(m_me[4], _cuConjf(right.m_me[1]))), _cuCmulf(m_me[5], _cuConjf(right.m_me[2])));
        ret.m_me[4] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(right.m_me[3])), _cuCmulf(m_me[4], _cuConjf(right.m_me[4]))), _cuCmulf(m_me[5], _cuConjf(right.m_me[5])));
        ret.m_me[5] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], _cuConjf(c6r)), _cuCmulf(m_me[4], _cuConjf(c7r))), _cuCmulf(m_me[5], _cuConjf(c8r)));
        return ret;
    }

    // this = this^dagger * right
    // (A^†*B)(i,j) = sum_k conj(A(k,i))*B(k,j)
    __device__ __inline__ void DaggerMul(const deviceSU3_12& right)
    {
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        // Row 0: col i=0 -> conj(A(k,0)) = {conj(m_me[0]), conj(m_me[3]), conj(c6)}
        CLGComplex res0 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[0]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[3])), _cuCmulf(_cuConjf(c6), c6r));
        CLGComplex res1 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[1]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[4])), _cuCmulf(_cuConjf(c6), c7r));
        CLGComplex res2 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[2]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[5])), _cuCmulf(_cuConjf(c6), c8r));
        // Row 1: col i=1 -> conj(A(k,1)) = {conj(m_me[1]), conj(m_me[4]), conj(c7)}
        CLGComplex res3 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[0]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[3])), _cuCmulf(_cuConjf(c7), c6r));
        CLGComplex res4 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[1]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[4])), _cuCmulf(_cuConjf(c7), c7r));
        CLGComplex res5 = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[2]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[5])), _cuCmulf(_cuConjf(c7), c8r));

        m_me[0] = res0; m_me[1] = res1; m_me[2] = res2;
        m_me[3] = res3; m_me[4] = res4; m_me[5] = res5;
    }

    // return this^dagger * right
    __device__ __inline__ deviceSU3_12 DaggerMulC(const deviceSU3_12& right) const
    {
        deviceSU3_12 ret;
        CLGComplex c6r, c7r, c8r;
        right._reconstructCol2(c6r, c7r, c8r);
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        // Row 0: col i=0 -> conj(A(k,0)) = {conj(m_me[0]), conj(m_me[3]), conj(c6)}
        ret.m_me[0] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[0]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[3])), _cuCmulf(_cuConjf(c6), c6r));
        ret.m_me[1] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[1]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[4])), _cuCmulf(_cuConjf(c6), c7r));
        ret.m_me[2] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), right.m_me[2]), _cuCmulf(_cuConjf(m_me[3]), right.m_me[5])), _cuCmulf(_cuConjf(c6), c8r));
        // Row 1: col i=1 -> conj(A(k,1)) = {conj(m_me[1]), conj(m_me[4]), conj(c7)}
        ret.m_me[3] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[0]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[3])), _cuCmulf(_cuConjf(c7), c6r));
        ret.m_me[4] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[1]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[4])), _cuCmulf(_cuConjf(c7), c7r));
        ret.m_me[5] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), right.m_me[2]), _cuCmulf(_cuConjf(m_me[4]), right.m_me[5])), _cuCmulf(_cuConjf(c7), c8r));
        return ret;
    }

    //--- Vector multiply ---

    // A*v: result[i] = sum_k A(i,k)*v[k]
    __device__ __inline__ deviceSU3Vector MulVector(const deviceSU3Vector& v) const
    {
        deviceSU3Vector ret;
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        ret.m_ve[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], v.m_ve[0]), _cuCmulf(m_me[1], v.m_ve[1])), _cuCmulf(m_me[2], v.m_ve[2]));
        ret.m_ve[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], v.m_ve[0]), _cuCmulf(m_me[4], v.m_ve[1])), _cuCmulf(m_me[5], v.m_ve[2]));
        ret.m_ve[2] = _cuCaddf(_cuCaddf(_cuCmulf(c6, v.m_ve[0]), _cuCmulf(c7, v.m_ve[1])), _cuCmulf(c8, v.m_ve[2]));
        return ret;
    }

    // A^†*v: result[i] = sum_k conj(A(k,i))*v[k]
    __device__ __inline__ deviceSU3Vector DagMulVector(const deviceSU3Vector& v) const
    {
        deviceSU3Vector ret;
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        ret.m_ve[0] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[0]), v.m_ve[0]), _cuCmulf(_cuConjf(m_me[3]), v.m_ve[1])), _cuCmulf(_cuConjf(c6), v.m_ve[2]));
        ret.m_ve[1] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[1]), v.m_ve[0]), _cuCmulf(_cuConjf(m_me[4]), v.m_ve[1])), _cuCmulf(_cuConjf(c7), v.m_ve[2]));
        ret.m_ve[2] = _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[2]), v.m_ve[0]), _cuCmulf(_cuConjf(m_me[5]), v.m_ve[1])), _cuCmulf(_cuConjf(c8), v.m_ve[2]));
        return ret;
    }

    __device__ __inline__ deviceWilsonVectorSU3 MulWilsonVector(const deviceWilsonVectorSU3& v) const
    {
        deviceWilsonVectorSU3 ret;
        ret.m_d[0] = MulVector(v.m_d[0]);
        ret.m_d[1] = MulVector(v.m_d[1]);
        ret.m_d[2] = MulVector(v.m_d[2]);
        ret.m_d[3] = MulVector(v.m_d[3]);
        return ret;
    }

    __device__ __inline__ deviceWilsonVectorSU3 DagMulWilsonVector(const deviceWilsonVectorSU3& v) const
    {
        deviceWilsonVectorSU3 ret;
        ret.m_d[0] = DagMulVector(v.m_d[0]);
        ret.m_d[1] = DagMulVector(v.m_d[1]);
        ret.m_d[2] = DagMulVector(v.m_d[2]);
        ret.m_d[3] = DagMulVector(v.m_d[3]);
        return ret;
    }

    //--- Conjugate transpose (SU(3) -> SU(3)) ---

    __device__ __inline__ void Dagger()
    {
        // Swap [0]<->[0], [1]<->[3], [2]<->[6]
        // But we only have 6 elements. After dagger:
        // new[0] = conj(old[0]), new[1] = conj(old[3]), new[2] = conj(old[6])
        // new[3] = conj(old[1]), new[4] = conj(old[4]), new[5] = conj(old[7])
        // old[6] = conj(cross(old_col0, old_col1)), old[7] = conj(cross(...))
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        CLGComplex new0 = _cuConjf(m_me[0]);
        CLGComplex new1 = _cuConjf(m_me[3]);
        CLGComplex new2 = _cuConjf(c6);
        CLGComplex new3 = _cuConjf(m_me[1]);
        CLGComplex new4 = _cuConjf(m_me[4]);
        CLGComplex new5 = _cuConjf(c7);

        m_me[0] = new0;
        m_me[1] = new1;
        m_me[2] = new2;
        m_me[3] = new3;
        m_me[4] = new4;
        m_me[5] = new5;
    }

    __device__ __inline__ deviceSU3_12 DaggerC() const
    {
        deviceSU3_12 ret;
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        ret.m_me[0] = _cuConjf(m_me[0]);
        ret.m_me[1] = _cuConjf(m_me[3]);
        ret.m_me[2] = _cuConjf(c6);
        ret.m_me[3] = _cuConjf(m_me[1]);
        ret.m_me[4] = _cuConjf(m_me[4]);
        ret.m_me[5] = _cuConjf(c7);
        return ret;
    }

    //--- Transpose (not conjugate, NOT SU(3) preserving in general, but for SU(3) it's fine since Transpose = Dagger for unitary) ---
    // Actually Transpose is NOT the same as Dagger. For SU(3), Transpose gives U^T which is NOT SU(3) in general.
    // However, the C version returns deviceSU3 so it's safe.

    __device__ __inline__ deviceSU3 Transpose() const
    {
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);

        deviceSU3 ret;
        ret.m_me[0] = m_me[0];
        ret.m_me[1] = m_me[3];
        ret.m_me[2] = c6;
        ret.m_me[3] = m_me[1];
        ret.m_me[4] = m_me[4];
        ret.m_me[5] = c7;
        ret.m_me[6] = m_me[2];
        ret.m_me[7] = m_me[5];
        ret.m_me[8] = c8;
        return ret;
    }

    //--- Projection (SU(3) -> SU(3), no-op for SU(3) input) ---

    __device__ __inline__ void Norm()
    {
        // For SU(3) input, this is a no-op.
        // For non-SU(3) input, we need to expand, project, compress.
        // Since we assume SU(3) input, do nothing.
    }

    __device__ __inline__ void SU3Proj()
    {
        // Same as Norm - no-op for SU(3) input.
    }

    __device__ __inline__ void U3Proj()
    {
        // No-op for SU(3) input.
    }

    __device__ __inline__ void Proj(BYTE ite = 4)
    {
        // No-op for SU(3) input.
    }

    __device__ __inline__ void CabbiboMarinariProj()
    {
        // No-op for SU(3) input.
    }

    //--- Inverse (SU(3) -> SU(3)) ---

    __device__ __inline__ deviceSU3_12 Inverse() const
    {
        // For SU(3), U^{-1} = U^dagger
        return DaggerC();
    }

    //--- Exponential maps (SU(3) algebra -> SU(3)) ---

    // exp(a * this) for anti-Hermitian this, analytic closed-form
    __device__ __inline__ deviceSU3_12 QuickExp(Real a) const
    {
        deviceSU3 full = toSU3();
        deviceSU3 result = full.QuickExp(a);
        return deviceSU3_12(result);
    }

    // exp(i*a*this) for anti-Hermitian this, using hep-lat/0311018
    __device__ __inline__ deviceSU3_12 StrictExpTA(Real a) const
    {
        deviceSU3 full = toSU3();
        deviceSU3 result = full.StrictExpTA(a);
        return deviceSU3_12(result);
    }

    //--- Scalar queries (work on any matrix, return scalars) ---

    __device__ __inline__ CLGComplex Tr() const
    {
        // tr = m0 + m4 + m8, where m8 = conj(cross product)
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        return _cuCaddf(m_me[0], _cuCaddf(m_me[4], c8));
    }

    __device__ __inline__ Real ReTr() const
    {
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        return m_me[0].x + m_me[4].x + c8.x;
    }

    __device__ __inline__ Real ImTr() const
    {
        CLGComplex c6, c7, c8;
        _reconstructCol2(c6, c7, c8);
        return m_me[0].y + m_me[4].y + c8.y;
    }

    __device__ __inline__ CLGComplex Det() const
    {
        // For SU(3), det = 1
        return _onec;
    }

    __device__ __inline__ cuDoubleComplex DoubleDet() const
    {
        return make_cuDoubleComplex(1.0, 0.0);
    }

    __device__ __inline__ static Real TrIm(const deviceSU3_12& a, const deviceSU3_12& b)
    {
        CLGComplex c6a, c7a, c8a;
        a._reconstructCol2(c6a, c7a, c8a);
        CLGComplex c6b, c7b, c8b;
        b._reconstructCol2(c6b, c7b, c8b);

        Real ret = a.m_me[0].y * b.m_me[0].y + a.m_me[1].y * b.m_me[3].y + a.m_me[2].y * c6b.y;
        ret += a.m_me[3].y * b.m_me[1].y + a.m_me[4].y * b.m_me[4].y + a.m_me[5].y * c7b.y;
        ret += c6a.y * b.m_me[2].y + c7a.y * b.m_me[5].y + c8a.y * c8b.y;
        return ret;
    }

#pragma endregion operators_su3_preserving

#pragma region operators_non_su3_return_deviceSU3

    // These operations don't preserve SU(3), so they return deviceSU3 (full 18 reals).

    __device__ __inline__ deviceSU3 AddC(const deviceSU3_12& right) const
    {
        deviceSU3 l = toSU3();
        deviceSU3 r = right.toSU3();
        l.Add(r);
        return l;
    }

    __device__ __inline__ deviceSU3 AddDaggerC(const deviceSU3_12& right) const
    {
        deviceSU3 l = toSU3();
        deviceSU3 r = right.toSU3();
        l.AddDagger(r);
        return l;
    }

    __device__ __inline__ deviceSU3 SubC(const deviceSU3_12& right) const
    {
        deviceSU3 l = toSU3();
        deviceSU3 r = right.toSU3();
        l.Sub(r);
        return l;
    }

    __device__ __inline__ deviceSU3 SubDaggerC(const deviceSU3_12& right) const
    {
        deviceSU3 l = toSU3();
        deviceSU3 r = right.toSU3();
        l.SubDagger(r);
        return l;
    }

    __device__ __inline__ deviceSU3 AddRealC(Real right) const
    {
        deviceSU3 ret = toSU3();
        ret.AddReal(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 AddCompC(const CLGComplex& right) const
    {
        deviceSU3 ret = toSU3();
        ret.AddComp(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 SubRealC(Real right) const
    {
        deviceSU3 ret = toSU3();
        ret.SubReal(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 SubCompC(const CLGComplex& right) const
    {
        deviceSU3 ret = toSU3();
        ret.SubComp(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 MulRealC(Real right) const
    {
        deviceSU3 ret = toSU3();
        ret.MulReal(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 MulCompC(const CLGComplex& right) const
    {
        deviceSU3 ret = toSU3();
        ret.MulComp(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 DivCompC(const CLGComplex& right) const
    {
        deviceSU3 ret = toSU3();
        ret.DivComp(right);
        return ret;
    }

    __device__ __inline__ deviceSU3 OppositeC() const
    {
        deviceSU3 ret = toSU3();
        ret.Opposite();
        return ret;
    }

    __device__ __inline__ deviceSU3 ReC() const
    {
        deviceSU3 ret = toSU3();
        ret.Re();
        return ret;
    }

    __device__ __inline__ deviceSU3 ImC() const
    {
        deviceSU3 ret = toSU3();
        ret.Im();
        return ret;
    }

    __device__ __inline__ deviceSU3 TaC() const
    {
        deviceSU3 ret = toSU3();
        ret.Ta();
        return ret;
    }

    __device__ __inline__ deviceSU3 ThC() const
    {
        deviceSU3 ret = toSU3();
        ret.Th();
        return ret;
    }

    __device__ __inline__ deviceSU3 iIm2C() const
    {
        deviceSU3 ret = toSU3();
        ret.iIm2();
        return ret;
    }

    __device__ __inline__ deviceSU3 Re2C() const
    {
        deviceSU3 ret = toSU3();
        ret.Re2();
        return ret;
    }

    __device__ __inline__ deviceSU3 Im2C() const
    {
        deviceSU3 ret = toSU3();
        return ret.Im2C();
    }

    //--- General exponential/log/power (return deviceSU3) ---

    __device__ __inline__ deviceSU3 Exp(const CLGComplex& a, BYTE uiPrecision) const
    {
        deviceSU3 full = toSU3();
        return full.Exp(a, uiPrecision);
    }

    __device__ __inline__ deviceSU3 ExpReal(Real a, BYTE uiPrecision) const
    {
        deviceSU3 full = toSU3();
        return full.ExpReal(a, uiPrecision);
    }

    __device__ __inline__ deviceSU3 StrictExp() const
    {
        deviceSU3 full = toSU3();
        return full.StrictExp();
    }

    __device__ __inline__ deviceSU3 Power(Real fPower) const
    {
        deviceSU3 full = toSU3();
        return full.Power(fPower);
    }

    __device__ __inline__ deviceSU3 Log() const
    {
        deviceSU3 full = toSU3();
        return full.Log();
    }

    __device__ __inline__ void CalculateEigenValues(CLGComplex& c1, CLGComplex& c2, CLGComplex& c3) const
    {
        deviceSU3 full = toSU3();
        full.CalculateEigenValues(c1, c2, c3);
    }

    __device__ __inline__ deviceSU3 EigenVectors(
        const CLGComplex& c1, const CLGComplex& c2, const CLGComplex& c3) const
    {
        deviceSU3 full = toSU3();
        return full.EigenVectors(c1, c2, c3);
    }

#pragma endregion operators_non_su3_return_deviceSU3

#pragma region static_utility

    static __device__ __inline__ void deviceUVW(const deviceSU3_12& Q, const deviceSU3_12& Q2, DOUBLE& u, DOUBLE& v, DOUBLE& w)
    {
        deviceSU3 qFull = Q.toSU3();
        deviceSU3 q2Full = Q2.toSU3();
        deviceSU3::deviceUVW(qFull, q2Full, u, v, w);
    }

    static __device__ __inline__ void deviceSqrtf012(DOUBLE& f0, DOUBLE& f1, DOUBLE& f2, const DOUBLE& u, const DOUBLE& v, const DOUBLE& w)
    {
        deviceSU3::deviceSqrtf012(f0, f1, f2, u, v, w);
    }

    static __device__ __inline__ CLGComplex Determinent(const CLGComplex* u)
    {
        return deviceSU3::Determinent(u);
    }

#pragma endregion static_utility

    //--- Unsafe mutating operations (not supported) ---
    // The following operations from deviceSU3 are NOT supported on deviceSU3_12
    // because they don't preserve SU(3) and mutate in-place:
    //   Add, AddDagger, Sub, SubDagger, AddReal, AddComp, SubReal, SubComp, AddId
    //   MulReal, MulComp, DivComp, Opposite, Re, Im, Ta, Th, iIm2, Re2
    // Use the C-suffix versions (AddC, SubC, etc.) which return deviceSU3 instead.

    // Conversion free functions
    static __device__ __inline__ deviceSU3_12 deviceSU3_to_deviceSU3_12(const deviceSU3& full)
    {
        return deviceSU3_12(full);
    }

    static __device__ __inline__ deviceSU3 deviceSU3_12_to_deviceSU3(const deviceSU3_12& compact)
    {
        return compact.toSU3();
    }
};

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _SU3_12_H_

//=============================================================================
// END OF FILE
//=============================================================================
