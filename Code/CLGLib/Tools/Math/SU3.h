//=============================================================================
// FILENAME : SU3.h
// 
// DESCRIPTION:
// This is helper functions for calculate D=3 complex matrix
// NOTE: It is not nessary a SU3
//
// The SU3 Matrix is
// 0 1 2
// 3 4 5
// 6 7 8
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _SU3_H_
#define _SU3_H_

#define __LINE_MUL(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(m_me[a], right.m_me[d]), _cuCmulf(m_me[b], right.m_me[ee])), _cuCmulf(m_me[c], right.m_me[ff]))

#define __LINE_MULND(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(m_me[a], _cuConjf(right.m_me[d])), _cuCmulf(m_me[b], _cuConjf(right.m_me[ee]))), _cuCmulf(m_me[c], _cuConjf(right.m_me[ff])))

#define __LINE_MULDN(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[a]), right.m_me[d]), _cuCmulf(_cuConjf(m_me[b]), right.m_me[ee])), _cuCmulf(_cuConjf(m_me[c]), _cuConjf(right.m_me[ff]))

// 1.0f / _sqrt(3)
#define InvSqrt3 (F(0.5773502691896258))
// 2.0f / _sqrt(3)
#define InvSqrt3_2 (F(1.1547005383792517))

#if _CLG_DOUBLEFLOAT
#define __SU3MATRIX_ALIGN 256
#else
#define __SU3MATRIX_ALIGN 128
#endif

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    struct alignas(__SU3MATRIX_ALIGN) deviceSU3
    {
        __device__ deviceSU3()
        {

        }

        __device__ deviceSU3(const deviceSU3& other)
        {
            memcpy(m_me, other.m_me, sizeof(_Complex) * 9);
        }

        __device__ void DebugPrint() const
        {
            printf("%2.3f %s %2.3f I, %2.3f %s %2.3f I, %2.3f %s %2.3f I;\n%2.3f %s %2.3f I, %2.3f %s %2.3f I, %2.3f %s %2.3f I;\n%2.3f %s %2.3f I, %2.3f %s %2.3f I, %2.3f %s %2.3f I;\n",
                m_me[0].x,
                m_me[0].y < 0 ? "" : "+",
                m_me[0].y,

                m_me[1].x,
                m_me[1].y < 0 ? "" : "+",
                m_me[1].y,

                m_me[2].x,
                m_me[2].y < 0 ? "" : "+",
                m_me[2].y,

                m_me[3].x,
                m_me[3].y < 0 ? "" : "+",
                m_me[3].y,

                m_me[4].x,
                m_me[4].y < 0 ? "" : "+",
                m_me[4].y,

                m_me[5].x,
                m_me[5].y < 0 ? "" : "+",
                m_me[5].y,

                m_me[6].x,
                m_me[6].y < 0 ? "" : "+",
                m_me[6].y,

                m_me[7].x,
                m_me[7].y < 0 ? "" : "+",
                m_me[7].y,

                m_me[8].x,
                m_me[8].y < 0 ? "" : "+",
                m_me[8].y
            );
        }

#pragma region create

        /**
        * ret = 0
        */
        __device__ __inline__ static deviceSU3 makeSU3Zero()
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[4] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[5] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[6] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[7] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[8] = _make_cuComplex(F(0.0), F(0.0));
            return ret;
        }

        /**
        * ret = I
        */
        __device__ __inline__ static deviceSU3 makeSU3Id()
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
            ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[4] = _make_cuComplex(F(1.0), F(0.0));
            ret.m_me[5] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[6] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[7] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[8] = _make_cuComplex(F(1.0), F(0.0));
            return ret;
        }

        /**
        * can be called only after CLatticeData is created
        * ret = random
        */
        __device__ __inline__ static deviceSU3 makeSU3Random(UINT fatIndex)
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[1] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[2] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[3] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[4] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[5] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[6] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[7] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[8] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.Norm();
            return ret;
        }

        /**
        * Here, we keep same with Bridge++, that H(P)/D.O.F. = 0.5
        * can be called only after CLatticeData is created
        * ret = r_a T_a, r_a is random real number, T_a are generators
        * r_a T_a = r1 T1 + r2 T2 + ...
        *
        *     r3 +r8/sqrt3     r1-ir2        r4-ir5
        * =   r1+ir2       -r3+r8/sqrt3      r6-ir7
        *     r4+ir5           r6+ir7      -2r8/sqrt3
        *
        */
        __device__ __inline__ static deviceSU3 makeSU3RandomGenerator(UINT fatIndex)
        {
            Real r1 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r2 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r3 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r4 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r5 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r6 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r7 = _deviceRandomGaussF(fatIndex) * InvSqrt2;
            Real r8 = _deviceRandomGaussF(fatIndex) * InvSqrt2;

            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(r3 + r8 * InvSqrt3, F(0.0));
            ret.m_me[1] = _make_cuComplex(r1, -r2);
            ret.m_me[2] = _make_cuComplex(r4, -r5);
            ret.m_me[3] = _make_cuComplex(r1, r2);
            ret.m_me[4] = _make_cuComplex(-r3 + r8 * InvSqrt3, F(0.0));
            ret.m_me[5] = _make_cuComplex(r6, -r7);
            ret.m_me[6] = _make_cuComplex(r4, r5);
            ret.m_me[7] = _make_cuComplex(r6, r7);
            ret.m_me[8] = _make_cuComplex(-r8 * InvSqrt3_2, F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3 makeSU3SumGenerator(Real fDivide)
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(fDivide + fDivide * InvSqrt3, F(0.0));
            ret.m_me[1] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[2] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[3] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[4] = _make_cuComplex(-fDivide + fDivide * InvSqrt3, F(0.0));
            ret.m_me[5] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[6] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[7] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[8] = _make_cuComplex(-fDivide * InvSqrt3_2, F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3* makeSU3Generator(UINT uiGenerator)
        {
            deviceSU3* ret = new deviceSU3();
            Real half = F(1.0) / F(0.5);
            ret->makeSU3Zero();
            switch (uiGenerator)
            {
            case 0:
            {
                /**
                *     0     1     0
                * =   1     0     0
                *     0     0     0
                */
                ret->m_me[1] = _make_cuComplex(half, F(0.0));
                ret->m_me[3] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 1:
            {
                /**
                *     0     -I     0
                * =   I     0     0
                *     0     0     0
                */
                ret->m_me[1] = _make_cuComplex(F(0.0), -half);
                ret->m_me[3] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 2:
            {
                /**
                *     1     0     0
                * =   0    -1     0
                *     0     0     0
                */
                ret->m_me[0] = _make_cuComplex(half, F(0.0));
                ret->m_me[4] = _make_cuComplex(-half, F(0.0));
            }
            break;
            case 3:
            {
                /**
                *     0     0     1
                * =   0     0     0
                *     1     0     0
                */
                ret->m_me[2] = _make_cuComplex(half, F(0.0));
                ret->m_me[6] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 4:
            {
                /**
                *     0     0    -i
                * =   0     0     0
                *     i     0     0
                */
                ret->m_me[2] = _make_cuComplex(F(0.0), -half);
                ret->m_me[6] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 5:
            {
                /**
                *     0     0     0
                * =   0     0     1
                *     0     1     0
                */
                ret->m_me[5] = _make_cuComplex(half, F(0.0));
                ret->m_me[7] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 6:
            {
                /**
                *     0     0     0
                * =   0     0    -i
                *     0     i     0
                */
                ret->m_me[5] = _make_cuComplex(F(0.0), -half);
                ret->m_me[7] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 7:
            {
                /**
                *     1     0     0
                * =   0     1     0
                *     0     0     -2
                */
                ret->m_me[0] = _make_cuComplex(InvSqrt3 * half, F(0.0));
                ret->m_me[4] = _make_cuComplex(InvSqrt3 * half, F(0.0));
                ret->m_me[8] = _make_cuComplex(-InvSqrt3_2 * half, F(0.0));
            }
            break;
            }
            return ret;
        }

#pragma endregion create

#pragma region operators

        __device__ __inline__ void Add(const deviceSU3& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right.m_me[0]);
            m_me[1] = _cuCaddf(m_me[1], right.m_me[1]);
            m_me[2] = _cuCaddf(m_me[2], right.m_me[2]);
            m_me[3] = _cuCaddf(m_me[3], right.m_me[3]);
            m_me[4] = _cuCaddf(m_me[4], right.m_me[4]);
            m_me[5] = _cuCaddf(m_me[5], right.m_me[5]);
            m_me[6] = _cuCaddf(m_me[6], right.m_me[6]);
            m_me[7] = _cuCaddf(m_me[7], right.m_me[7]);
            m_me[8] = _cuCaddf(m_me[8], right.m_me[8]);
        }

        __device__ __inline__ void Sub(const deviceSU3& right)
        {
            m_me[0] = _cuCsubf(m_me[0], right.m_me[0]);
            m_me[1] = _cuCsubf(m_me[1], right.m_me[1]);
            m_me[2] = _cuCsubf(m_me[2], right.m_me[2]);
            m_me[3] = _cuCsubf(m_me[3], right.m_me[3]);
            m_me[4] = _cuCsubf(m_me[4], right.m_me[4]);
            m_me[5] = _cuCsubf(m_me[5], right.m_me[5]);
            m_me[6] = _cuCsubf(m_me[6], right.m_me[6]);
            m_me[7] = _cuCsubf(m_me[7], right.m_me[7]);
            m_me[8] = _cuCsubf(m_me[8], right.m_me[8]);
        }

        __device__ __inline__ void AddReal(const Real& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right);
            m_me[4] = _cuCaddf(m_me[4], right);
            m_me[8] = _cuCaddf(m_me[8], right);
        }

        __device__ __inline__ void AddComp(const _Complex& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right);
            m_me[4] = _cuCaddf(m_me[4], right);
            m_me[8] = _cuCaddf(m_me[8], right);
        }

        __device__ __inline__ void SubReal(const Real& right)
        {
            m_me[0] = _cuCsubf(m_me[0], right);
            m_me[4] = _cuCsubf(m_me[4], right);
            m_me[8] = _cuCsubf(m_me[8], right);
        }

        __device__ __inline__ void SubComp(const _Complex& right)
        {
            m_me[0] = _cuCsubf(m_me[0], right);
            m_me[4] = _cuCsubf(m_me[4], right);
            m_me[8] = _cuCsubf(m_me[8], right);
        }

        /**
        * left = left * right
        */
        __device__ __inline__ void Mul(const deviceSU3& right)
        {
            _Complex res[9];
            res[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
            res[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
            res[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

            res[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
            res[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
            res[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

            res[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
            res[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
            res[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
            memcpy(m_me, res, sizeof(_Complex) * 9);
        }

        __device__ __inline__ deviceSU3 MulC(const deviceSU3& right) const
        {
            deviceSU3 ret;
            ret.m_me[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
            ret.m_me[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
            ret.m_me[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

            ret.m_me[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
            ret.m_me[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
            ret.m_me[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

            ret.m_me[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
            ret.m_me[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
            ret.m_me[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
            return ret;
        }

        __device__ __inline__ void MulReal(Real right)
        {
            m_me[0] = _cuCmulf(m_me[0], right);
            m_me[1] = _cuCmulf(m_me[1], right);
            m_me[2] = _cuCmulf(m_me[2], right);
            m_me[3] = _cuCmulf(m_me[3], right);
            m_me[4] = _cuCmulf(m_me[4], right);
            m_me[5] = _cuCmulf(m_me[5], right);
            m_me[6] = _cuCmulf(m_me[6], right);
            m_me[7] = _cuCmulf(m_me[7], right);
            m_me[8] = _cuCmulf(m_me[8], right);
        }

        __device__ __inline__ void MulComp(const _Complex& right)
        {
            m_me[0] = _cuCmulf(m_me[0], right);
            m_me[1] = _cuCmulf(m_me[1], right);
            m_me[2] = _cuCmulf(m_me[2], right);
            m_me[3] = _cuCmulf(m_me[3], right);
            m_me[4] = _cuCmulf(m_me[4], right);
            m_me[5] = _cuCmulf(m_me[5], right);
            m_me[6] = _cuCmulf(m_me[6], right);
            m_me[7] = _cuCmulf(m_me[7], right);
            m_me[8] = _cuCmulf(m_me[8], right);
        }

        __device__ __inline__ void MulDagger(const deviceSU3& right)
        {
            _Complex res[9];
            res[0] = __LINE_MULND(0, 1, 2, 0, 1, 2);
            res[1] = __LINE_MULND(0, 1, 2, 3, 4, 5);
            res[2] = __LINE_MULND(0, 1, 2, 6, 7, 8);

            res[3] = __LINE_MULND(3, 4, 5, 0, 1, 2);
            res[4] = __LINE_MULND(3, 4, 5, 3, 4, 5);
            res[5] = __LINE_MULND(3, 4, 5, 6, 7, 8);

            res[6] = __LINE_MULND(6, 7, 8, 0, 1, 2);
            res[7] = __LINE_MULND(6, 7, 8, 3, 4, 5);
            res[8] = __LINE_MULND(6, 7, 8, 6, 7, 8);
            memcpy(m_me, res, sizeof(_Complex) * 9);
        }

        __device__ __inline__ deviceSU3 MulDaggerc(const deviceSU3& right) const
        {
            deviceSU3 ret;
            ret.m_me[0] = __LINE_MULND(0, 1, 2, 0, 1, 2);
            ret.m_me[1] = __LINE_MULND(0, 1, 2, 3, 4, 5);
            ret.m_me[2] = __LINE_MULND(0, 1, 2, 6, 7, 8);

            ret.m_me[3] = __LINE_MULND(3, 4, 5, 0, 1, 2);
            ret.m_me[4] = __LINE_MULND(3, 4, 5, 3, 4, 5);
            ret.m_me[5] = __LINE_MULND(3, 4, 5, 6, 7, 8);

            ret.m_me[6] = __LINE_MULND(6, 7, 8, 0, 1, 2);
            ret.m_me[7] = __LINE_MULND(6, 7, 8, 3, 4, 5);
            ret.m_me[8] = __LINE_MULND(6, 7, 8, 6, 7, 8);
            return ret;
        }

        __device__ __inline__ void DaggerMul(const deviceSU3& right)
        {
            _Complex res[9];
            res[0] = __LINE_MULND(0, 3, 6, 0, 3, 6);
            res[1] = __LINE_MULND(0, 3, 6, 1, 4, 7);
            res[2] = __LINE_MULND(0, 3, 6, 2, 5, 8);

            res[3] = __LINE_MULND(1, 4, 7, 0, 3, 6);
            res[4] = __LINE_MULND(1, 4, 7, 1, 4, 7);
            res[5] = __LINE_MULND(1, 4, 7, 2, 5, 8);

            res[6] = __LINE_MULND(2, 5, 8, 0, 3, 6);
            res[7] = __LINE_MULND(2, 5, 8, 1, 4, 7);
            res[8] = __LINE_MULND(2, 5, 8, 2, 5, 8);
            memcpy(m_me, res, sizeof(_Complex) * 9);
        }

        __device__ __inline__ deviceSU3 DaggerMulc(const deviceSU3& right) const
        {
            deviceSU3 ret;
            ret.m_me[0] = __LINE_MULND(0, 3, 6, 0, 3, 6);
            ret.m_me[1] = __LINE_MULND(0, 3, 6, 1, 4, 7);
            ret.m_me[2] = __LINE_MULND(0, 3, 6, 2, 5, 8);

            ret.m_me[3] = __LINE_MULND(1, 4, 7, 0, 3, 6);
            ret.m_me[4] = __LINE_MULND(1, 4, 7, 1, 4, 7);
            ret.m_me[5] = __LINE_MULND(1, 4, 7, 2, 5, 8);

            ret.m_me[6] = __LINE_MULND(2, 5, 8, 0, 3, 6);
            ret.m_me[7] = __LINE_MULND(2, 5, 8, 1, 4, 7);
            ret.m_me[8] = __LINE_MULND(2, 5, 8, 2, 5, 8);
            return ret;
        }

        __device__ __inline__ deviceSU3Vector MulVector(const deviceSU3Vector& v) const
        {
            deviceSU3Vector ret;
            ret.m_ve[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], v.m_ve[0]), _cuCmulf(m_me[1], v.m_ve[1])), _cuCmulf(m_me[2], v.m_ve[2]));
            ret.m_ve[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], v.m_ve[0]), _cuCmulf(m_me[4], v.m_ve[1])), _cuCmulf(m_me[5], v.m_ve[2]));
            ret.m_ve[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[6], v.m_ve[0]), _cuCmulf(m_me[7], v.m_ve[1])), _cuCmulf(m_me[8], v.m_ve[2]));
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

        __device__ __inline__ deviceSU3 AddC(const deviceSU3& right) const { deviceSU3 ret(*this); ret.Add(right); return ret; }
        __device__ __inline__ deviceSU3 AddCompC(const _Complex& right) const { deviceSU3 ret(*this); ret.AddComp(right); return ret; }
        __device__ __inline__ deviceSU3 AddRealC(const Real& right) const { deviceSU3 ret(*this); ret.AddReal(right); return ret; }

        __device__ __inline__ deviceSU3 SubC(const deviceSU3& right) const { deviceSU3 ret(*this); ret.Sub(right); return ret; }
        __device__ __inline__ deviceSU3 SubCompC(const _Complex& right) const { deviceSU3 ret(*this); ret.SubComp(right); return ret; }
        __device__ __inline__ deviceSU3 SubRealC(const Real& right) const { deviceSU3 ret(*this); ret.SubReal(right); return ret; }

        __device__ __inline__ deviceSU3 MulCompC(const _Complex& right) const { deviceSU3 ret(*this); ret.MulComp(right); return ret; }
        __device__ __inline__ deviceSU3 MulRealC(const Real& right) const { deviceSU3 ret(*this); ret.MulReal(right); return ret; }

#pragma endregion operators

#pragma region useful functions

        /**
        * ret = ||U||, where U is a d=3 matrix
        * we do NOT need ||SU3|| because ||SU3||=1
        */
        __device__ static __inline__ _Complex Determinent(const _Complex* u)
        {
            return _cuCsubf(
                _cuCaddf(
                    _cuCaddf(
                        _cuCmulf(_cuCmulf(u[0], u[4]), u[8]),
                        _cuCmulf(_cuCmulf(u[1], u[5]), u[6])
                    ),
                    _cuCmulf(_cuCmulf(u[3], u[7]), u[2])
                ),
                _cuCaddf(
                    _cuCaddf(
                        _cuCmulf(_cuCmulf(u[2], u[4]), u[6]),
                        _cuCmulf(_cuCmulf(u[1], u[3]), u[8])
                    ),
                    _cuCmulf(_cuCmulf(u[0], u[5]), u[7])
                )
            );
        }

        __device__ __inline__ _Complex Tr() const
        {
            return _cuCaddf(m_me[0], _cuCaddf(m_me[4], m_me[8]));
        }

        /**
        * Re[Tr[U]]
        */
        __device__ __inline__ Real ReTr() const
        {
            return m_me[0].x + m_me[4].x + m_me[8].x;
        }

        /**
        * res = Conjugate[Transpose[U]]
        */
        __device__ __inline__ void Dagger()
        {
            _Complex res[9];
            res[0] = _cuConjf(m_me[0]);
            res[1] = _cuConjf(m_me[3]);
            res[2] = _cuConjf(m_me[6]);

            res[3] = _cuConjf(m_me[1]);
            res[4] = _cuConjf(m_me[4]);
            res[5] = _cuConjf(m_me[7]);

            res[6] = _cuConjf(m_me[2]);
            res[7] = _cuConjf(m_me[5]);
            res[8] = _cuConjf(m_me[8]);

            memcpy(m_me, res, sizeof(_Complex) * 9);
        }

        __device__ __inline__ deviceSU3 DaggerC() const
        {
            deviceSU3 tmp;
            tmp.m_me[0] = _cuConjf(m_me[0]);
            tmp.m_me[1] = _cuConjf(m_me[3]);
            tmp.m_me[2] = _cuConjf(m_me[6]);

            tmp.m_me[3] = _cuConjf(m_me[1]);
            tmp.m_me[4] = _cuConjf(m_me[4]);
            tmp.m_me[5] = _cuConjf(m_me[7]);

            tmp.m_me[6] = _cuConjf(m_me[2]);
            tmp.m_me[7] = _cuConjf(m_me[5]);
            tmp.m_me[8] = _cuConjf(m_me[8]);
            return tmp;
        }

        /**
        * A=(U-U^dagger)/2
        * res = A - tr(A)/3
        */
        __device__ __inline__ void Ta()
        {
            //0 1 2
            //3 4 5
            //6 7 8
            //new [0], [4], [8] = 2img I / 2
            //we DO NOT calculate it here
            //m_me[0] = _make_cuComplex(0.0f, _cuCimagf(m_me[0]));
            //m_me[4] = _make_cuComplex(0.0f, _cuCimagf(m_me[4]));
            //m_me[8] = _make_cuComplex(0.0f, _cuCimagf(m_me[8]));

            //new [1] = [1] - conj([3])
            //new [3] = [3] - conj([1]) = -conj(new [1])
            _Complex new1 = _cuCsubf(m_me[1], _cuConjf(m_me[3]));
            _Complex new2 = _cuCsubf(m_me[2], _cuConjf(m_me[6]));
            _Complex new5 = _cuCsubf(m_me[5], _cuConjf(m_me[7]));
            m_me[1] = _make_cuComplex(F(0.5) * _cuCrealf(new1), F(0.5) * _cuCimagf(new1));
            m_me[3] = _make_cuComplex(-_cuCrealf(m_me[1]), _cuCimagf(m_me[1]));
            m_me[2] = _make_cuComplex(F(0.5) * _cuCrealf(new2), F(0.5) * _cuCimagf(new2));
            m_me[6] = _make_cuComplex(-_cuCrealf(m_me[2]), _cuCimagf(m_me[2]));
            m_me[5] = _make_cuComplex(F(0.5) * _cuCrealf(new5), F(0.5) * _cuCimagf(new5));
            m_me[7] = _make_cuComplex(-_cuCrealf(m_me[5]), _cuCimagf(m_me[5]));

            Real fImgTr = (m_me[0].y + m_me[4].y + m_me[8].y) / F(3.0);
            m_me[0] = _make_cuComplex(F(0.0), m_me[0].y - fImgTr);
            m_me[4] = _make_cuComplex(F(0.0), m_me[4].y - fImgTr);
            m_me[8] = _make_cuComplex(F(0.0), m_me[8].y - fImgTr);
        }

        /**
        * make any matrix to SU3
        */
        __device__ __inline__ void Norm()
        {
            Real fDiv1 = 1.0f / _sqrt(cuCabsSqf(m_me[0]) + cuCabsSqf(m_me[1]) + cuCabsSqf(m_me[2]));
            m_me[0] = _cuCmulf(m_me[0], fDiv1);
            m_me[1] = _cuCmulf(m_me[1], fDiv1);
            m_me[2] = _cuCmulf(m_me[2], fDiv1);

            //it is the name in Bridge++
            _Complex sp1 = _cuConjf(
                _cuCaddf(_cuCaddf(
                    _cuCmulf(m_me[0], _cuConjf(m_me[3]))
                    , _cuCmulf(m_me[1], _cuConjf(m_me[4])))
                    , _cuCmulf(m_me[2], _cuConjf(m_me[5])))
            );

            m_me[3] = _cuCsubf(m_me[3], _cuCmulf(sp1, m_me[0]));
            m_me[4] = _cuCsubf(m_me[4], _cuCmulf(sp1, m_me[1]));
            m_me[5] = _cuCsubf(m_me[5], _cuCmulf(sp1, m_me[2]));
            Real fDiv2 = F(1.0) / _sqrt(cuCabsSqf(m_me[3]) + cuCabsSqf(m_me[4]) + cuCabsSqf(m_me[5]));
            m_me[3] = _cuCmulf(m_me[3], fDiv2);
            m_me[4] = _cuCmulf(m_me[4], fDiv2);
            m_me[5] = _cuCmulf(m_me[5], fDiv2);

            m_me[6] = _cuConjf(_cuCsubf(_cuCmulf(m_me[1], m_me[5]), _cuCmulf(m_me[2], m_me[4])));
            m_me[7] = _cuConjf(_cuCsubf(_cuCmulf(m_me[2], m_me[3]), _cuCmulf(m_me[0], m_me[5])));
            m_me[8] = _cuConjf(_cuCsubf(_cuCmulf(m_me[0], m_me[4]), _cuCmulf(m_me[1], m_me[3])));
        }

        /**
        * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
        *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
        */
        __device__ __inline__ deviceSU3 Exp(const _Complex& a, UINT uiPrecision) const
        {
            deviceSU3 identity = makeSU3Id();
            deviceSU3 tmp;

            /**
            * tmp = U
            * loop1: tmp = 1+aU/N
            * loop2: tmp = 1+(aU/(N-1))(1+aU/N)
            * loop3: tmp = 1+(aU/(N-2))(1+aU/(N-1)(1+aU/N))
            * ...
            * loopN: tmp = 1+aU (1 + aU/2 (1+...))
            */
            for (UINT i = 0; i < uiPrecision; ++i)
            {
                Real exp_factor = F(1.0) / (uiPrecision - i);
                _Complex alpha = _cuCmulf(a, exp_factor);
                //aU/(N-i) = this x alpha
                deviceSU3 aUoN = MulCompC(alpha);
                if (0 == i)
                {
                    tmp = aUoN;
                }
                else
                {
                    tmp.Mul(aUoN);
                }
                tmp = identity.AddC(tmp);
            }
            tmp.Norm();
            return tmp;
        }

#pragma endregion

        union
        {
            deviceSU3Vector m_v[3];
            _Complex m_me[9];
        };
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
