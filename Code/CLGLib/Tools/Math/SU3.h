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

#define __LINE_MULDN(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[a]), right.m_me[d]), _cuCmulf(_cuConjf(m_me[b]), right.m_me[ee])), _cuCmulf(_cuConjf(m_me[c]), right.m_me[ff]))

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
            printf("=%1.7f%s%1.7fi, %1.7f%s%1.7fi, %1.7f%s%1.7fi;\n %1.7f%s%1.7fi, %1.7f%s%1.7fi, %1.7f%s%1.7fi;\n %1.7f%s%1.7fi, %1.7f%s%1.7fi, %1.7f%s%1.7fi;\n",
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
        * ret = i r_a T_a, r_a is random real number, T_a are generators
        * r_a T_a = r1 T1 + r2 T2 + ...
        *
        *     r3 +r8/sqrt3     r1-ir2        r4-ir5
        * =   r1+ir2       -r3+r8/sqrt3      r6-ir7
        *     r4+ir5           r6+ir7      -2r8/sqrt3
        *
        */
        __device__ __inline__ static deviceSU3 makeSU3RandomGenerator(UINT fatIndex)
        {
            Real r1 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r2 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r3 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r4 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r5 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r6 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r7 = _deviceRandomGaussFSqrt2(fatIndex);
            Real r8 = _deviceRandomGaussFSqrt2(fatIndex);

            deviceSU3 ret;
            //we directly generate i ra Ta instead of ra Ta
            ret.m_me[0] = _make_cuComplex(F(0.0), r8 * InvSqrt3 + r3);
            ret.m_me[1] = _make_cuComplex(r2, r1);
            ret.m_me[2] = _make_cuComplex(r5, r4);
            ret.m_me[3] = _make_cuComplex(-r2, r1);
            ret.m_me[4] = _make_cuComplex(F(0.0), r8 * InvSqrt3 - r3);
            ret.m_me[5] = _make_cuComplex(r7, r6);
            ret.m_me[6] = _make_cuComplex(-r5, r4);
            ret.m_me[7] = _make_cuComplex(-r7, r6);
            ret.m_me[8] = _make_cuComplex(F(0.0), -r8 * InvSqrt3_2);

            //ret.m_me[0] = _make_cuComplex(r3 + r8 * InvSqrt3, F(0.0));
            //ret.m_me[1] = _make_cuComplex(r1, -r2);
            //ret.m_me[2] = _make_cuComplex(r4, -r5);
            //ret.m_me[3] = _make_cuComplex(r1, r2);
            //ret.m_me[4] = _make_cuComplex(-r3 + r8 * InvSqrt3, F(0.0));
            //ret.m_me[5] = _make_cuComplex(r6, -r7);
            //ret.m_me[6] = _make_cuComplex(r4, r5);
            //ret.m_me[7] = _make_cuComplex(r6, r7);
            //ret.m_me[8] = _make_cuComplex(-r8 * InvSqrt3_2, F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3 makeSU3SumGenerator(Real fDivide)
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(fDivide * InvSqrt3 + fDivide, F(0.0));
            ret.m_me[1] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[2] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[3] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[4] = _make_cuComplex(fDivide * InvSqrt3 - fDivide, F(0.0));
            ret.m_me[5] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[6] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[7] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[8] = _make_cuComplex(-fDivide * InvSqrt3_2, F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3 makeSU3Generator(UINT uiGenerator)
        {
            deviceSU3 ret = makeSU3Zero();
            Real half = F(0.5);
            switch (uiGenerator)
            {
            case 0:
            {
                /**
                *     0     1     0
                * =   1     0     0
                *     0     0     0
                */
                ret.m_me[1] = _make_cuComplex(half, F(0.0));
                ret.m_me[3] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 1:
            {
                /**
                *     0     -I     0
                * =   I     0     0
                *     0     0     0
                */
                ret.m_me[1] = _make_cuComplex(F(0.0), -half);
                ret.m_me[3] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 2:
            {
                /**
                *     1     0     0
                * =   0    -1     0
                *     0     0     0
                */
                ret.m_me[0] = _make_cuComplex(half, F(0.0));
                ret.m_me[4] = _make_cuComplex(-half, F(0.0));
            }
            break;
            case 3:
            {
                /**
                *     0     0     1
                * =   0     0     0
                *     1     0     0
                */
                ret.m_me[2] = _make_cuComplex(half, F(0.0));
                ret.m_me[6] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 4:
            {
                /**
                *     0     0    -i
                * =   0     0     0
                *     i     0     0
                */
                ret.m_me[2] = _make_cuComplex(F(0.0), -half);
                ret.m_me[6] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 5:
            {
                /**
                *     0     0     0
                * =   0     0     1
                *     0     1     0
                */
                ret.m_me[5] = _make_cuComplex(half, F(0.0));
                ret.m_me[7] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 6:
            {
                /**
                *     0     0     0
                * =   0     0    -i
                *     0     i     0
                */
                ret.m_me[5] = _make_cuComplex(F(0.0), -half);
                ret.m_me[7] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 7:
            {
                /**
                *     1     0     0
                * =   0     1     0
                *     0     0     -2
                */
                ret.m_me[0] = _make_cuComplex(InvSqrt3 * half, F(0.0));
                ret.m_me[4] = _make_cuComplex(InvSqrt3 * half, F(0.0));
                ret.m_me[8] = _make_cuComplex(-InvSqrt3_2 * half, F(0.0));
            }
            break;
            default:
            case 8:
            {
                /**
                *     1     0     0
                * =   0     1     0
                *     0     0     -2
                */
                ret.m_me[1] = _make_cuComplex(half, -half);
                ret.m_me[3] = _make_cuComplex(half, half);
                ret.m_me[0] = _make_cuComplex(InvSqrt3 * half + half, F(0.0));
                ret.m_me[4] = _make_cuComplex(InvSqrt3 * half - half, F(0.0));
                ret.m_me[2] = _make_cuComplex(half, -half);
                ret.m_me[6] = _make_cuComplex(half, half);
                ret.m_me[5] = _make_cuComplex(half, -half);
                ret.m_me[7] = _make_cuComplex(half, half);
                ret.m_me[8] = _make_cuComplex(-InvSqrt3_2 * half, F(0.0));
            }
            break;
            }
            return ret;
        }

        /**
        * This is not right left^{dagger}
        * Let left = L1,L2,L3 right = R1, R2, R3
        * li = sum _{spin}Li
        * ri = sum _{spin}Ri
        * res = 
        *  l1^* r1, l2^* r1, l3^*r1
        *  l1^* r2, l2^* r2, l3^*r2
        *  l1^* r3, l2^* r3, l3^*r3
        */
        __device__ __inline__ static deviceSU3 makeSU3Contract(const deviceWilsonVectorSU3& left, const deviceWilsonVectorSU3& right)
        {
            deviceSU3 ret;
            ret.m_me[0] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[0]), right.m_d[0].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[0]), right.m_d[1].m_ve[0])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[0]), right.m_d[2].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[0]), right.m_d[3].m_ve[0])
                )
            );

            ret.m_me[1] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[1]), right.m_d[0].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[1]), right.m_d[1].m_ve[0])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[1]), right.m_d[2].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[1]), right.m_d[3].m_ve[0])
                )
            );

            ret.m_me[2] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[2]), right.m_d[0].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[2]), right.m_d[1].m_ve[0])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[2]), right.m_d[2].m_ve[0]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[2]), right.m_d[3].m_ve[0])
                )
            );

            ret.m_me[3] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[0]), right.m_d[0].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[0]), right.m_d[1].m_ve[1])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[0]), right.m_d[2].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[0]), right.m_d[3].m_ve[1])
                )
            );

            ret.m_me[4] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[1]), right.m_d[0].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[1]), right.m_d[1].m_ve[1])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[1]), right.m_d[2].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[1]), right.m_d[3].m_ve[1])
                )
            );

            ret.m_me[5] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[2]), right.m_d[0].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[2]), right.m_d[1].m_ve[1])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[2]), right.m_d[2].m_ve[1]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[2]), right.m_d[3].m_ve[1])
                )
            );

            ret.m_me[6] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[0]), right.m_d[0].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[0]), right.m_d[1].m_ve[2])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[0]), right.m_d[2].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[0]), right.m_d[3].m_ve[2])
                )
            );

            ret.m_me[7] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[1]), right.m_d[0].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[1]), right.m_d[1].m_ve[2])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[1]), right.m_d[2].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[1]), right.m_d[3].m_ve[2])
                )
            );

            ret.m_me[8] = _cuCaddf(
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[0].m_ve[2]), right.m_d[0].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[1].m_ve[2]), right.m_d[1].m_ve[2])
                ),
                _cuCaddf(
                    _cuCmulf(_cuConjf(left.m_d[2].m_ve[2]), right.m_d[2].m_ve[2]),
                    _cuCmulf(_cuConjf(left.m_d[3].m_ve[2]), right.m_d[3].m_ve[2])
                )
            );

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

        __device__ __inline__ void AddReal(Real right)
        {
            m_me[0] = cuCaddf_cr(m_me[0], right);
            m_me[4] = cuCaddf_cr(m_me[4], right);
            m_me[8] = cuCaddf_cr(m_me[8], right);
        }

        __device__ __inline__ void AddComp(const _Complex& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right);
            m_me[4] = _cuCaddf(m_me[4], right);
            m_me[8] = _cuCaddf(m_me[8], right);
        }

        __device__ __inline__ void SubReal(Real right)
        {
            m_me[0] = cuCsubf_cr(m_me[0], right);
            m_me[4] = cuCsubf_cr(m_me[4], right);
            m_me[8] = cuCsubf_cr(m_me[8], right);
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
            m_me[0] = cuCmulf_cr(m_me[0], right);
            m_me[1] = cuCmulf_cr(m_me[1], right);
            m_me[2] = cuCmulf_cr(m_me[2], right);
            m_me[3] = cuCmulf_cr(m_me[3], right);
            m_me[4] = cuCmulf_cr(m_me[4], right);
            m_me[5] = cuCmulf_cr(m_me[5], right);
            m_me[6] = cuCmulf_cr(m_me[6], right);
            m_me[7] = cuCmulf_cr(m_me[7], right);
            m_me[8] = cuCmulf_cr(m_me[8], right);
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

        __device__ __inline__ deviceSU3 MulDaggerC(const deviceSU3& right) const
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
            res[0] = __LINE_MULDN(0, 3, 6, 0, 3, 6);
            res[1] = __LINE_MULDN(0, 3, 6, 1, 4, 7);
            res[2] = __LINE_MULDN(0, 3, 6, 2, 5, 8);

            res[3] = __LINE_MULDN(1, 4, 7, 0, 3, 6);
            res[4] = __LINE_MULDN(1, 4, 7, 1, 4, 7);
            res[5] = __LINE_MULDN(1, 4, 7, 2, 5, 8);

            res[6] = __LINE_MULDN(2, 5, 8, 0, 3, 6);
            res[7] = __LINE_MULDN(2, 5, 8, 1, 4, 7);
            res[8] = __LINE_MULDN(2, 5, 8, 2, 5, 8);
            memcpy(m_me, res, sizeof(_Complex) * 9);
        }

        __device__ __inline__ deviceSU3 DaggerMulC(const deviceSU3& right) const
        {
            deviceSU3 ret;
            ret.m_me[0] = __LINE_MULDN(0, 3, 6, 0, 3, 6);
            ret.m_me[1] = __LINE_MULDN(0, 3, 6, 1, 4, 7);
            ret.m_me[2] = __LINE_MULDN(0, 3, 6, 2, 5, 8);

            ret.m_me[3] = __LINE_MULDN(1, 4, 7, 0, 3, 6);
            ret.m_me[4] = __LINE_MULDN(1, 4, 7, 1, 4, 7);
            ret.m_me[5] = __LINE_MULDN(1, 4, 7, 2, 5, 8);

            ret.m_me[6] = __LINE_MULDN(2, 5, 8, 0, 3, 6);
            ret.m_me[7] = __LINE_MULDN(2, 5, 8, 1, 4, 7);
            ret.m_me[8] = __LINE_MULDN(2, 5, 8, 2, 5, 8);
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

        __device__ __inline__ Real ImTr() const
        {
            return m_me[0].y + m_me[4].y + m_me[8].y;
        }

        __device__ __inline__ void Re()
        {
            m_me[0] = _make_cuComplex(m_me[0].x, F(0.0));
            m_me[1] = _make_cuComplex(m_me[1].x, F(0.0));
            m_me[2] = _make_cuComplex(m_me[2].x, F(0.0));
            m_me[3] = _make_cuComplex(m_me[3].x, F(0.0));
            m_me[4] = _make_cuComplex(m_me[4].x, F(0.0));
            m_me[5] = _make_cuComplex(m_me[5].x, F(0.0));
            m_me[6] = _make_cuComplex(m_me[6].x, F(0.0));
            m_me[7] = _make_cuComplex(m_me[7].x, F(0.0));
            m_me[8] = _make_cuComplex(m_me[8].x, F(0.0));
        }

        __device__ __inline__ deviceSU3 ReC() const { deviceSU3 ret(*this); ret.Re(); return ret; }

        __device__ __inline__ void Im()
        {
            m_me[0] = _make_cuComplex(m_me[0].y, F(0.0));
            m_me[1] = _make_cuComplex(m_me[1].y, F(0.0));
            m_me[2] = _make_cuComplex(m_me[2].y, F(0.0));
            m_me[3] = _make_cuComplex(m_me[3].y, F(0.0));
            m_me[4] = _make_cuComplex(m_me[4].y, F(0.0));
            m_me[5] = _make_cuComplex(m_me[5].y, F(0.0));
            m_me[6] = _make_cuComplex(m_me[6].y, F(0.0));
            m_me[7] = _make_cuComplex(m_me[7].y, F(0.0));
            m_me[8] = _make_cuComplex(m_me[8].y, F(0.0));
        }

        __device__ __inline__ deviceSU3 ImC() const { deviceSU3 ret(*this); ret.Im(); return ret; }

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


        __device__ __inline deviceSU3 Transpose() const 
        {

            deviceSU3 tmp;
            tmp.m_me[0] = m_me[0];
            tmp.m_me[1] = m_me[3];
            tmp.m_me[2] = m_me[6];

            tmp.m_me[3] = m_me[1];
            tmp.m_me[4] = m_me[4];
            tmp.m_me[5] = m_me[7];

            tmp.m_me[6] = m_me[2];
            tmp.m_me[7] = m_me[5];
            tmp.m_me[8] = m_me[8];
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

            Real fImgTr = __div(m_me[0].y + m_me[4].y + m_me[8].y, F(3.0));
            m_me[0] = _make_cuComplex(F(0.0), m_me[0].y - fImgTr);
            m_me[4] = _make_cuComplex(F(0.0), m_me[4].y - fImgTr);
            m_me[8] = _make_cuComplex(F(0.0), m_me[8].y - fImgTr);
        }

        /**
        * return -i(U-U^dagger) = ((-iU)+(-iU)dagger)
        */
        __device__ __inline__ deviceSU3 Im2C() const
        {
            deviceSU3 ret;
            //0 1 2
            //3 4 5
            //6 7 8
            // -i(x-y*)=i(y*-x)=i(yr-iyi-xr-ixi)=(yi+xi)+i(yr-xr)
            // new 1 = -i(1-3*)
            // new 2 = -i(2-6*)
            // new 5 = -i(5-7*)
            // (new 1)*= i(1*-3) = -i(3-1*) = new 3
            // new 0 = -i(0-0*) = 2Im(0)

            ret.m_me[1] = _make_cuComplex(m_me[1].y + m_me[3].y, m_me[3].x - m_me[1].x);
            ret.m_me[2] = _make_cuComplex(m_me[2].y + m_me[6].y, m_me[6].x - m_me[2].x);
            ret.m_me[5] = _make_cuComplex(m_me[5].y + m_me[7].y, m_me[7].x - m_me[5].x);
            ret.m_me[3] = _cuConjf(ret.m_me[1]);
            ret.m_me[6] = _cuConjf(ret.m_me[2]);
            ret.m_me[7] = _cuConjf(ret.m_me[5]);

            ret.m_me[0] = _make_cuComplex(F(2.0) * m_me[0].y, F(0.0));
            ret.m_me[4] = _make_cuComplex(F(2.0) * m_me[4].y, F(0.0));
            ret.m_me[8] = _make_cuComplex(F(2.0) * m_me[8].y, F(0.0));
            return ret;
        }

        /**
        * make any matrix to SU3
        */
        __device__ __inline__ void Norm()
        {
            Real fDiv1 = __div(F(1.0), _sqrt(__cuCabsSqf(m_me[0]) + __cuCabsSqf(m_me[1]) + __cuCabsSqf(m_me[2])));
            m_me[0] = cuCmulf_cr(m_me[0], fDiv1);
            m_me[1] = cuCmulf_cr(m_me[1], fDiv1);
            m_me[2] = cuCmulf_cr(m_me[2], fDiv1);

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
            Real fDiv2 = __div(F(1.0), _sqrt(__cuCabsSqf(m_me[3]) + __cuCabsSqf(m_me[4]) + __cuCabsSqf(m_me[5])));
            m_me[3] = cuCmulf_cr(m_me[3], fDiv2);
            m_me[4] = cuCmulf_cr(m_me[4], fDiv2);
            m_me[5] = cuCmulf_cr(m_me[5], fDiv2);

            m_me[6] = _cuConjf(_cuCsubf(_cuCmulf(m_me[1], m_me[5]), _cuCmulf(m_me[2], m_me[4])));
            m_me[7] = _cuConjf(_cuCsubf(_cuCmulf(m_me[2], m_me[3]), _cuCmulf(m_me[0], m_me[5])));
            m_me[8] = _cuConjf(_cuCsubf(_cuCmulf(m_me[0], m_me[4]), _cuCmulf(m_me[1], m_me[3])));
        }

        /**
        * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
        *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
        */
        __device__ __inline__ deviceSU3 Exp(const _Complex& a, UINT uiPrecision, UBOOL bNormed = FALSE) const
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
                Real exp_factor = __div(F(1.0), uiPrecision - i);
                _Complex alpha = cuCmulf_cr(a, exp_factor);
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
            if (bNormed)
            {
                tmp.Norm();
            }
            
            return tmp;
        }

        __device__ __inline__ deviceSU3 ExpReal(Real a, UINT uiPrecision, UBOOL bNormed = FALSE) const
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
                _Complex alpha = _make_cuComplex(__div(a, (uiPrecision - i)), F(0.0));
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
            if (bNormed)
            {
                tmp.Norm();
            }
            return tmp;
        }

        /**
        * For anti hermitian matrix only
        * see SU3Norm.nb
        */
        __device__ __inline__ deviceSU3 QuickExp(Real a) const
        {
            Real aOver12 = a * F(0.08333333333333333333333333333333);
            Real aOver6 = a * F(0.16666666666666666666666666666667);

            //y1 = a(Im(m11-m22))/12, y2 = a(Im(m11-m33))/12, y3= a(Im(m22-m33))/12
            Real y1 = aOver12 * (m_me[0].y - m_me[4].y);
            Real y2 = aOver12 * (m_me[0].y - m_me[8].y);
            Real y3 = aOver6 * (m_me[4].y - m_me[8].y);

            //l1 = y1^2 + (|m12|/4)^2
            //l2 = y2^2 + (|m13|/4)^2
            //l3 = y3^2 + (|m23|/2)^2
            Real a2 = a * a;
            Real aSqOver16 = a2 * F(0.0625);
            Real aSqOver4 = a2 * F(0.25);
            Real l1 = y1 * y1 + aSqOver16 * __cuCabsSqf(m_me[1]); 
            Real l2 = y2 * y2 + aSqOver16 * __cuCabsSqf(m_me[2]);
            Real l3 = y2 * y2 + aSqOver4 * __cuCabsSqf(m_me[5]);

            //dn1 = 1/ (1 + y1^2 +  (|m12|/4)^2)
            Real dn1 = __div(F(1.0), F(1.0) + l1);
            Real dn2 = __div(F(1.0), F(1.0) + l2);
            Real dn3 = __div(F(1.0), F(1.0) + l3);

            //dn1(1 - l1) = -l1 x dn1 + dn1
            Real c1r = dn1 - l1 * dn1;
            Real c2r = dn2 - l2 * dn2;
            Real c3r = dn3 - l3 * dn3;
            Real c1i = F(2.0) * y1 * dn1;
            Real c2i = F(2.0) * y2 * dn2;
            Real c3i = F(2.0) * y3 * dn3;

            Real halfa = F(0.5) * a;
            Real halfadn1 = halfa * dn1;
            Real halfadn2 = halfa * dn2;
            Real adn3 = a * dn3;
            Real a1r = halfadn1 * m_me[1].x;
            Real a2r = halfadn2 * m_me[2].x;
            Real a3r = adn3 * m_me[5].x;
            Real a1i = halfadn1 * m_me[1].y;
            Real a2i = halfadn2 * m_me[2].y;
            Real a3i = adn3 * m_me[5].y;

            deviceSU3 midMatrix;
            Real abs2a = a2r * a2r + a2i * a2i;
            Real abs2c = c2r * c2r - c2i * c2i;
            //q1 = c2i - c2r c3i - c2i c3r
            //q2 = c2r - c2i c3i + c2r c3r
            Real q1 = c2i - c2r * c3i - c2i * c3r;
            Real q2 = c2r - c2i * c3i + c2r * c3r; 
            Real c2ic2r = F(2.0) * c2i * c2r;
            //(abs2c-abs2a c3r)+I(2 c2i c2r+abs2a c3i)
            midMatrix.m_me[0] = _make_cuComplex(
                abs2c - abs2a * c3r,
                c2ic2r + abs2a * c3i
            );

            //-a2i a3i-a2r a3r+I (a2r a3i-a2i a3r)
            midMatrix.m_me[1] = _make_cuComplex(
                -a2i * a3i - a2r * a3r,
                 a2r * a3i - a2i * a3r
            );

            //(a2r q2-a2i q1)+(a2i q2+a2r q1)I
            Real a2rq1 = a2r * q1;
            Real a2iq1 = a2i * q1;
            Real a2rq2 = a2r * q2;
            Real a2iq2 = a2i * q2;
            midMatrix.m_me[2] = _make_cuComplex(
                a2rq2 - a2iq1,
                a2iq2 + a2rq1
            );

            midMatrix.m_me[3] = _make_cuComplex(midMatrix.m_me[1].x, -midMatrix.m_me[1].y);
            midMatrix.m_me[4] = _make_cuComplex(c3r, c3i);

            //a3i c2i+a3r c2r+I (-a3r c2i+a3i c2r)
            Real a3ic2i = a3i * c2i;
            Real a3rc2r = a3r * c2r;
            Real a3rc2i = a3r * c2i;
            Real a3ic2r = a3i * c2r;
            midMatrix.m_me[5] = _make_cuComplex(
                a3ic2i + a3rc2r,
                a3ic2r - a3rc2i
            );

            //-a2r q2-a2i q1+(a2i q2-a2r q1)I
            midMatrix.m_me[6] = _make_cuComplex(
                -a2rq2 - a2iq1,
                 a2iq2 - a2rq1
            );

            //a3i c2i-a3r c2r+I (a3r c2i+a3i c2r)
            midMatrix.m_me[7] = _make_cuComplex(
                a3ic2i - a3rc2r,
                a3rc2i + a3ic2r
            );

            //-abs2a-2 c2i c2r c3i+ c3r abs2c+I(- c3i abs2c-2 c2i c2r c3r)
            midMatrix.m_me[8] = _make_cuComplex(
                c3r * abs2c - abs2a - c3i * c2ic2r,
               -c3i * abs2c - c3r * c2ic2r
            );

            deviceSU3 u1Matrix = deviceSU3::makeSU3Id();
            u1Matrix.m_me[0] = _make_cuComplex(c1r, c1i);
            u1Matrix.m_me[4] = _make_cuComplex(c1r, -c1i);
            u1Matrix.m_me[1] = _make_cuComplex(a1r, a1i);
            u1Matrix.m_me[3] = _make_cuComplex(-a1r, a1i);

            midMatrix.Mul(u1Matrix);
            u1Matrix.Mul(midMatrix);

            return u1Matrix;
        }

#pragma endregion

        //union
        //{
        //    deviceSU3Vector m_v[3];
        //    _Complex m_me[9];
        //};
        _Complex m_me[16]; //Only the first 9 elememt is using, 16 is for padding and align
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
