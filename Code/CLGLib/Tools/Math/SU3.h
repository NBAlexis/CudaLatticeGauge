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

#define __LINE_MULARRAY(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(left[a], right[d]), _cuCmulf(left[b], right[ee])), _cuCmulf(left[c], right[ff]))

#define __LINE_MULND(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(m_me[a], _cuConjf(right.m_me[d])), _cuCmulf(m_me[b], _cuConjf(right.m_me[ee]))), _cuCmulf(m_me[c], _cuConjf(right.m_me[ff])))

#define __LINE_MULDN(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[a]), right.m_me[d]), _cuCmulf(_cuConjf(m_me[b]), right.m_me[ee])), _cuCmulf(_cuConjf(m_me[c]), right.m_me[ff]))

//self^dagger.self
#define __LINE_MULDN_ME(a, b, c, d, ee, ff) _cuCaddf(_cuCaddf(_cuCmulf(_cuConjf(m_me[a]), m_me[d]), _cuCmulf(_cuConjf(m_me[b]), m_me[ee])), _cuCmulf(_cuConjf(m_me[c]), m_me[ff]))


//x * a + y * b
//= (xr * ar - xi * ai) + (yr * br - yi * bi)
//+((xi * ar + xr * ai) + (yi * br + yr * bi))i
#define __TWO_COMPLEX_MULT_4(xr, xi, yr, yi, ar, ai, br, bi) \
    _make_cuComplex( ((xr) * (ar) - (xi) * (ai)) + ((yr) * (br) - (yi) * (bi)), \
                     ((xi) * (ar) + (xr) * (ai)) + ((yi) * (br) + (yr) * (bi)) );



// C = (A - a) (B - b)
//tmpM1.m_me[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], m_me[0]), _cuCmulf(m_me[1], m_me[3])), _cuCmulf(m_me[2], m_me[6]));
#define __LINE_MUL_POLY(res, a, b, c, d, e, f) res = _cuCaddf(_cuCaddf(_cuCmulf(a, d), _cuCmulf(b, e)), _cuCmulf(c, f));


#if _CLG_DOUBLEFLOAT
#define __SU3MATRIX_ALIGN 256
#else
#define __SU3MATRIX_ALIGN 128
#endif

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#if _CLG_DOUBLEFLOAT
    //alignas 256 will crash ptax, so we manully align
    struct /*alignas(__SU3MATRIX_ALIGN)*/ deviceSU3
#else
    struct /*alignas(__SU3MATRIX_ALIGN)*/ deviceSU3
#endif
    {
        __device__ deviceSU3()
        {

        }

        __device__ deviceSU3(const deviceSU3& other)
        {
            memcpy(m_me, other.m_me, sizeof(CLGComplex) * 9);
        }

        __device__ deviceSU3(const CLGComplex* __restrict__ other)
        {
            memcpy(m_me, other, sizeof(CLGComplex) * 9);
        }

        __device__ void DebugPrint() const
        {
            printf("={{%1.7f%s%1.7f I, %1.7f%s%1.7f I, %1.7f%s%1.7f I},\n {%1.7f%s%1.7f I, %1.7f%s%1.7f I, %1.7f%s%1.7f I},\n {%1.7f%s%1.7f I, %1.7f%s%1.7f I, %1.7f%s%1.7f I}};\n",
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
            const Real r1 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r2 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r3 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r4 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r5 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r6 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r7 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r8 = _deviceRandomGaussFSqrt2(fatIndex);

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
            const Real half = F(0.5);
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
         *        m11 I  m12   m13
         * ret =  -m12*  m22 I m23
         *        -m13* -m23*  -(m11+m22)I
         */
        __device__ __inline__ static deviceSU3 makeSU3TA(
            const CLGComplex& m12, 
            const CLGComplex& m13, 
            const CLGComplex& m23,
            Real m11, Real m22)
        {
            deviceSU3 ret;
            ret.m_me[0] = _make_cuComplex(F(0.0), m11);;
            ret.m_me[1] = m12;
            ret.m_me[2] = m13;
            ret.m_me[3] = _make_cuComplex(-m12.x, m12.y);
            ret.m_me[4] = _make_cuComplex(F(0.0), m22);
            ret.m_me[5] = m23;
            ret.m_me[6] = _make_cuComplex(-m13.x, m13.y);
            ret.m_me[7] = _make_cuComplex(-m23.x, m23.y);
            ret.m_me[8] = _make_cuComplex(F(0.0), - m11 - m22);;
            return ret;
        }

        /**
        * This is right left^{dagger} for column SU3 vector with spin index summed
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

        /**
         * Another version of makeSU3Contract without spin sum
         * if v = {{},{},{}}, column vector
         * This is right.left^+
         */
        __device__ __inline__ static deviceSU3 makeSU3ContractV(const deviceSU3Vector& left, const deviceSU3Vector& right)
        {
            deviceSU3 ret;
            ret.m_me[0] = _cuCmulf(_cuConjf(left.m_ve[0]), right.m_ve[0]);
            ret.m_me[1] = _cuCmulf(_cuConjf(left.m_ve[1]), right.m_ve[0]);
            ret.m_me[2] = _cuCmulf(_cuConjf(left.m_ve[2]), right.m_ve[0]);
            ret.m_me[3] = _cuCmulf(_cuConjf(left.m_ve[0]), right.m_ve[1]);
            ret.m_me[4] = _cuCmulf(_cuConjf(left.m_ve[1]), right.m_ve[1]);
            ret.m_me[5] = _cuCmulf(_cuConjf(left.m_ve[2]), right.m_ve[1]);
            ret.m_me[6] = _cuCmulf(_cuConjf(left.m_ve[0]), right.m_ve[2]);
            ret.m_me[7] = _cuCmulf(_cuConjf(left.m_ve[1]), right.m_ve[2]);
            ret.m_me[8] = _cuCmulf(_cuConjf(left.m_ve[2]), right.m_ve[2]);
            
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

        __device__ __inline__ void AddDagger(const deviceSU3& right)
        {
            m_me[0] = _cuCaddf(m_me[0], _cuConjf(right.m_me[0]));
            m_me[1] = _cuCaddf(m_me[1], _cuConjf(right.m_me[3]));
            m_me[2] = _cuCaddf(m_me[2], _cuConjf(right.m_me[6]));
            m_me[3] = _cuCaddf(m_me[3], _cuConjf(right.m_me[1]));
            m_me[4] = _cuCaddf(m_me[4], _cuConjf(right.m_me[4]));
            m_me[5] = _cuCaddf(m_me[5], _cuConjf(right.m_me[7]));
            m_me[6] = _cuCaddf(m_me[6], _cuConjf(right.m_me[2]));
            m_me[7] = _cuCaddf(m_me[7], _cuConjf(right.m_me[5]));
            m_me[8] = _cuCaddf(m_me[8], _cuConjf(right.m_me[8]));
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

        __device__ __inline__ void SubDagger(const deviceSU3& right)
        {
            m_me[0] = _cuCsubf(m_me[0], _cuConjf(right.m_me[0]));
            m_me[1] = _cuCsubf(m_me[1], _cuConjf(right.m_me[3]));
            m_me[2] = _cuCsubf(m_me[2], _cuConjf(right.m_me[6]));
            m_me[3] = _cuCsubf(m_me[3], _cuConjf(right.m_me[1]));
            m_me[4] = _cuCsubf(m_me[4], _cuConjf(right.m_me[4]));
            m_me[5] = _cuCsubf(m_me[5], _cuConjf(right.m_me[7]));
            m_me[6] = _cuCsubf(m_me[6], _cuConjf(right.m_me[2]));
            m_me[7] = _cuCsubf(m_me[7], _cuConjf(right.m_me[5]));
            m_me[8] = _cuCsubf(m_me[8], _cuConjf(right.m_me[8]));
        }

        __device__ __inline__ void AddReal(Real right)
        {
            m_me[0] = cuCaddf_cr(m_me[0], right);
            m_me[4] = cuCaddf_cr(m_me[4], right);
            m_me[8] = cuCaddf_cr(m_me[8], right);
        }

        __device__ __inline__ void AddComp(const CLGComplex& right)
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

        __device__ __inline__ void SubComp(const CLGComplex& right)
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
            CLGComplex res[9];
            res[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
            res[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
            res[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

            res[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
            res[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
            res[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

            res[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
            res[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
            res[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
            memcpy(m_me, res, sizeof(CLGComplex) * 9);
        }

        /**
        * low reg count version of multiply
        */
        __device__ __inline__ static void ArrayMul(CLGComplex* left, const CLGComplex* __restrict__ right)
        {
            CLGComplex tmp[3];
            tmp[0] = __LINE_MULARRAY(0, 1, 2, 0, 3, 6);
            tmp[1] = __LINE_MULARRAY(0, 1, 2, 1, 4, 7);
            tmp[2] = __LINE_MULARRAY(0, 1, 2, 2, 5, 8);
            memcpy(left, tmp, sizeof(CLGComplex) * 3);

            tmp[0] = __LINE_MULARRAY(3, 4, 5, 0, 3, 6);
            tmp[1] = __LINE_MULARRAY(3, 4, 5, 1, 4, 7);
            tmp[2] = __LINE_MULARRAY(3, 4, 5, 2, 5, 8);
            memcpy(left + 3, tmp, sizeof(CLGComplex) * 3);

            tmp[0] = __LINE_MULARRAY(6, 7, 8, 0, 3, 6);
            tmp[1] = __LINE_MULARRAY(6, 7, 8, 1, 4, 7);
            tmp[2] = __LINE_MULARRAY(6, 7, 8, 2, 5, 8);
            memcpy(left + 6, tmp, sizeof(CLGComplex) * 3);
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

        __device__ __inline__ void MulComp(const CLGComplex& right)
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

        __device__ __inline__ void DivComp(const CLGComplex& right)
        {
            m_me[0] = _cuCdivf(m_me[0], right);
            m_me[1] = _cuCdivf(m_me[1], right);
            m_me[2] = _cuCdivf(m_me[2], right);
            m_me[3] = _cuCdivf(m_me[3], right);
            m_me[4] = _cuCdivf(m_me[4], right);
            m_me[5] = _cuCdivf(m_me[5], right);
            m_me[6] = _cuCdivf(m_me[6], right);
            m_me[7] = _cuCdivf(m_me[7], right);
            m_me[8] = _cuCdivf(m_me[8], right);
        }

        __device__ __inline__ void MulDagger(const deviceSU3& right)
        {
            CLGComplex res[9];
            res[0] = __LINE_MULND(0, 1, 2, 0, 1, 2);
            res[1] = __LINE_MULND(0, 1, 2, 3, 4, 5);
            res[2] = __LINE_MULND(0, 1, 2, 6, 7, 8);

            res[3] = __LINE_MULND(3, 4, 5, 0, 1, 2);
            res[4] = __LINE_MULND(3, 4, 5, 3, 4, 5);
            res[5] = __LINE_MULND(3, 4, 5, 6, 7, 8);

            res[6] = __LINE_MULND(6, 7, 8, 0, 1, 2);
            res[7] = __LINE_MULND(6, 7, 8, 3, 4, 5);
            res[8] = __LINE_MULND(6, 7, 8, 6, 7, 8);
            memcpy(m_me, res, sizeof(CLGComplex) * 9);
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
            CLGComplex res[9];
            res[0] = __LINE_MULDN(0, 3, 6, 0, 3, 6);
            res[1] = __LINE_MULDN(0, 3, 6, 1, 4, 7);
            res[2] = __LINE_MULDN(0, 3, 6, 2, 5, 8);

            res[3] = __LINE_MULDN(1, 4, 7, 0, 3, 6);
            res[4] = __LINE_MULDN(1, 4, 7, 1, 4, 7);
            res[5] = __LINE_MULDN(1, 4, 7, 2, 5, 8);

            res[6] = __LINE_MULDN(2, 5, 8, 0, 3, 6);
            res[7] = __LINE_MULDN(2, 5, 8, 1, 4, 7);
            res[8] = __LINE_MULDN(2, 5, 8, 2, 5, 8);
            memcpy(m_me, res, sizeof(CLGComplex) * 9);
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
        __device__ __inline__ deviceSU3 AddCompC(const CLGComplex& right) const { deviceSU3 ret(*this); ret.AddComp(right); return ret; }
        __device__ __inline__ deviceSU3 AddRealC(const Real& right) const { deviceSU3 ret(*this); ret.AddReal(right); return ret; }

        __device__ __inline__ deviceSU3 SubC(const deviceSU3& right) const { deviceSU3 ret(*this); ret.Sub(right); return ret; }
        __device__ __inline__ deviceSU3 SubCompC(const CLGComplex& right) const { deviceSU3 ret(*this); ret.SubComp(right); return ret; }
        __device__ __inline__ deviceSU3 SubRealC(const Real& right) const { deviceSU3 ret(*this); ret.SubReal(right); return ret; }

        __device__ __inline__ deviceSU3 MulCompC(const CLGComplex& right) const { deviceSU3 ret(*this); ret.MulComp(right); return ret; }
        __device__ __inline__ deviceSU3 MulRealC(const Real& right) const { deviceSU3 ret(*this); ret.MulReal(right); return ret; }
        __device__ __inline__ deviceSU3 DivCompC(const CLGComplex& right) const { deviceSU3 ret(*this); ret.DivComp(right); return ret; }

        //the reduce regcount version
        __device__ __inline__ void MulRealCArray(CLGComplex* res, const Real& right) const 
        { 
            res[0] = _make_cuComplex(m_me[0].x * right, m_me[0].y * right);
            res[1] = _make_cuComplex(m_me[1].x * right, m_me[1].y * right);
            res[2] = _make_cuComplex(m_me[2].x * right, m_me[2].y * right);
            res[3] = _make_cuComplex(m_me[3].x * right, m_me[3].y * right);
            res[4] = _make_cuComplex(m_me[4].x * right, m_me[4].y * right);
            res[5] = _make_cuComplex(m_me[5].x * right, m_me[5].y * right);
            res[6] = _make_cuComplex(m_me[6].x * right, m_me[6].y * right);
            res[7] = _make_cuComplex(m_me[7].x * right, m_me[7].y * right);
            res[8] = _make_cuComplex(m_me[8].x * right, m_me[8].y * right);
        }

#pragma endregion operators

#pragma region useful functions

        /**
        * ret = ||U||, where U is a d=3 matrix
        * we do NOT need ||SU3|| because ||SU3||=1
        */
        __device__ static __inline__ CLGComplex Determinent(const CLGComplex* u)
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

        __device__ __inline__ CLGComplex Tr() const
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
            CLGComplex res[9];
            res[0] = _cuConjf(m_me[0]);
            res[1] = _cuConjf(m_me[3]);
            res[2] = _cuConjf(m_me[6]);

            res[3] = _cuConjf(m_me[1]);
            res[4] = _cuConjf(m_me[4]);
            res[5] = _cuConjf(m_me[7]);

            res[6] = _cuConjf(m_me[2]);
            res[7] = _cuConjf(m_me[5]);
            res[8] = _cuConjf(m_me[8]);

            memcpy(m_me, res, sizeof(CLGComplex) * 9);
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
            const CLGComplex new1 = _cuCsubf(m_me[1], _cuConjf(m_me[3]));
            const CLGComplex new2 = _cuCsubf(m_me[2], _cuConjf(m_me[6]));
            const CLGComplex new5 = _cuCsubf(m_me[5], _cuConjf(m_me[7]));
            m_me[1] = _make_cuComplex(F(0.5) * _cuCrealf(new1), F(0.5) * _cuCimagf(new1));
            m_me[3] = _make_cuComplex(-_cuCrealf(m_me[1]), _cuCimagf(m_me[1]));
            m_me[2] = _make_cuComplex(F(0.5) * _cuCrealf(new2), F(0.5) * _cuCimagf(new2));
            m_me[6] = _make_cuComplex(-_cuCrealf(m_me[2]), _cuCimagf(m_me[2]));
            m_me[5] = _make_cuComplex(F(0.5) * _cuCrealf(new5), F(0.5) * _cuCimagf(new5));
            m_me[7] = _make_cuComplex(-_cuCrealf(m_me[5]), _cuCimagf(m_me[5]));

            const Real fImgTr = __div(m_me[0].y + m_me[4].y + m_me[8].y, F(3.0));
            m_me[0] = _make_cuComplex(F(0.0), m_me[0].y - fImgTr);
            m_me[4] = _make_cuComplex(F(0.0), m_me[4].y - fImgTr);
            m_me[8] = _make_cuComplex(F(0.0), m_me[8].y - fImgTr);
        }

        /**
        * Return Tr[Im[a].Im[b]] 
        * It is used in topological charge
        */
        __device__ __inline__ static Real TrIm(const deviceSU3& a, const deviceSU3&b)
        {
            //0,1,2 * 0,3,6
            Real ret = a.m_me[0].y * b.m_me[0].y + a.m_me[1].y * b.m_me[3].y + a.m_me[2].y * b.m_me[6].y;
            //3,4,5 * 1,4,7
            ret += a.m_me[3].y * b.m_me[1].y + a.m_me[4].y * b.m_me[4].y + a.m_me[5].y * b.m_me[7].y;
            //6,7,8 * 2,5,8
            ret += a.m_me[6].y * b.m_me[2].y + a.m_me[7].y * b.m_me[5].y + a.m_me[8].y * b.m_me[8].y;

            return ret;
        }

        /**
        * res = U - U^dagger
        * It is like a matrix Im(M)
        * res = 2i Im(M), so is called iIm2
        */
        __device__ __inline__ void iIm2()
        {
            //0 1 2
            //3 4 5
            //6 7 8

            //new [1] = [1] - conj([3])
            //new [3] = [3] - conj([1]) = -conj(new [1])
            const CLGComplex new1 = _cuCsubf(m_me[1], _cuConjf(m_me[3]));
            const CLGComplex new2 = _cuCsubf(m_me[2], _cuConjf(m_me[6]));
            const CLGComplex new5 = _cuCsubf(m_me[5], _cuConjf(m_me[7]));
            m_me[1] = _make_cuComplex(_cuCrealf(new1), _cuCimagf(new1));
            m_me[3] = _make_cuComplex(-_cuCrealf(m_me[1]), _cuCimagf(m_me[1]));
            m_me[2] = _make_cuComplex(_cuCrealf(new2), _cuCimagf(new2));
            m_me[6] = _make_cuComplex(-_cuCrealf(m_me[2]), _cuCimagf(m_me[2]));
            m_me[5] = _make_cuComplex(_cuCrealf(new5), _cuCimagf(new5));
            m_me[7] = _make_cuComplex(-_cuCrealf(m_me[5]), _cuCimagf(m_me[5]));

            m_me[0] = _make_cuComplex(F(0.0), F(2.0) * m_me[0].y);
            m_me[4] = _make_cuComplex(F(0.0), F(2.0) * m_me[4].y);
            m_me[8] = _make_cuComplex(F(0.0), F(2.0) * m_me[8].y);
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
            const Real fDiv1 = __div(F(1.0), _sqrt(__cuCabsSqf(m_me[0]) + __cuCabsSqf(m_me[1]) + __cuCabsSqf(m_me[2])));
            m_me[0] = cuCmulf_cr(m_me[0], fDiv1);
            m_me[1] = cuCmulf_cr(m_me[1], fDiv1);
            m_me[2] = cuCmulf_cr(m_me[2], fDiv1);

            //it is the name in Bridge++
            const CLGComplex sp1 = _cuConjf(
                _cuCaddf(_cuCaddf(
                    _cuCmulf(m_me[0], _cuConjf(m_me[3]))
                    , _cuCmulf(m_me[1], _cuConjf(m_me[4])))
                    , _cuCmulf(m_me[2], _cuConjf(m_me[5])))
            );

            m_me[3] = _cuCsubf(m_me[3], _cuCmulf(sp1, m_me[0]));
            m_me[4] = _cuCsubf(m_me[4], _cuCmulf(sp1, m_me[1]));
            m_me[5] = _cuCsubf(m_me[5], _cuCmulf(sp1, m_me[2]));
            const Real fDiv2 = __div(F(1.0), _sqrt(__cuCabsSqf(m_me[3]) + __cuCabsSqf(m_me[4]) + __cuCabsSqf(m_me[5])));
            m_me[3] = cuCmulf_cr(m_me[3], fDiv2);
            m_me[4] = cuCmulf_cr(m_me[4], fDiv2);
            m_me[5] = cuCmulf_cr(m_me[5], fDiv2);

            m_me[6] = _cuConjf(_cuCsubf(_cuCmulf(m_me[1], m_me[5]), _cuCmulf(m_me[2], m_me[4])));
            m_me[7] = _cuConjf(_cuCsubf(_cuCmulf(m_me[2], m_me[3]), _cuCmulf(m_me[0], m_me[5])));
            m_me[8] = _cuConjf(_cuCsubf(_cuCmulf(m_me[0], m_me[4]), _cuCmulf(m_me[1], m_me[3])));
        }

        /**
         * HYP projection
         * When Det[U] is large, it need large ite...
         */
        __device__ __inline__ void Proj(BYTE ite = 4)
        {
            //tr(me^+ me) = m0^* m0 + m1^* m1 + m2^* m2 + ...
            Real fDiv =
                __cuCabsSqf(m_me[0])
                + __cuCabsSqf(m_me[1])
                + __cuCabsSqf(m_me[2])
                + __cuCabsSqf(m_me[3])
                + __cuCabsSqf(m_me[4])
                + __cuCabsSqf(m_me[5])
                + __cuCabsSqf(m_me[6])
                + __cuCabsSqf(m_me[7])
                + __cuCabsSqf(m_me[8]);
            //1 / sqrt( tr(me^+ me) / 3 )
            fDiv = __rcp(_sqrt(fDiv * F(0.3333333333333333)));
            m_me[0] = cuCmulf_cr(m_me[0], fDiv);
            m_me[1] = cuCmulf_cr(m_me[1], fDiv);
            m_me[2] = cuCmulf_cr(m_me[2], fDiv);
            m_me[3] = cuCmulf_cr(m_me[3], fDiv);
            m_me[4] = cuCmulf_cr(m_me[4], fDiv);
            m_me[5] = cuCmulf_cr(m_me[5], fDiv);
            m_me[6] = cuCmulf_cr(m_me[6], fDiv);
            m_me[7] = cuCmulf_cr(m_me[7], fDiv);
            m_me[8] = cuCmulf_cr(m_me[8], fDiv);

            CLGComplex x[9];
            CLGComplex tmp[3];
            for (BYTE byIt = 0; byIt < ite; ++byIt)
            {
                //x = (-1/2) me^+ me
                x[0] = cuCmulf_cr(__LINE_MULDN_ME(0, 3, 6, 0, 3, 6), F(-0.5));
                x[1] = cuCmulf_cr(__LINE_MULDN_ME(0, 3, 6, 1, 4, 7), F(-0.5));
                x[2] = cuCmulf_cr(__LINE_MULDN_ME(0, 3, 6, 2, 5, 8), F(-0.5));
                x[3] = cuCmulf_cr(__LINE_MULDN_ME(1, 4, 7, 0, 3, 6), F(-0.5));
                x[4] = cuCmulf_cr(__LINE_MULDN_ME(1, 4, 7, 1, 4, 7), F(-0.5));
                x[5] = cuCmulf_cr(__LINE_MULDN_ME(1, 4, 7, 2, 5, 8), F(-0.5));
                x[6] = cuCmulf_cr(__LINE_MULDN_ME(2, 5, 8, 0, 3, 6), F(-0.5));
                x[7] = cuCmulf_cr(__LINE_MULDN_ME(2, 5, 8, 1, 4, 7), F(-0.5));
                x[8] = cuCmulf_cr(__LINE_MULDN_ME(2, 5, 8, 2, 5, 8), F(-0.5));

                //x += (3/2)
                x[0].x = x[0].x + F(1.5);
                x[4].x = x[4].x + F(1.5);
                x[8].x = x[8].x + F(1.5);

                //x = me.x
                tmp[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], x[0]), _cuCmulf(m_me[1], x[3])), _cuCmulf(m_me[2], x[6]));
                tmp[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], x[1]), _cuCmulf(m_me[1], x[4])), _cuCmulf(m_me[2], x[7]));
                tmp[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[0], x[2]), _cuCmulf(m_me[1], x[5])), _cuCmulf(m_me[2], x[8]));
                //we do not need m_me[0,1,2] anymore
                memcpy(m_me, tmp, sizeof(CLGComplex) * 3);
                tmp[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], x[0]), _cuCmulf(m_me[4], x[3])), _cuCmulf(m_me[5], x[6]));
                tmp[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], x[1]), _cuCmulf(m_me[4], x[4])), _cuCmulf(m_me[5], x[7]));
                tmp[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[3], x[2]), _cuCmulf(m_me[4], x[5])), _cuCmulf(m_me[5], x[8]));
                memcpy(m_me + 3, tmp, sizeof(CLGComplex) * 3);
                tmp[0] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[6], x[0]), _cuCmulf(m_me[7], x[3])), _cuCmulf(m_me[8], x[6]));
                tmp[1] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[6], x[1]), _cuCmulf(m_me[7], x[4])), _cuCmulf(m_me[8], x[7]));
                tmp[2] = _cuCaddf(_cuCaddf(_cuCmulf(m_me[6], x[2]), _cuCmulf(m_me[7], x[5])), _cuCmulf(m_me[8], x[8]));
                memcpy(m_me + 6, tmp, sizeof(CLGComplex) * 3);

                //me = x
                //coef = det(me)
                CLGComplex coef = Determinent(m_me);
                //coef = 1 - (i/3)Im(coef)
                coef = _make_cuComplex(F(1.0), -F(0.3333333333333333) * coef.y);
                //me = coef me
                m_me[0] = _cuCmulf(coef, m_me[0]);
                m_me[1] = _cuCmulf(coef, m_me[1]);
                m_me[2] = _cuCmulf(coef, m_me[2]);
                m_me[3] = _cuCmulf(coef, m_me[3]);
                m_me[4] = _cuCmulf(coef, m_me[4]);
                m_me[5] = _cuCmulf(coef, m_me[5]);
                m_me[6] = _cuCmulf(coef, m_me[6]);
                m_me[7] = _cuCmulf(coef, m_me[7]);
                m_me[8] = _cuCmulf(coef, m_me[8]);
            }
        }

        /**
         * Cabbibo-Marinari Projection
         * Only one iteration can make it SU3, however, 6 iteration makes it large trace
         * NOTE: It has been tested using Mathematica, the trace is not increasing...
         *
         * Note for iteration:
         * u1 = m.cmproj()
         * m1 = u1.m
         * u2 = m1.cmproj().u1
         * m3 = u2.m
         * u3 = m3.cmproj().u2
         * ...
         * retr(u3.m)>retr(u2.m)>retr(u1.m), where u1,u2,u3 are su3, m is random matrix
         * finally, u = dagger(u3) so that retr(u3^+ . m) or retr(u3.m^+) is maximized
         */
        __device__ __inline__ void CabbiboMarinariProj(/*BYTE ite = 1*/)
        {
            //for (BYTE iteration = 0; iteration < ite; ++iteration)
            //{
                CLGComplex a11 = _cuCaddf(_cuConjf(m_me[0]), m_me[4]);
                CLGComplex b11 = _cuCaddf(_cuConjf(m_me[0]), m_me[8]);
                CLGComplex c22 = _cuCaddf(_cuConjf(m_me[4]), m_me[8]);
                CLGComplex a12 = _cuCsubf(_cuConjf(m_me[3]), m_me[1]);
                CLGComplex b13 = _cuCsubf(_cuConjf(m_me[6]), m_me[2]);
                CLGComplex c23 = _cuCsubf(_cuConjf(m_me[7]), m_me[5]);
                //CLGComplex a12 = _cuCsubf(_cuConjf(m_me[1]), m_me[3]);
                //CLGComplex b13 = _cuCsubf(_cuConjf(m_me[2]), m_me[6]);
                //CLGComplex c23 = _cuCsubf(_cuConjf(m_me[5]), m_me[7]);
                Real fNorm = __rcp(_sqrt(__cuCabsSqf(a11) + __cuCabsSqf(a12)));
                a11 = cuCmulf_cr(a11, fNorm);
                a12 = cuCmulf_cr(a12, fNorm);
                fNorm = __rcp(_sqrt(__cuCabsSqf(b11) + __cuCabsSqf(b13)));
                b11 = cuCmulf_cr(b11, fNorm);
                b13 = cuCmulf_cr(b13, fNorm);
                fNorm = __rcp(_sqrt(__cuCabsSqf(c22) + __cuCabsSqf(c23)));
                c22 = cuCmulf_cr(c22, fNorm);
                c23 = cuCmulf_cr(c23, fNorm);

                /**
                 * ({
                {a11 b11,
                a12 c22 - a11 b13 Conjugate[c23], 
                a11 b13 Conjugate[c22] + a12 c23},
                {-b11 Conjugate[a12], 
                c22 Conjugate[a11] + b13 Conjugate[a12] Conjugate[c23], 
                c23 Conjugate[a11] - b13 Conjugate[a12] Conjugate[c22]},
                {-Conjugate[b13],
                -Conjugate[b11] Conjugate[c23], 
                Conjugate[b11] Conjugate[c22]}
                })
                 */
                m_me[0] = _cuCmulf(a11, b11);
                m_me[1] = _cuCsubf(_cuCmulf(a12, c22), _cuCmulf(_cuCmulf(a11, b13), _cuConjf(c23)));
                m_me[2] = _cuCaddf(_cuCmulf(_cuCmulf(a11, b13), _cuConjf(c22)), _cuCmulf(a12, c23));

                m_me[3].x = -b11.x * a12.x - b11.y * a12.y;
                m_me[3].y = a12.y * b11.x - a12.x * b11.y;
                m_me[4] = _cuCaddf(_cuCmulf(c22, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c23)));
                m_me[5] = _cuCsubf(_cuCmulf(c23, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c22)));

                m_me[6].x = -b13.x; m_me[6].y = b13.y;
                m_me[7].x = b11.y * c23.y - b11.x * c23.x;
                m_me[7].y = b11.x * c23.y + b11.y * c23.x;
                m_me[8].x = b11.x * c22.x - b11.y * c22.y;
                m_me[8].y = - b11.x * c22.y - b11.y * c22.x;
            //}
        }

        /**
        * U = U + 1
        */
        __device__ __inline__ void AddId()
        {
            m_me[0].x += F(1.0);
            m_me[4].x += F(1.0);
            m_me[8].x += F(1.0);
        }

        /**
        * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
        *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
        */
        __device__ __inline__ deviceSU3 Exp(const CLGComplex& a, BYTE uiPrecision) const
        {
            deviceSU3 tmp;

            /**
            * tmp = U
            * loop1: tmp = 1+aU/N
            * loop2: tmp = 1+(aU/(N-1))(1+aU/N)
            * loop3: tmp = 1+(aU/(N-2))(1+aU/(N-1)(1+aU/N))
            * ...
            * loopN: tmp = 1+aU (1 + aU/2 (1+...))
            */
            for (BYTE i = 0; i < uiPrecision; ++i)
            {
                const Real exp_factor = __rcp(uiPrecision - i);
                CLGComplex alpha = cuCmulf_cr(a, exp_factor);
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
                tmp.AddId();
            }
            //if (bNormed)
            //{
            //    tmp.Norm();
            //}
            
            return tmp;
        }

        __device__ __inline__ deviceSU3 ExpReal(Real a, BYTE uiPrecision) const
        {
            CLGComplex tmp[9];

            /**
            * tmp = U
            * loop1: tmp = 1+aU/N
            * loop2: tmp = 1+(aU/(N-1))(1+aU/N)
            * loop3: tmp = 1+(aU/(N-2))(1+aU/(N-1)(1+aU/N))
            * ...
            * loopN: tmp = 1+aU (1 + aU/2 (1+...))
            */
            for (BYTE i = 0; i < uiPrecision; ++i)
            {
                //CLGComplex alpha = _make_cuComplex(__div(a, (uiPrecision - i)), F(0.0));
                //aU/(N-i) = this x alpha
                //deviceSU3 aUoN = MulCompC(_make_cuComplex(__div(a, (uiPrecision - i)), F(0.0)));
                if (0 == i)
                {
                    MulRealCArray(tmp, __div(a, (uiPrecision - i)));
                }
                else
                {
                    CLGComplex tmp2[9];
                    MulRealCArray(tmp2, __div(a, (uiPrecision - i)));
                    ArrayMul(tmp, tmp2);
                }
                tmp[0].x += F(1.0);
                tmp[4].x += F(1.0);
                tmp[8].x += F(1.0);
            }
            return deviceSU3(tmp);
        }

        //for release, let the compliler do this
#if _CLG_DEBUG
#define aOver12 tmpReal[0]
#define aOver6 tmpReal[1]
#define y1 tmpReal[2]
#define y2 tmpReal[0] //aOver12 not used
#define y3 tmpReal[1] //aOver6 not used
#define a2 tmpReal[3]
#define aSqOver16 tmpReal[4]
#define aSqOver4 tmpReal[3] //a2 not used
#define l1 tmpReal[5]
#define l2 tmpReal[4] //aSqOver16 not used
#define l3 tmpReal[3] //aSqOver4 not used
#define dn1 tmpReal[6]
#define dn2 tmpReal[7]
#define dn3 tmpReal[8] //now y1(2),y2(0),y3(1),l1(5),l2(4),l3(3) all not used
#define halfa tmpReal[0]
#define halfadn1 tmpReal[1]
#define halfadn2 tmpReal[0] //halfa not used
#define adn3 tmpReal[2]

#define abs2a tmpReal[0]
#define abs2c tmpReal[1]
#define c2ic2r tmpReal[2]

#define q1 tmpReal[0]
#define q2 tmpReal[1]
#define a2rq1 tmpReal[2]
#define a2iq1 tmpReal[0]
#define a2rq2 tmpReal[3]
#define a2iq2 tmpReal[1]

#define a3ic2i tmpReal[0]
#define a3rc2r tmpReal[1]
#define a3rc2i tmpReal[2]
#define a3ic2r tmpReal[3]
#endif
        /**
        * For anti hermitian matrix only
        * see SU3Norm.nb
        */
        __device__ __inline__ deviceSU3 QuickExp(Real a) const
        {
#if _CLG_DEBUG
            Real tmpReal[9];
            aOver12 = a * F(0.08333333333333333333333333333333);
            aOver6 = a * F(0.16666666666666666666666666666667);

            //y1 = a(Im(m11-m22))/12, y2 = a(Im(m11-m33))/12, y3= a(Im(m22-m33))/6
            y1 = aOver12 * (m_me[0].y - m_me[4].y);
            y2 = aOver12 * (m_me[0].y - m_me[8].y);
            y3 = aOver6 * (m_me[4].y - m_me[8].y);

            //l1 = y1^2 + (|m12|/4)^2
            //l2 = y2^2 + (|m13|/4)^2
            //l3 = y3^2 + (|m23|/2)^2
            a2 = a * a;
            aSqOver16 = a2 * F(0.0625);
            aSqOver4 = a2 * F(0.25);
            l1 = y1 * y1 + aSqOver16 * __cuCabsSqf(m_me[1]);
            l2 = y2 * y2 + aSqOver16 * __cuCabsSqf(m_me[2]);
            l3 = y3 * y3 + aSqOver4 * __cuCabsSqf(m_me[5]);

            //dn1 = 1/ (1 + y1^2 +  (|m12|/4)^2)
            dn1 = __div(F(1.0), F(1.0) + l1);
            dn2 = __div(F(1.0), F(1.0) + l2);
            dn3 = __div(F(1.0), F(1.0) + l3);

            //dn1(1 - l1) = -l1 x dn1 + dn1
            Real c1r = dn1 - l1 * dn1;
            Real c2r = dn2 - l2 * dn2;
            Real c3r = dn3 - l3 * dn3;
            Real c1i = F(2.0) * y1 * dn1;
            Real c2i = F(2.0) * y2 * dn2;
            Real c3i = F(2.0) * y3 * dn3;

            halfa = F(0.5) * a;
            halfadn1 = halfa * dn1;
            halfadn2 = halfa * dn2;
            adn3 = a * dn3;
            Real a1r = halfadn1 * m_me[1].x;
            Real a2r = halfadn2 * m_me[2].x;
            Real a3r = adn3 * m_me[5].x;
            Real a1i = halfadn1 * m_me[1].y;
            Real a2i = halfadn2 * m_me[2].y;
            Real a3i = adn3 * m_me[5].y;
#else
            const Real aOver12 = a * F(0.08333333333333333333333333333333);
            const Real aOver6 = a * F(0.16666666666666666666666666666667);

            //y1 = a(Im(m11-m22))/12, y2 = a(Im(m11-m33))/12, y3= a(Im(m22-m33))/6
            const Real y1 = aOver12 * (m_me[0].y - m_me[4].y);
            const Real y2 = aOver12 * (m_me[0].y - m_me[8].y);
            const Real y3 = aOver6 * (m_me[4].y - m_me[8].y);

            //l1 = y1^2 + (|m12|/4)^2
            //l2 = y2^2 + (|m13|/4)^2
            //l3 = y3^2 + (|m23|/2)^2
            const Real a2 = a * a;
            const Real aSqOver16 = a2 * F(0.0625);
            const Real aSqOver4 = a2 * F(0.25);
            const Real l1 = y1 * y1 + aSqOver16 * __cuCabsSqf(m_me[1]);
            const Real l2 = y2 * y2 + aSqOver16 * __cuCabsSqf(m_me[2]);
            const Real l3 = y3 * y3 + aSqOver4 * __cuCabsSqf(m_me[5]);

            //dn1 = 1/ (1 + y1^2 +  (|m12|/4)^2)
            const Real dn1 = __div(F(1.0), F(1.0) + l1);
            const Real dn2 = __div(F(1.0), F(1.0) + l2);
            const Real dn3 = __div(F(1.0), F(1.0) + l3);

            //dn1(1 - l1) = -l1 x dn1 + dn1
            const Real c1r = dn1 - l1 * dn1;
            const Real c2r = dn2 - l2 * dn2;
            const Real c3r = dn3 - l3 * dn3;
            const Real c1i = F(2.0) * y1 * dn1;
            const Real c2i = F(2.0) * y2 * dn2;
            const Real c3i = F(2.0) * y3 * dn3;

            const Real halfa = F(0.5) * a;
            const Real halfadn1 = halfa * dn1;
            const Real halfadn2 = halfa * dn2;
            const Real adn3 = a * dn3;
            const Real a1r = halfadn1 * m_me[1].x;
            const Real a2r = halfadn2 * m_me[2].x;
            const Real a3r = adn3 * m_me[5].x;
            const Real a1i = halfadn1 * m_me[1].y;
            const Real a2i = halfadn2 * m_me[2].y;
            const Real a3i = adn3 * m_me[5].y;
#endif

            //6 Temp Reals
            //deviceSU3 midMatrix;
            //to reduce regcount
            CLGComplex midMatrix[8];
#if _CLG_DEBUG
            abs2a = a2r * a2r + a2i * a2i;
            abs2c = c2r * c2r - c2i * c2i;

            c2ic2r = F(2.0) * c2i * c2r;
            //(abs2c-abs2a c3r)+I(2 c2i c2r+abs2a c3i)
            midMatrix[0] = _make_cuComplex(
                abs2c - abs2a * c3r,
                c2ic2r + abs2a * c3i
            );

            //-abs2a-2 c2i c2r c3i+ c3r abs2c+I(- c3i abs2c-2 c2i c2r c3r)
            midMatrix[7] = _make_cuComplex(
                c3r * abs2c - abs2a - c3i * c2ic2r,
                -c3i * abs2c - c3r * c2ic2r
            );

            //-a2i a3i-a2r a3r+I (a2r a3i-a2i a3r)
            midMatrix[1] = _make_cuComplex(
                -a2i * a3i - a2r * a3r,
                a2r * a3i - a2i * a3r
            );

            //q1 = c2i - c2r c3i - c2i c3r
            //q2 = c2r - c2i c3i + c2r c3r
            q1 = c2i - c2r * c3i - c2i * c3r;
            q2 = c2r - c2i * c3i + c2r * c3r;

            //(a2r q2-a2i q1)+(a2i q2+a2r q1)I
            a2rq1 = a2r * q1;
            a2iq1 = a2i * q1;
            a2rq2 = a2r * q2;
            a2iq2 = a2i * q2;
            midMatrix[2] = _make_cuComplex(
                a2rq2 - a2iq1,
                a2iq2 + a2rq1
            );

            //-a2r q2-a2i q1+(a2i q2-a2r q1)I
            midMatrix[5] = _make_cuComplex(
                -a2rq2 - a2iq1,
                a2iq2 - a2rq1
            );

            //midMatrix[3] = _make_cuComplex(midMatrix[1].x, -midMatrix[1].y);
            midMatrix[3] = _make_cuComplex(c3r, c3i);

            //a3i c2i+a3r c2r+I (-a3r c2i+a3i c2r)
            a3ic2i = a3i * c2i;
            a3rc2r = a3r * c2r;
            a3rc2i = a3r * c2i;
            a3ic2r = a3i * c2r;
            midMatrix[4] = _make_cuComplex(
                a3ic2i + a3rc2r,
                a3ic2r - a3rc2i
            );

            //a3i c2i-a3r c2r+I (a3r c2i+a3i c2r)
            midMatrix[6] = _make_cuComplex(
                a3ic2i - a3rc2r,
                a3rc2i + a3ic2r
            );
#else
            const Real abs2a = a2r * a2r + a2i * a2i;
            const Real abs2c = c2r * c2r - c2i * c2i;

            const Real c2ic2r = F(2.0) * c2i * c2r;
            //(abs2c-abs2a c3r)+I(2 c2i c2r+abs2a c3i)
            midMatrix[0] = _make_cuComplex(
                abs2c - abs2a * c3r,
                c2ic2r + abs2a * c3i
            );

            //-abs2a-2 c2i c2r c3i+ c3r abs2c+I(- c3i abs2c-2 c2i c2r c3r)
            midMatrix[7] = _make_cuComplex(
                c3r * abs2c - abs2a - c3i * c2ic2r,
                -c3i * abs2c - c3r * c2ic2r
            );

            //-a2i a3i-a2r a3r+I (a2r a3i-a2i a3r)
            midMatrix[1] = _make_cuComplex(
                -a2i * a3i - a2r * a3r,
                 a2r * a3i - a2i * a3r
            );

            //q1 = c2i - c2r c3i - c2i c3r
            //q2 = c2r - c2i c3i + c2r c3r
            const Real q1 = c2i - c2r * c3i - c2i * c3r;
            const Real q2 = c2r - c2i * c3i + c2r * c3r;

            //(a2r q2-a2i q1)+(a2i q2+a2r q1)I
            const Real a2rq1 = a2r * q1;
            const Real a2iq1 = a2i * q1;
            const Real a2rq2 = a2r * q2;
            const Real a2iq2 = a2i * q2;
            midMatrix[2] = _make_cuComplex(
                a2rq2 - a2iq1,
                a2iq2 + a2rq1
            );

            //-a2r q2-a2i q1+(a2i q2-a2r q1)I
            midMatrix[5] = _make_cuComplex(
               -a2rq2 - a2iq1,
                a2iq2 - a2rq1
            );

            //midMatrix[3] = _make_cuComplex(midMatrix[1].x, -midMatrix[1].y);
            midMatrix[3] = _make_cuComplex(c3r, c3i);

            //a3i c2i+a3r c2r+I (-a3r c2i+a3i c2r)
            const Real a3ic2i = a3i * c2i;
            const Real a3rc2r = a3r * c2r;
            const Real a3rc2i = a3r * c2i;
            const Real a3ic2r = a3i * c2r;
            midMatrix[4] = _make_cuComplex(
                a3ic2i + a3rc2r,
                a3ic2r - a3rc2i
            );

            //a3i c2i-a3r c2r+I (a3r c2i+a3i c2r)
            midMatrix[6] = _make_cuComplex(
                a3ic2i - a3rc2r,
                a3rc2i + a3ic2r
            );
#endif

            //To reduce regcount
            //deviceSU3 u1Matrix = deviceSU3::makeSU3Id();
            //u1Matrix.m_me[0] = _make_cuComplex(c1r, c1i);
            //u1Matrix.m_me[4] = _make_cuComplex(c1r, -c1i);
            //u1Matrix.m_me[1] = _make_cuComplex(a1r, a1i);
            //u1Matrix.m_me[3] = _make_cuComplex(-a1r, a1i);

            //midMatrix.Mul(u1Matrix);
            //u1Matrix.Mul(midMatrix);
            CLGComplex res[9];

            //First step is midMatrix * u1Matrix
            //midMatrix is
            //0  1  2
            //1* 3  4
            //5  6  7
            //u1 Matrix is
            // c  a   0
            //-a* c*  0
            // 0  0   1

            //x = midMatrix[0], y = midMatrix[1], A = c, B = -a*
            res[0] = __TWO_COMPLEX_MULT_4(midMatrix[0].x, midMatrix[0].y, midMatrix[1].x, midMatrix[1].y, c1r, c1i, -a1r, a1i);
            //x = midMatrix[0], y = midMatrix[1], A = a, B = c*
            res[1] = __TWO_COMPLEX_MULT_4(midMatrix[0].x, midMatrix[0].y, midMatrix[1].x, midMatrix[1].y, a1r, a1i, c1r, -c1i);
            res[2] = midMatrix[2];
            //x = midMatrix[1]*, y = midMatrix[3], A = c, B = -a*
            res[3] = __TWO_COMPLEX_MULT_4(midMatrix[1].x, -midMatrix[1].y, midMatrix[3].x, midMatrix[3].y, c1r, c1i, -a1r, a1i);
            //x = midMatrix[1]*, y = midMatrix[3], A = a, B = c*
            res[4] = __TWO_COMPLEX_MULT_4(midMatrix[1].x, -midMatrix[1].y, midMatrix[3].x, midMatrix[3].y, a1r, a1i, c1r, -c1i);
            res[5] = midMatrix[4];

            //x = midMatrix[5], y = midMatrix[6], A = c, B = -a*
            res[6] = __TWO_COMPLEX_MULT_4(midMatrix[5].x, midMatrix[5].y, midMatrix[6].x, midMatrix[6].y, c1r, c1i, -a1r, a1i);
            //x = midMatrix[5], y = midMatrix[6], A = a, B = c*
            res[7] = __TWO_COMPLEX_MULT_4(midMatrix[5].x, midMatrix[5].y, midMatrix[6].x, midMatrix[6].y, a1r, a1i, c1r, -c1i);
            res[8] = midMatrix[7];

            //midMatrix is not using now, we reuse it
            //Now we calculate 
            //u1Matrix * res
            //u1 Matrix is
            // c  a   0
            //-a* c*  0
            // 0  0   1
            //res is
            //0  1  2
            //3  4  5
            //6  7  8
            //
            // c a 0 3
            midMatrix[0] = __TWO_COMPLEX_MULT_4(c1r, c1i, a1r, a1i, res[0].x, res[0].y, res[3].x, res[3].y);
            // c a 1 4
            midMatrix[1] = __TWO_COMPLEX_MULT_4(c1r, c1i, a1r, a1i, res[1].x, res[1].y, res[4].x, res[4].y);
            // c a 2 5
            midMatrix[2] = __TWO_COMPLEX_MULT_4(c1r, c1i, a1r, a1i, res[2].x, res[2].y, res[5].x, res[5].y);
            // -a* c* 0 3
            midMatrix[3] = __TWO_COMPLEX_MULT_4(-a1r, a1i, c1r, -c1i, res[0].x, res[0].y, res[3].x, res[3].y);
            // -a* c* 1 4
            midMatrix[4] = __TWO_COMPLEX_MULT_4(-a1r, a1i, c1r, -c1i, res[1].x, res[1].y, res[4].x, res[4].y);
            // -a* c* 2 5
            midMatrix[5] = __TWO_COMPLEX_MULT_4(-a1r, a1i, c1r, -c1i, res[2].x, res[2].y, res[5].x, res[5].y);
            memcpy(res, midMatrix, sizeof(CLGComplex) * 6);

            return deviceSU3(res);
        }

#if _CLG_DEBUG
#undef aOver12
#undef aOver6
#undef y1
#undef y2
#undef y3
#undef a2
#undef aSqOver16
#undef aSqOver4
#undef l1
#undef l2
#undef l3
#undef dn1
#undef dn2
#undef dn3
#undef halfa
#undef halfadn1
#undef halfadn2
#undef adn3

#undef abs2a
#undef abs2c
#undef c2ic2r

#undef q1
#undef q2
#undef a2rq1
#undef a2iq1
#undef a2rq2
#undef a2iq2

#undef a3ic2i
#undef a3rc2r
#undef a3rc2i
#undef a3ic2r
#endif

#pragma endregion

#pragma region power and log

        /**
         * Eigen values of any 3x3 matrix
         */
        __device__ __inline__ void CalculateEigenValues(CLGComplex& c1, CLGComplex& c2, CLGComplex& c3) const
        {
            const Real fUpperTraingle = __cuCabsSqf(m_me[1]) + __cuCabsSqf(m_me[2]) + __cuCabsSqf(m_me[5]);
            const Real fDownTraingle = __cuCabsSqf(m_me[3]) + __cuCabsSqf(m_me[6]) + __cuCabsSqf(m_me[7]);

            if (fUpperTraingle < _CLG_FLT_EPSILON
             || fDownTraingle < _CLG_FLT_EPSILON)
            {
                c1 = m_me[0];
                c2 = m_me[4];
                c3 = m_me[8];
                return;
            }

            //q=tr / 3
            const CLGComplex fQ = cuCmulf_cr(Tr(), F(0.33333333333333333333333));
            deviceSU3 cpy = SubCompC(fQ);
            const CLGComplex fP = __cuCsqrtf(cuCmulf_cr((cpy.MulC(cpy)).Tr(), F(0.1666666666666666666667)));
            cpy.DivComp(fP);          
            CLGComplex fDet = Determinent(cpy.m_me);

            //[ r+sqrt(r^2-4) ] ^(1/3)
            fDet = __cuCpowerf(_cuCaddf(fDet, 
                __cuCsqrtf(cuCsubf_cr(_cuCmulf(fDet, fDet), F(4.0)))), 
                F(0.333333333333333333333333));

            const CLGComplex _1over2Power13 = _make_cuComplex(
                F(1.2599210498948731647672106072782283505702514647015079800819), F(0.0));
            const CLGComplex _m1over2Power13 = 
                _make_cuComplex(F(0.629960524947436582383605303639114175285125732350753990041), 
                    F(1.0911236359717214035600726141898088813258733387403009407035));
            const CLGComplex _m1over2Power13Star = 
                _make_cuComplex(-F(0.629960524947436582383605303639114175285125732350753990041), 
                    F(1.0911236359717214035600726141898088813258733387403009407035));

            //q+p(det/r+r/det)

            c1 = _cuCaddf(fQ, _cuCmulf(fP,
                _cuCaddf(_cuCdivf(fDet, _m1over2Power13Star), _cuCdivf(_m1over2Power13Star, fDet))
                ));
            c2 = _cuCaddf(fQ, _cuCmulf(fP,
                _cuCaddf(_cuCdivf(fDet, _1over2Power13), _cuCdivf(_1over2Power13, fDet))
            ));
            c3 = _cuCsubf(fQ, _cuCmulf(fP,
                _cuCaddf(_cuCdivf(fDet, _m1over2Power13), _cuCdivf(_m1over2Power13, fDet))
            ));
        }

        __device__ __inline__ deviceSU3 EigenVectors(
            const CLGComplex& c1, 
            const CLGComplex& c2, 
            const CLGComplex& c3) const
        {
            deviceSU3 tmpM1;
            deviceSU3 retRes;
            
            CLGComplex right_m_lambda1 = _cuCsubf(m_me[0], c3);
            CLGComplex right_m_lambda2 = _cuCsubf(m_me[4], c3);
            CLGComplex right_m_lambda3 = _cuCsubf(m_me[8], c3);
            CLGComplex left_m_lambda = _cuCsubf(m_me[0], c2);
            __LINE_MUL_POLY(tmpM1.m_me[0], left_m_lambda, m_me[1], m_me[2], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[1], left_m_lambda, m_me[1], m_me[2], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[2], left_m_lambda, m_me[1], m_me[2], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[4], c2);
            __LINE_MUL_POLY(tmpM1.m_me[3], m_me[3], left_m_lambda, m_me[5], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[4], m_me[3], left_m_lambda, m_me[5], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[5], m_me[3], left_m_lambda, m_me[5], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[8], c2);
            __LINE_MUL_POLY(tmpM1.m_me[6], m_me[6], m_me[7], left_m_lambda, right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[7], m_me[6], m_me[7], left_m_lambda, m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[8], m_me[6], m_me[7], left_m_lambda, m_me[2], m_me[5], right_m_lambda3);

            //v1 = 0, 3, 6
            retRes.m_me[0] = _cuCaddf(_cuCaddf(tmpM1.m_me[0], tmpM1.m_me[1]), tmpM1.m_me[2]);
            retRes.m_me[3] = _cuCaddf(_cuCaddf(tmpM1.m_me[3], tmpM1.m_me[4]), tmpM1.m_me[5]);
            retRes.m_me[6] = _cuCaddf(_cuCaddf(tmpM1.m_me[6], tmpM1.m_me[7]), tmpM1.m_me[8]);
            if (__cuCabsSqf(retRes.m_me[0]) > _CLG_FLT_EPSILON)
            {
                retRes.m_me[3] = _cuCdivf(retRes.m_me[3], retRes.m_me[0]);
                retRes.m_me[6] = _cuCdivf(retRes.m_me[6], retRes.m_me[0]);
                retRes.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
            }
            Real fDenorm = __cuCabsSqf(retRes.m_me[0]) + __cuCabsSqf(retRes.m_me[3]) + __cuCabsSqf(retRes.m_me[6]);
            if (fDenorm > _CLG_FLT_EPSILON)
            {
                fDenorm = __rcp(_sqrt(fDenorm));
                retRes.m_me[0] = cuCmulf_cr(retRes.m_me[0], fDenorm);
                retRes.m_me[3] = cuCmulf_cr(retRes.m_me[3], fDenorm);
                retRes.m_me[6] = cuCmulf_cr(retRes.m_me[6], fDenorm);
            }
            else
            {
                retRes.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
                retRes.m_me[6] = _make_cuComplex(F(0.0), F(0.0));
            }


            right_m_lambda1 = _cuCsubf(m_me[0], c1);
            right_m_lambda2 = _cuCsubf(m_me[4], c1);
            right_m_lambda3 = _cuCsubf(m_me[8], c1);
            left_m_lambda = _cuCsubf(m_me[0], c3);
            __LINE_MUL_POLY(tmpM1.m_me[0], left_m_lambda, m_me[1], m_me[2], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[1], left_m_lambda, m_me[1], m_me[2], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[2], left_m_lambda, m_me[1], m_me[2], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[4], c3);
            __LINE_MUL_POLY(tmpM1.m_me[3], m_me[3], left_m_lambda, m_me[5], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[4], m_me[3], left_m_lambda, m_me[5], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[5], m_me[3], left_m_lambda, m_me[5], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[8], c3);
            __LINE_MUL_POLY(tmpM1.m_me[6], m_me[6], m_me[7], left_m_lambda, right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[7], m_me[6], m_me[7], left_m_lambda, m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[8], m_me[6], m_me[7], left_m_lambda, m_me[2], m_me[5], right_m_lambda3);
            retRes.m_me[1] = _cuCaddf(_cuCaddf(tmpM1.m_me[0], tmpM1.m_me[1]), tmpM1.m_me[2]);
            retRes.m_me[4] = _cuCaddf(_cuCaddf(tmpM1.m_me[3], tmpM1.m_me[4]), tmpM1.m_me[5]);
            retRes.m_me[7] = _cuCaddf(_cuCaddf(tmpM1.m_me[6], tmpM1.m_me[7]), tmpM1.m_me[8]);
            CLGComplex v1v2 = _cuCaddf(_cuCaddf(
                _cuCmulf(_cuConjf(retRes.m_me[0]), retRes.m_me[1]), 
                _cuCmulf(_cuConjf(retRes.m_me[3]), retRes.m_me[4])),
                _cuCmulf(_cuConjf(retRes.m_me[6]), retRes.m_me[7]));
            retRes.m_me[1] = _cuCsubf(retRes.m_me[1], _cuCmulf(v1v2, retRes.m_me[0]));
            retRes.m_me[4] = _cuCsubf(retRes.m_me[4], _cuCmulf(v1v2, retRes.m_me[3]));
            retRes.m_me[7] = _cuCsubf(retRes.m_me[7], _cuCmulf(v1v2, retRes.m_me[6]));
            if (__cuCabsSqf(retRes.m_me[4]) > _CLG_FLT_EPSILON)
            {
                retRes.m_me[1] = _cuCdivf(retRes.m_me[1], retRes.m_me[4]);
                retRes.m_me[7] = _cuCdivf(retRes.m_me[7], retRes.m_me[4]);
                retRes.m_me[4] = _make_cuComplex(F(1.0), F(0.0));
            }
            fDenorm = __cuCabsSqf(retRes.m_me[1]) + __cuCabsSqf(retRes.m_me[4]) + __cuCabsSqf(retRes.m_me[7]);
            if (fDenorm > _CLG_FLT_EPSILON)
            {
                fDenorm = __rcp(_sqrt(fDenorm));
                retRes.m_me[1] = cuCmulf_cr(retRes.m_me[1], fDenorm);
                retRes.m_me[4] = cuCmulf_cr(retRes.m_me[4], fDenorm);
                retRes.m_me[7] = cuCmulf_cr(retRes.m_me[7], fDenorm);
            }
            else
            {
                retRes.m_me[4] = _make_cuComplex(F(0.0), F(0.0));
                retRes.m_me[7] = _make_cuComplex(F(0.0), F(0.0));
            }


            right_m_lambda1 = _cuCsubf(m_me[0], c2);
            right_m_lambda2 = _cuCsubf(m_me[4], c2);
            right_m_lambda3 = _cuCsubf(m_me[8], c2);
            left_m_lambda = _cuCsubf(m_me[0], c1);
            __LINE_MUL_POLY(tmpM1.m_me[0], left_m_lambda, m_me[1], m_me[2], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[1], left_m_lambda, m_me[1], m_me[2], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[2], left_m_lambda, m_me[1], m_me[2], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[4], c1);
            __LINE_MUL_POLY(tmpM1.m_me[3], m_me[3], left_m_lambda, m_me[5], right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[4], m_me[3], left_m_lambda, m_me[5], m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[5], m_me[3], left_m_lambda, m_me[5], m_me[2], m_me[5], right_m_lambda3);
            left_m_lambda = _cuCsubf(m_me[8], c1);
            __LINE_MUL_POLY(tmpM1.m_me[6], m_me[6], m_me[7], left_m_lambda, right_m_lambda1, m_me[3], m_me[6]);
            __LINE_MUL_POLY(tmpM1.m_me[7], m_me[6], m_me[7], left_m_lambda, m_me[1], right_m_lambda2, m_me[7]);
            __LINE_MUL_POLY(tmpM1.m_me[8], m_me[6], m_me[7], left_m_lambda, m_me[2], m_me[5], right_m_lambda3);

            retRes.m_me[2] = _cuCaddf(_cuCaddf(tmpM1.m_me[0], tmpM1.m_me[1]), tmpM1.m_me[2]);
            retRes.m_me[5] = _cuCaddf(_cuCaddf(tmpM1.m_me[3], tmpM1.m_me[4]), tmpM1.m_me[5]);
            retRes.m_me[8] = _cuCaddf(_cuCaddf(tmpM1.m_me[6], tmpM1.m_me[7]), tmpM1.m_me[8]);
            v1v2 = _cuCaddf(_cuCaddf(
                _cuCmulf(_cuConjf(retRes.m_me[0]), retRes.m_me[2]),
                _cuCmulf(_cuConjf(retRes.m_me[3]), retRes.m_me[5])),
                _cuCmulf(_cuConjf(retRes.m_me[6]), retRes.m_me[8]));
            const CLGComplex v2v3 = _cuCaddf(_cuCaddf(
                _cuCmulf(_cuConjf(retRes.m_me[1]), retRes.m_me[2]),
                _cuCmulf(_cuConjf(retRes.m_me[4]), retRes.m_me[5])),
                _cuCmulf(_cuConjf(retRes.m_me[7]), retRes.m_me[8]));
            retRes.m_me[2] = _cuCsubf(retRes.m_me[2], _cuCmulf(v1v2, retRes.m_me[0]));
            retRes.m_me[5] = _cuCsubf(retRes.m_me[5], _cuCmulf(v1v2, retRes.m_me[3]));
            retRes.m_me[8] = _cuCsubf(retRes.m_me[8], _cuCmulf(v1v2, retRes.m_me[6]));
            retRes.m_me[2] = _cuCsubf(retRes.m_me[2], _cuCmulf(v2v3, retRes.m_me[1]));
            retRes.m_me[5] = _cuCsubf(retRes.m_me[5], _cuCmulf(v2v3, retRes.m_me[4]));
            retRes.m_me[8] = _cuCsubf(retRes.m_me[8], _cuCmulf(v2v3, retRes.m_me[7]));

            if (__cuCabsSqf(retRes.m_me[8]) > _CLG_FLT_EPSILON)
            {
                retRes.m_me[2] = _cuCdivf(retRes.m_me[2], retRes.m_me[8]);
                retRes.m_me[5] = _cuCdivf(retRes.m_me[5], retRes.m_me[8]);
                retRes.m_me[8] = _make_cuComplex(F(1.0), F(0.0));
            }
            fDenorm = __cuCabsSqf(retRes.m_me[2]) + __cuCabsSqf(retRes.m_me[5]) + __cuCabsSqf(retRes.m_me[8]);
            if (fDenorm > _CLG_FLT_EPSILON)
            {
                fDenorm = __rcp(_sqrt(fDenorm));
                retRes.m_me[2] = cuCmulf_cr(retRes.m_me[2], fDenorm);
                retRes.m_me[5] = cuCmulf_cr(retRes.m_me[5], fDenorm);
                retRes.m_me[8] = cuCmulf_cr(retRes.m_me[8], fDenorm);
            }
            else
            {
                retRes.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
                retRes.m_me[5] = _make_cuComplex(F(0.0), F(0.0));
            }

            return retRes;
        }

        __device__ __inline__ deviceSU3 Power(Real fPower) const
        {
            CLGComplex c1;
            CLGComplex c2;
            CLGComplex c3;
            CalculateEigenValues(c1, c2, c3);
            const deviceSU3 diagonal = EigenVectors(c1, c2, c3);

            deviceSU3 ret;
            ret.m_me[0] = __cuCpowerf(c1, fPower);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = _zeroc;
            ret.m_me[4] = __cuCpowerf(c2, fPower);
            ret.m_me[5] = _zeroc;
            ret.m_me[6] = _zeroc;
            ret.m_me[7] = _zeroc;
            ret.m_me[8] = __cuCpowerf(c3, fPower);
            ret.MulDagger(diagonal);
            ret = diagonal.MulC(ret);
            return ret;
        }

        __device__ __inline__ deviceSU3 Log() const
        {
            CLGComplex c1;
            CLGComplex c2;
            CLGComplex c3;
            CalculateEigenValues(c1, c2, c3);
            const deviceSU3 diagonal = EigenVectors(c1, c2, c3);
            deviceSU3 ret;
            ret.m_me[0] = __cuClogf(c1);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = _zeroc;
            ret.m_me[4] = __cuClogf(c2);
            ret.m_me[5] = _zeroc;
            ret.m_me[6] = _zeroc;
            ret.m_me[7] = _zeroc;
            ret.m_me[8] = __cuClogf(c3);
            ret.MulDagger(diagonal);
            ret = diagonal.MulC(ret);
            return ret;
        }

        __device__ __inline__ deviceSU3 StrictExp() const
        {
            const Real fAbsAll = 
              __cuCabsSqf(m_me[1]) + __cuCabsSqf(m_me[2]) + __cuCabsSqf(m_me[3])
            + __cuCabsSqf(m_me[4]) + __cuCabsSqf(m_me[5]) + __cuCabsSqf(m_me[6])
            + __cuCabsSqf(m_me[7]) + __cuCabsSqf(m_me[8]) + __cuCabsSqf(m_me[9]);

            //if A is small, there will be problems
            if (fAbsAll < F(0.0000001))
            {
                return QuickExp(F(1.0));
            }

            CLGComplex c1;
            CLGComplex c2;
            CLGComplex c3;
            CalculateEigenValues(c1, c2, c3);
            const deviceSU3 diagonal = EigenVectors(c1, c2, c3);

            deviceSU3 ret;
            ret.m_me[0] = __cuCexpf(c1);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = _zeroc;
            ret.m_me[4] = __cuCexpf(c2);
            ret.m_me[5] = _zeroc;
            ret.m_me[6] = _zeroc;
            ret.m_me[7] = _zeroc;
            ret.m_me[8] = __cuCexpf(c3);
            ret.MulDagger(diagonal);
            ret = diagonal.MulC(ret);
            return ret;
        }

#pragma endregion

        //union
        //{
        //    deviceSU3Vector m_v[3];
        //    CLGComplex m_me[9];
        //};
        CLGComplex m_me[16]; //Only the first 9 elememt is using, 16 is for padding and align
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
