//=============================================================================
// FILENAME : SU2.h
// 
// DESCRIPTION:
// This is helper functions for calculate D=2 complex matrix
// NOTE: It is not nessary a SU3
//
// The SU2 Matrix is
// 0 1
// 2 3
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _SU2_H_
#define _SU2_H_

#define __LINE_MUL2(a, b, c, d) _cuCaddf(_cuCmulf(m_me[a], right.m_me[c]), _cuCmulf(m_me[b], right.m_me[d]))

#define __LINE_MULARRAY2(a, b, c, d) _cuCaddf(_cuCmulf(left[a], right[c]), _cuCmulf(left[b], right[d]))

#define __LINE_MULND2(a, b, c, d) _cuCaddf(_cuCmulf(m_me[a], _cuConjf(right.m_me[c])), _cuCmulf(m_me[b], _cuConjf(right.m_me[d])))

#define __LINE_MULDN2(a, b, c, d) _cuCaddf(_cuCmulf(_cuConjf(m_me[a]), right.m_me[c]), _cuCmulf(_cuConjf(m_me[b]), right.m_me[d]))

//self^dagger.self
#define __LINE_MULDN_ME2(a, b, c, d) _cuCaddf(_cuCmulf(_cuConjf(m_me[a]), m_me[c]), _cuCmulf(_cuConjf(m_me[b]), m_me[d]))


__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    struct deviceSU2
    {
        __device__ deviceSU2()
        {

        }

        __device__ deviceSU2(const deviceSU2& other)
        {
            memcpy(m_me, other.m_me, sizeof(CLGComplex) * 4);
        }

        __device__ deviceSU2(const CLGComplex* __restrict__ other)
        {
            memcpy(m_me, other, sizeof(CLGComplex) * 4);
        }

        __device__ void DebugPrint(const char* header = NULL) const
        {
            printf("%s%s{{%1.7f%s%1.7f I, %1.7f%s%1.7f I},\n {%1.7f%s%1.7f I, %1.7f%s%1.7f I}};\n",
                NULL == header ? "" : header,
                NULL == header ? "" : "=",

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
                m_me[3].y
            );
        }

#pragma region create

        /**
        * ret = 0
        */
        __device__ __inline__ static deviceSU2 makeSU2Zero()
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
            return ret;
        }

        /**
        * ret = I
        */
        __device__ __inline__ static deviceSU2 makeSU2Id()
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
            ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
            return ret;
        }

        /**
        * can be called only after CLatticeData is created
        * ret = random
        */
        __device__ __inline__ static deviceSU2 makeSU2Random(UINT fatIndex)
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[1] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[2] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[3] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.Norm();
            return ret;
        }

        __device__ __inline__ static deviceSU2 makeSU2RandomAny(UINT fatIndex)
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[1] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[2] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            ret.m_me[3] = _make_cuComplex(_deviceRandomF(fatIndex) - F(0.5), _deviceRandomF(fatIndex) - F(0.5));
            return ret;
        }

        /**
        * Here, we keep same with Bridge++, that H(P)/D.O.F. = 0.5
        * can be called only after CLatticeData is created
        * ret = i r_a T_a, r_a is random real number, T_a are generators
        * r_a T_a = r1 T1 + r2 T2 + ...
        *
        *     r3         r1-ir2
        * =   r1+ir2     -r3
        *
        */
        __device__ __inline__ static deviceSU2 makeSU2RandomGenerator(UINT fatIndex)
        {
            const Real r1 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r2 = _deviceRandomGaussFSqrt2(fatIndex);
            const Real r3 = _deviceRandomGaussFSqrt2(fatIndex);

            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(F(0.0), r3);
            ret.m_me[1] = _make_cuComplex(r2, r1);
            ret.m_me[2] = _make_cuComplex(-r2, r1);
            ret.m_me[3] = _make_cuComplex(F(0.0), -r3);
            return ret;
        }

        __device__ __inline__ static deviceSU2 makeSU2SumGenerator(Real fDivide)
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex( fDivide, F(0.0));
            ret.m_me[1] = _make_cuComplex(fDivide, -fDivide);
            ret.m_me[2] = _make_cuComplex(fDivide, fDivide);
            ret.m_me[3] = _make_cuComplex(-fDivide, F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU2 makeSU3Generator(UINT uiGenerator)
        {
            deviceSU2 ret = makeSU2Zero();
            const Real half = F(1.0);
            switch (uiGenerator)
            {
            case 0:
            {
                /**
                *     0     1 
                * =   1     0  
                */
                ret.m_me[1] = _make_cuComplex(half, F(0.0));
                ret.m_me[2] = _make_cuComplex(half, F(0.0));
            }
            break;
            case 1:
            {
                /**
                *     0     -I  
                * =   I     0   
                */
                ret.m_me[1] = _make_cuComplex(F(0.0), -half);
                ret.m_me[2] = _make_cuComplex(F(0.0), half);
            }
            break;
            case 2:
            {
                /**
                *     1     0 
                * =   0    -1
                */
                ret.m_me[0] = _make_cuComplex(half, F(0.0));
                ret.m_me[3] = _make_cuComplex(-half, F(0.0));
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
            default:
            {
                ret.m_me[0] = _make_cuComplex(half, F(0.0));
                ret.m_me[1] = _make_cuComplex(half, -half);
                ret.m_me[2] = _make_cuComplex(half, half);
                ret.m_me[3] = _make_cuComplex(-half, F(0.0));
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
        __device__ __inline__ static deviceSU2 makeSU2TA(const CLGComplex& m12, Real m11)
        {
            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(F(0.0), m11);;
            ret.m_me[1] = m12;
            ret.m_me[2] = _make_cuComplex(-m12.x, m12.y);
            ret.m_me[3] = _make_cuComplex(F(0.0), -m11);
            return ret;
        }

        /**
         * Another version of makeSU3Contract without spin sum
         * if v = {{},{},{}}, column vector
         * This is right.left^+
         */
        __device__ __inline__ static deviceSU2 makeSU2ContractV(const deviceSU2Vector& left, const deviceSU2Vector& right)
        {
            deviceSU2 ret;
            ret.m_me[0] = _cuCmulf(_cuConjf(left.m_ve[0]), right.m_ve[0]);
            ret.m_me[1] = _cuCmulf(_cuConjf(left.m_ve[1]), right.m_ve[0]);
            ret.m_me[2] = _cuCmulf(_cuConjf(left.m_ve[0]), right.m_ve[1]);
            ret.m_me[3] = _cuCmulf(_cuConjf(left.m_ve[1]), right.m_ve[1]);           
            return ret;
        }

#pragma endregion create

#pragma region operators

        __device__ __inline__ void Add(const deviceSU2& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right.m_me[0]);
            m_me[1] = _cuCaddf(m_me[1], right.m_me[1]);
            m_me[2] = _cuCaddf(m_me[2], right.m_me[2]);
            m_me[3] = _cuCaddf(m_me[3], right.m_me[3]);
        }

        __device__ __inline__ void AddDagger(const deviceSU2& right)
        {
            m_me[0] = _cuCaddf(m_me[0], _cuConjf(right.m_me[0]));
            m_me[1] = _cuCaddf(m_me[1], _cuConjf(right.m_me[2]));
            m_me[2] = _cuCaddf(m_me[2], _cuConjf(right.m_me[1]));
            m_me[3] = _cuCaddf(m_me[3], _cuConjf(right.m_me[3]));
        }

        __device__ __inline__ void Sub(const deviceSU2& right)
        {
            m_me[0] = _cuCsubf(m_me[0], right.m_me[0]);
            m_me[1] = _cuCsubf(m_me[1], right.m_me[1]);
            m_me[2] = _cuCsubf(m_me[2], right.m_me[2]);
            m_me[3] = _cuCsubf(m_me[3], right.m_me[3]);
        }

        __device__ __inline__ void SubDagger(const deviceSU2& right)
        {
            m_me[0] = _cuCsubf(m_me[0], _cuConjf(right.m_me[0]));
            m_me[1] = _cuCsubf(m_me[1], _cuConjf(right.m_me[2]));
            m_me[2] = _cuCsubf(m_me[2], _cuConjf(right.m_me[1]));
            m_me[3] = _cuCsubf(m_me[3], _cuConjf(right.m_me[3]));
        }

        __device__ __inline__ void AddReal(Real right)
        {
            m_me[0] = cuCaddf_cr(m_me[0], right);
            m_me[3] = cuCaddf_cr(m_me[3], right);
        }

        __device__ __inline__ void AddComp(const CLGComplex& right)
        {
            m_me[0] = _cuCaddf(m_me[0], right);
            m_me[3] = _cuCaddf(m_me[3], right);
        }

        __device__ __inline__ void SubReal(Real right)
        {
            m_me[0] = cuCsubf_cr(m_me[0], right);
            m_me[3] = cuCsubf_cr(m_me[3], right);
        }

        __device__ __inline__ void SubComp(const CLGComplex& right)
        {
            m_me[0] = _cuCsubf(m_me[0], right);
            m_me[3] = _cuCsubf(m_me[3], right);
        }

        /**
        * left = left * right
        */
        __device__ __inline__ void Mul(const deviceSU2& right)
        {
            CLGComplex res[4];
            res[0] = __LINE_MUL2(0, 1, 0, 2);
            res[1] = __LINE_MUL2(0, 1, 1, 3);

            res[2] = __LINE_MUL2(2, 3, 0, 2);
            res[3] = __LINE_MUL2(2, 3, 1, 3);
            memcpy(m_me, res, sizeof(CLGComplex) * 4);
        }

        /**
        * low reg count version of multiply
        */
        __device__ __inline__ static void ArrayMul(CLGComplex* left, const CLGComplex* __restrict__ right)
        {
            CLGComplex tmp[2];
            tmp[0] = __LINE_MULARRAY2(0, 1, 0, 2);
            tmp[1] = __LINE_MULARRAY2(0, 1, 1, 3);
            memcpy(left, tmp, sizeof(CLGComplex) * 2);

            tmp[0] = __LINE_MULARRAY2(2, 3, 0, 2);
            tmp[1] = __LINE_MULARRAY2(2, 3, 1, 3);
            memcpy(left + 2, tmp, sizeof(CLGComplex) * 2);
        }

        __device__ __inline__ deviceSU2 MulC(const deviceSU2& right) const
        {
            deviceSU2 ret;
            ret.m_me[0] = __LINE_MUL2(0, 1, 0, 2);
            ret.m_me[1] = __LINE_MUL2(0, 1, 1, 3);

            ret.m_me[2] = __LINE_MUL2(2, 3, 0, 2);
            ret.m_me[3] = __LINE_MUL2(2, 3, 1, 3);

            return ret;
        }

        __device__ __inline__ void MulReal(Real right)
        {
            m_me[0] = cuCmulf_cr(m_me[0], right);
            m_me[1] = cuCmulf_cr(m_me[1], right);
            m_me[2] = cuCmulf_cr(m_me[2], right);
            m_me[3] = cuCmulf_cr(m_me[3], right);
        }

        __device__ __inline__ void MulComp(const CLGComplex& right)
        {
            m_me[0] = _cuCmulf(m_me[0], right);
            m_me[1] = _cuCmulf(m_me[1], right);
            m_me[2] = _cuCmulf(m_me[2], right);
            m_me[3] = _cuCmulf(m_me[3], right);
        }

        __device__ __inline__ void DivComp(const CLGComplex& right)
        {
            m_me[0] = _cuCdivf(m_me[0], right);
            m_me[1] = _cuCdivf(m_me[1], right);
            m_me[2] = _cuCdivf(m_me[2], right);
            m_me[3] = _cuCdivf(m_me[3], right);
        }

        __device__ __inline__ void MulDagger(const deviceSU2& right)
        {
            CLGComplex res[4];
            res[0] = __LINE_MULND2(0, 1, 0, 1);
            res[1] = __LINE_MULND2(0, 1, 2, 3);

            res[2] = __LINE_MULND2(2, 3, 0, 1);
            res[3] = __LINE_MULND2(2, 3, 2, 3);

            memcpy(m_me, res, sizeof(CLGComplex) * 4);
        }

        __device__ __inline__ deviceSU2 MulDaggerC(const deviceSU2& right) const
        {
            deviceSU2 ret;
            ret.m_me[0] = __LINE_MULND2(0, 1, 0, 1);
            ret.m_me[1] = __LINE_MULND2(0, 1, 2, 3);

            ret.m_me[2] = __LINE_MULND2(2, 3, 0, 1);
            ret.m_me[3] = __LINE_MULND2(2, 3, 2, 3);

            return ret;
        }

        __device__ __inline__ void DaggerMul(const deviceSU2& right)
        {
            CLGComplex res[4];
            res[0] = __LINE_MULDN2(0, 2, 0, 2);
            res[1] = __LINE_MULDN2(0, 2, 1, 3);

            res[2] = __LINE_MULDN2(1, 3, 0, 2);
            res[3] = __LINE_MULDN2(1, 3, 1, 3);

            memcpy(m_me, res, sizeof(CLGComplex) * 4);
        }

        __device__ __inline__ deviceSU2 DaggerMulC(const deviceSU2& right) const
        {
            deviceSU2 ret;
            ret.m_me[0] = __LINE_MULDN2(0, 2, 0, 2);
            ret.m_me[1] = __LINE_MULDN2(0, 2, 1, 3);

            ret.m_me[2] = __LINE_MULDN2(1, 3, 0, 2);
            ret.m_me[3] = __LINE_MULDN2(1, 3, 1, 3);

            return ret;
        }

        __device__ __inline__ deviceSU2Vector MulVector(const deviceSU2Vector& v) const
        {
            deviceSU2Vector ret;
            ret.m_ve[0] = _cuCaddf(_cuCmulf(m_me[0], v.m_ve[0]), _cuCmulf(m_me[1], v.m_ve[1]));
            ret.m_ve[1] = _cuCaddf(_cuCmulf(m_me[2], v.m_ve[0]), _cuCmulf(m_me[3], v.m_ve[1]));
            return ret;
        }

        __device__ __inline__ deviceSU2 AddC(const deviceSU2& right) const { deviceSU2 ret(*this); ret.Add(right); return ret; }
        __device__ __inline__ deviceSU2 AddCompC(const CLGComplex& right) const { deviceSU2 ret(*this); ret.AddComp(right); return ret; }
        __device__ __inline__ deviceSU2 AddRealC(const Real& right) const { deviceSU2 ret(*this); ret.AddReal(right); return ret; }

        __device__ __inline__ deviceSU2 SubC(const deviceSU2& right) const { deviceSU2 ret(*this); ret.Sub(right); return ret; }
        __device__ __inline__ deviceSU2 SubCompC(const CLGComplex& right) const { deviceSU2 ret(*this); ret.SubComp(right); return ret; }
        __device__ __inline__ deviceSU2 SubRealC(const Real& right) const { deviceSU2 ret(*this); ret.SubReal(right); return ret; }

        __device__ __inline__ deviceSU2 MulCompC(const CLGComplex& right) const { deviceSU2 ret(*this); ret.MulComp(right); return ret; }
        __device__ __inline__ deviceSU2 MulRealC(const Real& right) const { deviceSU2 ret(*this); ret.MulReal(right); return ret; }
        __device__ __inline__ deviceSU2 DivCompC(const CLGComplex& right) const { deviceSU2 ret(*this); ret.DivComp(right); return ret; }

        //the reduce regcount version
        __device__ __inline__ void MulRealCArray(CLGComplex* res, const Real& right) const 
        { 
            res[0] = _make_cuComplex(m_me[0].x * right, m_me[0].y * right);
            res[1] = _make_cuComplex(m_me[1].x * right, m_me[1].y * right);
            res[2] = _make_cuComplex(m_me[2].x * right, m_me[2].y * right);
            res[3] = _make_cuComplex(m_me[3].x * right, m_me[3].y * right);
        }

#pragma endregion operators

#pragma region useful functions

        /**
        * ret = ||U||, where U is a d=3 matrix
        * we do NOT need ||SU3|| because ||SU3||=1
        */
        __device__ static __inline__ CLGComplex Determinent(const CLGComplex* u)
        {
            return _cuCsubf(_cuCmulf(u[0], u[3]), _cuCmulf(u[1], u[2]));
        }

        __device__ __inline__ CLGComplex Tr() const
        {
            return _cuCaddf(m_me[0], m_me[3]);
        }

        /**
        * Re[Tr[U]]
        */
        __device__ __inline__ Real ReTr() const
        {
            return m_me[0].x + m_me[3].x;
        }

        __device__ __inline__ Real ImTr() const
        {
            return m_me[0].y + m_me[3].y;
        }

        __device__ __inline__ void Re()
        {
            m_me[0] = _make_cuComplex(m_me[0].x, F(0.0));
            m_me[1] = _make_cuComplex(m_me[1].x, F(0.0));
            m_me[2] = _make_cuComplex(m_me[2].x, F(0.0));
            m_me[3] = _make_cuComplex(m_me[3].x, F(0.0));
        }

        __device__ __inline__ deviceSU2 ReC() const { deviceSU2 ret(*this); ret.Re(); return ret; }

        __device__ __inline__ void Im()
        {
            m_me[0] = _make_cuComplex(m_me[0].y, F(0.0));
            m_me[1] = _make_cuComplex(m_me[1].y, F(0.0));
            m_me[2] = _make_cuComplex(m_me[2].y, F(0.0));
            m_me[3] = _make_cuComplex(m_me[3].y, F(0.0));
        }

        __device__ __inline__ deviceSU2 ImC() const { deviceSU2 ret(*this); ret.Im(); return ret; }

        /**
        * res = Conjugate[Transpose[U]]
        */
        __device__ __inline__ void Dagger()
        {
            CLGComplex res[4];
            res[0] = _cuConjf(m_me[0]);
            res[1] = _cuConjf(m_me[2]);
            res[2] = _cuConjf(m_me[1]);
            res[3] = _cuConjf(m_me[3]);
            memcpy(m_me, res, sizeof(CLGComplex) * 4);
        }

        __device__ __inline__ deviceSU2 DaggerC() const
        {
            deviceSU2 tmp;
            tmp.m_me[0] = _cuConjf(m_me[0]);
            tmp.m_me[1] = _cuConjf(m_me[2]);
            tmp.m_me[2] = _cuConjf(m_me[1]);
            tmp.m_me[3] = _cuConjf(m_me[3]);
            return tmp;
        }

        __device__ __inline deviceSU2 Transpose() const
        {

            deviceSU2 tmp;
            tmp.m_me[0] = m_me[0];
            tmp.m_me[1] = m_me[2];
            tmp.m_me[2] = m_me[1];
            tmp.m_me[3] = m_me[3];
            return tmp;
        }

        /**
        * A=(U-U^dagger)/2
        * res = A - tr(A)/2
        */
        __device__ __inline__ void Ta()
        {
            m_me[1].x = F(0.5) * (m_me[1].x - m_me[2].x);
            m_me[1].y = F(0.5) * (m_me[1].y + m_me[2].y);
            m_me[2].x = -m_me[1].x;
            m_me[2].y = m_me[1].y;
            Real sub = F(0.5) * (m_me[0].y + m_me[3].y);
            m_me[0].x = 0;
            m_me[0].y = m_me[0].y - sub;
            m_me[3].x = 0;
            m_me[3].y = m_me[3].y - sub;
        }

        /**
        * Return Tr[Im[a].Im[b]] 
        * It is used in topological charge
        */
        //__device__ __inline__ static Real TrIm(const deviceSU2& a, const deviceSU2&b)
        //{
        //    //0,1,2 * 0,3,6
        //    Real ret = a.m_me[0].y * b.m_me[0].y + a.m_me[1].y * b.m_me[3].y + a.m_me[2].y * b.m_me[6].y;
        //    //3,4,5 * 1,4,7
        //    ret += a.m_me[3].y * b.m_me[1].y + a.m_me[4].y * b.m_me[4].y + a.m_me[5].y * b.m_me[7].y;
        //    //6,7,8 * 2,5,8
        //    ret += a.m_me[6].y * b.m_me[2].y + a.m_me[7].y * b.m_me[5].y + a.m_me[8].y * b.m_me[8].y;

        //    return ret;
        //}

        /**
        * res = U - U^dagger
        * It is like a matrix Im(M)
        * res = 2i Im(M), so is called iIm2
        */
        //__device__ __inline__ void iIm2()
        //{
        //    //0 1 2
        //    //3 4 5
        //    //6 7 8

        //    //new [1] = [1] - conj([3])
        //    //new [3] = [3] - conj([1]) = -conj(new [1])
        //    const CLGComplex new1 = _cuCsubf(m_me[1], _cuConjf(m_me[3]));
        //    const CLGComplex new2 = _cuCsubf(m_me[2], _cuConjf(m_me[6]));
        //    const CLGComplex new5 = _cuCsubf(m_me[5], _cuConjf(m_me[7]));
        //    m_me[1] = _make_cuComplex(_cuCrealf(new1), _cuCimagf(new1));
        //    m_me[3] = _make_cuComplex(-_cuCrealf(m_me[1]), _cuCimagf(m_me[1]));
        //    m_me[2] = _make_cuComplex(_cuCrealf(new2), _cuCimagf(new2));
        //    m_me[6] = _make_cuComplex(-_cuCrealf(m_me[2]), _cuCimagf(m_me[2]));
        //    m_me[5] = _make_cuComplex(_cuCrealf(new5), _cuCimagf(new5));
        //    m_me[7] = _make_cuComplex(-_cuCrealf(m_me[5]), _cuCimagf(m_me[5]));

        //    m_me[0] = _make_cuComplex(F(0.0), F(2.0) * m_me[0].y);
        //    m_me[4] = _make_cuComplex(F(0.0), F(2.0) * m_me[4].y);
        //    m_me[8] = _make_cuComplex(F(0.0), F(2.0) * m_me[8].y);
        //}

        /**
        * return -i(U-U^dagger) = ((-iU)+(-iU)dagger)
        */
        //__device__ __inline__ deviceSU2 Im2C() const
        //{
        //    deviceSU2 ret;
        //    //0 1 2
        //    //3 4 5
        //    //6 7 8
        //    // -i(x-y*)=i(y*-x)=i(yr-iyi-xr-ixi)=(yi+xi)+i(yr-xr)
        //    // new 1 = -i(1-3*)
        //    // new 2 = -i(2-6*)
        //    // new 5 = -i(5-7*)
        //    // (new 1)*= i(1*-3) = -i(3-1*) = new 3
        //    // new 0 = -i(0-0*) = 2Im(0)

        //    ret.m_me[1] = _make_cuComplex(m_me[1].y + m_me[3].y, m_me[3].x - m_me[1].x);
        //    ret.m_me[2] = _make_cuComplex(m_me[2].y + m_me[6].y, m_me[6].x - m_me[2].x);
        //    ret.m_me[5] = _make_cuComplex(m_me[5].y + m_me[7].y, m_me[7].x - m_me[5].x);
        //    ret.m_me[3] = _cuConjf(ret.m_me[1]);
        //    ret.m_me[6] = _cuConjf(ret.m_me[2]);
        //    ret.m_me[7] = _cuConjf(ret.m_me[5]);

        //    ret.m_me[0] = _make_cuComplex(F(2.0) * m_me[0].y, F(0.0));
        //    ret.m_me[4] = _make_cuComplex(F(2.0) * m_me[4].y, F(0.0));
        //    ret.m_me[8] = _make_cuComplex(F(2.0) * m_me[8].y, F(0.0));
        //    return ret;
        //}

        /**
        * make any matrix to SU3
        */
        __device__ __inline__ void Norm()
        {
            const Real fDiv1 = __div(F(1.0), _sqrt(__cuCabsSqf(m_me[0]) + __cuCabsSqf(m_me[1])));
            m_me[0] = cuCmulf_cr(m_me[0], fDiv1);
            m_me[1] = cuCmulf_cr(m_me[1], fDiv1);
            m_me[2].x = -m_me[1].x;
            m_me[2].y =  m_me[1].y;
            m_me[3] = _cuConjf(m_me[0]);
        }

        /**
         * HYP projection
         * When Det[U] is large, it need large ite...
         */
        __device__ __inline__ void Proj(BYTE ite = 3)
        {
            //tr(me^+ me) = m0^* m0 + m1^* m1 + m2^* m2 + ...
            //Real fDiv =
            //    __cuCabsSqf(m_me[0])
            //    + __cuCabsSqf(m_me[1])
            //    + __cuCabsSqf(m_me[2])
            //    + __cuCabsSqf(m_me[3]);
            ////1 / sqrt( tr(me^+ me) / 3 )
            //fDiv = __rcp(_sqrt(fDiv * F(0.5)));
            //m_me[0] = cuCmulf_cr(m_me[0], fDiv);
            //m_me[1] = cuCmulf_cr(m_me[1], fDiv);
            //m_me[2] = cuCmulf_cr(m_me[2], fDiv);
            //m_me[3] = cuCmulf_cr(m_me[3], fDiv);

            CLGComplex x[4];
            CLGComplex tmp[2];
            for (BYTE byIt = 0; byIt < ite; ++byIt)
            {
                //x = (-1/2) me^+ me
                x[0] = cuCmulf_cr(__LINE_MULDN_ME2(0, 2, 0, 2), F(-0.5));
                x[1] = cuCmulf_cr(__LINE_MULDN_ME2(0, 2, 1, 3), F(-0.5));
                x[2] = cuCmulf_cr(__LINE_MULDN_ME2(1, 3, 0, 2), F(-0.5));
                x[3] = cuCmulf_cr(__LINE_MULDN_ME2(1, 3, 1, 3), F(-0.5));

                Real fDiv = __rcp(-(x[0].x + x[3].x));
                x[0] = cuCmulf_cr(x[0], fDiv);
                x[1] = cuCmulf_cr(x[1], fDiv);
                x[2] = cuCmulf_cr(x[2], fDiv);
                x[3] = cuCmulf_cr(x[3], fDiv);

                fDiv = _sqrt(fDiv);
                m_me[0] = cuCmulf_cr(m_me[0], fDiv);
                m_me[1] = cuCmulf_cr(m_me[1], fDiv);
                m_me[2] = cuCmulf_cr(m_me[2], fDiv);
                m_me[3] = cuCmulf_cr(m_me[3], fDiv);

                //x += (3/2)
                x[0].x = x[0].x + F(1.5);
                x[3].x = x[3].x + F(1.5);

                //x = me.x
                tmp[0] = _cuCaddf(_cuCmulf(m_me[0], x[0]), _cuCmulf(m_me[1], x[2]));
                tmp[1] = _cuCaddf(_cuCmulf(m_me[0], x[1]), _cuCmulf(m_me[1], x[3]));
                memcpy(m_me, tmp, sizeof(CLGComplex) * 2);

                tmp[0] = _cuCaddf(_cuCmulf(m_me[2], x[0]), _cuCmulf(m_me[3], x[2]));
                tmp[1] = _cuCaddf(_cuCmulf(m_me[2], x[1]), _cuCmulf(m_me[3], x[3]));
                memcpy(m_me + 2, tmp, sizeof(CLGComplex) * 2);

                //me = x
                //coef = det(me)
                const Real fArg = -__cuCargf(Determinent(m_me)) * F(0.5);
                const CLGComplex coef = _make_cuComplex(_cos(fArg), _sin(fArg));
                m_me[0] = _cuCmulf(coef, m_me[0]);
                m_me[1] = _cuCmulf(coef, m_me[1]);
                m_me[2] = _cuCmulf(coef, m_me[2]);
                m_me[3] = _cuCmulf(coef, m_me[3]);
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
        //__device__ __inline__ void CabbiboMarinariProj(/*BYTE ite = 1*/)
        //{
        //    //for (BYTE iteration = 0; iteration < ite; ++iteration)
        //    //{
        //        CLGComplex a11 = _cuCaddf(_cuConjf(m_me[0]), m_me[4]);
        //        CLGComplex b11 = _cuCaddf(_cuConjf(m_me[0]), m_me[8]);
        //        CLGComplex c22 = _cuCaddf(_cuConjf(m_me[4]), m_me[8]);
        //        CLGComplex a12 = _cuCsubf(_cuConjf(m_me[3]), m_me[1]);
        //        CLGComplex b13 = _cuCsubf(_cuConjf(m_me[6]), m_me[2]);
        //        CLGComplex c23 = _cuCsubf(_cuConjf(m_me[7]), m_me[5]);
        //        //CLGComplex a12 = _cuCsubf(_cuConjf(m_me[1]), m_me[3]);
        //        //CLGComplex b13 = _cuCsubf(_cuConjf(m_me[2]), m_me[6]);
        //        //CLGComplex c23 = _cuCsubf(_cuConjf(m_me[5]), m_me[7]);
        //        Real fNorm = __rcp(_sqrt(__cuCabsSqf(a11) + __cuCabsSqf(a12)));
        //        a11 = cuCmulf_cr(a11, fNorm);
        //        a12 = cuCmulf_cr(a12, fNorm);
        //        fNorm = __rcp(_sqrt(__cuCabsSqf(b11) + __cuCabsSqf(b13)));
        //        b11 = cuCmulf_cr(b11, fNorm);
        //        b13 = cuCmulf_cr(b13, fNorm);
        //        fNorm = __rcp(_sqrt(__cuCabsSqf(c22) + __cuCabsSqf(c23)));
        //        c22 = cuCmulf_cr(c22, fNorm);
        //        c23 = cuCmulf_cr(c23, fNorm);

        //        /**
        //         * ({
        //        {a11 b11,
        //        a12 c22 - a11 b13 Conjugate[c23], 
        //        a11 b13 Conjugate[c22] + a12 c23},
        //        {-b11 Conjugate[a12], 
        //        c22 Conjugate[a11] + b13 Conjugate[a12] Conjugate[c23], 
        //        c23 Conjugate[a11] - b13 Conjugate[a12] Conjugate[c22]},
        //        {-Conjugate[b13],
        //        -Conjugate[b11] Conjugate[c23], 
        //        Conjugate[b11] Conjugate[c22]}
        //        })
        //         */
        //        m_me[0] = _cuCmulf(a11, b11);
        //        m_me[1] = _cuCsubf(_cuCmulf(a12, c22), _cuCmulf(_cuCmulf(a11, b13), _cuConjf(c23)));
        //        m_me[2] = _cuCaddf(_cuCmulf(_cuCmulf(a11, b13), _cuConjf(c22)), _cuCmulf(a12, c23));

        //        m_me[3].x = -b11.x * a12.x - b11.y * a12.y;
        //        m_me[3].y = a12.y * b11.x - a12.x * b11.y;
        //        m_me[4] = _cuCaddf(_cuCmulf(c22, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c23)));
        //        m_me[5] = _cuCsubf(_cuCmulf(c23, _cuConjf(a11)), _cuCmulf(_cuCmulf(b13, _cuConjf(a12)), _cuConjf(c22)));

        //        m_me[6].x = -b13.x; m_me[6].y = b13.y;
        //        m_me[7].x = b11.y * c23.y - b11.x * c23.x;
        //        m_me[7].y = b11.x * c23.y + b11.y * c23.x;
        //        m_me[8].x = b11.x * c22.x - b11.y * c22.y;
        //        m_me[8].y = - b11.x * c22.y - b11.y * c22.x;
        //    //}
        //}

        /**
        * U = U + 1
        */
        __device__ __inline__ void AddId()
        {
            m_me[0].x += F(1.0);
            m_me[3].x += F(1.0);
        }

        /**
        * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
        *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
        */
        __device__ __inline__ deviceSU2 Exp(const CLGComplex& a, BYTE uiPrecision) const
        {
            deviceSU2 tmp;

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
                deviceSU2 aUoN = MulCompC(alpha);
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

        __device__ __inline__ deviceSU2 ExpReal(Real a, BYTE uiPrecision) const
        {
            CLGComplex tmp[4];

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
                    CLGComplex tmp2[4];
                    MulRealCArray(tmp2, __div(a, (uiPrecision - i)));
                    ArrayMul(tmp, tmp2);
                }
                tmp[0].x += F(1.0);
                tmp[3].x += F(1.0);
            }
            return deviceSU2(tmp);
        }

        /**
        * For anti hermitian matrix only,
        * 
        * assume M = {{I aa, bb + cc I}, {-bb + cc I, -aa I}} = I (cc sx + bb sy + aa sz)
        * 
        * res = {{Cos[l] + (I aa Sin[l])/l, ((bb + I cc) Sin[l])/l}, {-(((bb - I cc) Sin[l])/l), Cos[l] - (I aa Sin[l])/l}}
        */
        __device__ __inline__ deviceSU2 QuickExp(Real a) const
        {
            const Real aa = a * m_me[0].y;
            const Real bb = a * m_me[1].x;
            const Real cc = a * m_me[1].y;

            const Real L = _sqrt(aa * aa + bb * bb + cc * cc);
            const Real cs = _cos(L);
            const Real sn = _sin(L) / L;
            const Real asn = aa * sn;
            const Real bsn = bb * sn;
            const Real csn = cc * sn;

            deviceSU2 ret;
            ret.m_me[0] = _make_cuComplex(cs, asn);
            ret.m_me[1] = _make_cuComplex(bsn, csn);
            ret.m_me[2] = _make_cuComplex(-bsn, csn);
            ret.m_me[3] = _make_cuComplex(cs, -asn);
            return ret;
        }

#pragma endregion

#pragma region power and log

        /**
         * Eigen values of any 2x2 matrix
         * 
         * x1 = (a + d) / 2
         * x2 = sqrt((a-d)^2 + 4bc) / 2
         * 
         * lambda1 = x1 + x2
         * lambda2 = x1 - x2
         * 
         * v1 = {b, lambda1 - a}
         * v2 = {lambda2 - d, c}
         * 
         */
        //__device__ __inline__ void CalculateEigenValues(CLGComplex& c1, CLGComplex& c2, CLGComplex& c3, UBOOL bDebug = FALSE) const
        __device__ __inline__ void CalculateEigenValues(CLGComplex& c1, CLGComplex& c2, deviceSU2& vectors) const
        {
            CLGComplex x1 = _cuCaddf(m_me[0], m_me[3]);
            x1.x = x1.x * F(0.5);
            x1.y = x1.y * F(0.5);

            CLGComplex x2 = _cuCsubf(m_me[0], m_me[3]);
            CLGComplex x3 = _cuCmulf(m_me[1], m_me[2]);
            x3.x = x3.x * F(4.0);
            x3.y = x3.y * F(4.0);
            x2 = _cuCaddf(_cuCmulf(x2, x2), x3);
            x2 = __cuCsqrtf(x2);
            x2.x = x2.x * F(0.5);
            x2.y = x2.y * F(0.5);

            c1 = _cuCaddf(x1, x2);
            c2 = _cuCsubf(x1, x2);
            
            //v1 = 0, 2
            vectors.m_me[0] = m_me[1];
            vectors.m_me[2] = _cuCsubf(c1, m_me[0]);

            //v2 = 1, 3
            vectors.m_me[1] = _cuCsubf(c2, m_me[3]);
            vectors.m_me[3] = m_me[2];

            const Real fdiv1 = __rcp(_sqrt(__cuCabsSqf(vectors.m_me[0]) + __cuCabsSqf(vectors.m_me[2])));
            const Real fdiv2 = __rcp(_sqrt(__cuCabsSqf(vectors.m_me[1]) + __cuCabsSqf(vectors.m_me[3])));

            vectors.m_me[0].x = vectors.m_me[0].x * fdiv1;
            vectors.m_me[0].y = vectors.m_me[0].y * fdiv1;

            vectors.m_me[2].x = vectors.m_me[2].x * fdiv1;
            vectors.m_me[2].y = vectors.m_me[2].y * fdiv1;

            vectors.m_me[1].x = vectors.m_me[1].x * fdiv2;
            vectors.m_me[1].y = vectors.m_me[1].y * fdiv2;

            vectors.m_me[3].x = vectors.m_me[3].x * fdiv2;
            vectors.m_me[3].y = vectors.m_me[3].y * fdiv2;
        }


        __device__ __inline__ deviceSU2 Power(Real fPower) const
        {
            CLGComplex c1;
            CLGComplex c2;
            deviceSU2 v;
            CalculateEigenValues(c1, c2, v);

            deviceSU2 ret;
            ret.m_me[0] = __cuCpowerf(c1, fPower);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = __cuCpowerf(c2, fPower);
            ret.MulDagger(v);
            ret = v.MulC(ret);
            return ret;
        }

        __device__ __inline__ deviceSU2 Log() const
        {
            CLGComplex c1;
            CLGComplex c2;
            deviceSU2 v;
            CalculateEigenValues(c1, c2, v);

            deviceSU2 ret;
            ret.m_me[0] = __cuClogf(c1);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = __cuClogf(c2);
            ret.MulDagger(v);
            ret = v.MulC(ret);
            return ret;
        }

        //__device__ __inline__ deviceSU3 StrictExp(UBOOL bDebug = FALSE) const
        __device__ __inline__ deviceSU2 StrictExp() const
        {
            CLGComplex c1;
            CLGComplex c2;
            deviceSU2 v;
            CalculateEigenValues(c1, c2, v);

            deviceSU2 ret;
            ret.m_me[0] = __cuCexpf(c1);
            ret.m_me[1] = _zeroc;
            ret.m_me[2] = _zeroc;
            ret.m_me[3] = __cuCexpf(c2);
            ret.MulDagger(v);
            ret = v.MulC(ret);
            return ret;
        }

#pragma endregion

        CLGComplex m_me[4];
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _SU2_H_

//=============================================================================
// END OF FILE
//=============================================================================
