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

#define __LINE_MUL(a, b, c, d, ee, ff) cuCaddf(cuCaddf(cuCmulf(m_me[a], right.m_me[d]), cuCmulf(m_me[b], right.m_me[ee])), cuCmulf(m_me[c], right.m_me[ff]))
// 1.0f / sqrt(3)
#define InvSqrt3 (0.5773502691896258f)
// 2.0f / sqrt(3)
#define InvSqrt3_2 (1.1547005383792517f)

__BEGIN_NAMESPACE

struct CLGAPI deviceSU3
{
    __device__ deviceSU3()
    {
        
    }

    __device__ deviceSU3(const deviceSU3& other)
    {
        memcpy(m_me, other.m_me, sizeof(cuComplex) * 9);
    }

    #pragma region create

    /**
    * ret = 0
    */
    __device__ __inline__ static deviceSU3 makeSU3Zero()
    {
        deviceSU3 ret;
        ret.m_me[0] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[1] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[2] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[3] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[4] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[5] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[6] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[7] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[8] = make_cuComplex(0.0f, 0.0f);
        return ret;
    }

    /**
    * ret = I
    */
    __device__ __inline__ static deviceSU3 makeSU3Id()
    {
        deviceSU3 ret;
        ret.m_me[0] = make_cuComplex(1.0f, 0.0f);
        ret.m_me[1] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[2] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[3] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[4] = make_cuComplex(1.0f, 0.0f);
        ret.m_me[5] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[6] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[7] = make_cuComplex(0.0f, 0.0f);
        ret.m_me[8] = make_cuComplex(1.0f, 0.0f);
        return ret;
    }

    /**
    * can be called only after CLatticeData is created
    * ret = random
    */
    __device__ __inline__ static deviceSU3 makeSU3Random(CRandomSchrage* randomGenerator, UINT fatIndex)
    {
        deviceSU3 ret;
        ret.m_me[0] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[1] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[2] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[3] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[4] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[5] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[6] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[7] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.m_me[8] = make_cuComplex(randomGenerator->_deviceRandomF(fatIndex) - 0.5f, randomGenerator->_deviceRandomF(fatIndex) - 0.5f);
        ret.Norm();
        return ret;
    }

    /**
    * can be called only after CLatticeData is created
    * ret = i r_a T_a, r_a is random real number, T_a are generators
    * r_a T_a = r1 T1 + r2 T2 + ...
    *
    *     r3 +r8/sqrt3     r1-ir2        r4-ir5
    * =   r1+ir2       -r3+r8/sqrt3      r6-ir7
    *     r4+ir5           r6+ir7      -2r8/sqrt3
    *
    */
    __device__ __inline__ static deviceSU3 makeSU3RandomGenerator(CRandomSchrage* randomGenerator, UINT fatIndex)
    {
        FLOAT r1 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r2 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r3 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r4 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r5 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r6 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r7 = randomGenerator->_deviceRandomGaussF(fatIndex);
        FLOAT r8 = randomGenerator->_deviceRandomGaussF(fatIndex);

        deviceSU3 ret;
        ret.m_me[0] = make_cuComplex(0.0f, r3 + r8 * InvSqrt3);
        ret.m_me[1] = make_cuComplex(r2, r1);
        ret.m_me[2] = make_cuComplex(r5, r4);
        ret.m_me[3] = make_cuComplex(-r2, r1);
        ret.m_me[4] = make_cuComplex(0.0f, -r3 + r8 * InvSqrt3);
        ret.m_me[5] = make_cuComplex(r7, r6);
        ret.m_me[6] = make_cuComplex(-r5, r4);
        ret.m_me[7] = make_cuComplex(-r7, r6);
        ret.m_me[8] = make_cuComplex(0.0f, - r8 * InvSqrt3_2);
        return ret;
    }

    #pragma endregion create

    #pragma region operators

    __device__ __inline__ void Add(const deviceSU3& right)
    {
        m_me[0] = cuCaddf(m_me[0], right.m_me[0]);
        m_me[1] = cuCaddf(m_me[1], right.m_me[1]);
        m_me[2] = cuCaddf(m_me[2], right.m_me[2]);
        m_me[3] = cuCaddf(m_me[3], right.m_me[3]);
        m_me[4] = cuCaddf(m_me[4], right.m_me[4]);
        m_me[5] = cuCaddf(m_me[5], right.m_me[5]);
        m_me[6] = cuCaddf(m_me[6], right.m_me[6]);
        m_me[7] = cuCaddf(m_me[7], right.m_me[7]);
        m_me[8] = cuCaddf(m_me[8], right.m_me[8]);
    }

    __device__ __inline__ void Sub(const deviceSU3& right)
    {
        m_me[0] = cuCsubf(m_me[0], right.m_me[0]);
        m_me[1] = cuCsubf(m_me[1], right.m_me[1]);
        m_me[2] = cuCsubf(m_me[2], right.m_me[2]);
        m_me[3] = cuCsubf(m_me[3], right.m_me[3]);
        m_me[4] = cuCsubf(m_me[4], right.m_me[4]);
        m_me[5] = cuCsubf(m_me[5], right.m_me[5]);
        m_me[6] = cuCsubf(m_me[6], right.m_me[6]);
        m_me[7] = cuCsubf(m_me[7], right.m_me[7]);
        m_me[8] = cuCsubf(m_me[8], right.m_me[8]);
    }

    __device__ __inline__ void Add(const cuComplex& right)
    {
        m_me[0] = cuCaddf(m_me[0], right);
        m_me[4] = cuCaddf(m_me[4], right);
        m_me[8] = cuCaddf(m_me[8], right);
    }

    __device__ __inline__ void Sub(const cuComplex& right)
    {
        m_me[0] = cuCsubf(m_me[0], right);
        m_me[4] = cuCsubf(m_me[4], right);
        m_me[8] = cuCsubf(m_me[8], right);
    }

    /**
    * left = left * right
    */
    __device__ __inline__ void Mul(const deviceSU3& right)
    {
        cuComplex res[9];
        res[0] = __LINE_MUL(0, 1, 2, 0, 3, 6);
        res[1] = __LINE_MUL(0, 1, 2, 1, 4, 7);
        res[2] = __LINE_MUL(0, 1, 2, 2, 5, 8);

        res[3] = __LINE_MUL(3, 4, 5, 0, 3, 6);
        res[4] = __LINE_MUL(3, 4, 5, 1, 4, 7);
        res[5] = __LINE_MUL(3, 4, 5, 2, 5, 8);

        res[6] = __LINE_MUL(6, 7, 8, 0, 3, 6);
        res[7] = __LINE_MUL(6, 7, 8, 1, 4, 7);
        res[8] = __LINE_MUL(6, 7, 8, 2, 5, 8);
        memcpy(m_me, res, sizeof(cuComplex) * 9);
    }

    __device__ __inline void Mul(const cuComplex& right)
    {
        m_me[0] = cuCmulf(m_me[0], right);
        m_me[1] = cuCmulf(m_me[1], right);
        m_me[2] = cuCmulf(m_me[2], right);
        m_me[3] = cuCmulf(m_me[3], right);
        m_me[4] = cuCmulf(m_me[4], right);
        m_me[5] = cuCmulf(m_me[5], right);
        m_me[6] = cuCmulf(m_me[6], right);
        m_me[7] = cuCmulf(m_me[7], right);
        m_me[8] = cuCmulf(m_me[8], right);
    }

    /**
    * left = left * right
    */
    __device__ __inline__ void Scale(const cuComplex& right)
    {
        m_me[0] = cuCmulf(m_me[0], right);
        m_me[1] = cuCmulf(m_me[1], right);
        m_me[2] = cuCmulf(m_me[2], right);
        m_me[3] = cuCmulf(m_me[3], right);
        m_me[4] = cuCmulf(m_me[4], right);
        m_me[5] = cuCmulf(m_me[5], right);
        m_me[6] = cuCmulf(m_me[6], right);
        m_me[7] = cuCmulf(m_me[7], right);
        m_me[8] = cuCmulf(m_me[8], right);
    }

    __device__ __inline__ deviceSU3 Addc(const deviceSU3& right) const  { deviceSU3 ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceSU3 Subc(const deviceSU3& right) const  { deviceSU3 ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceSU3 Addc(const cuComplex& right) const  { deviceSU3 ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceSU3 Subc(const cuComplex& right) const  { deviceSU3 ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceSU3 Mulc(const deviceSU3& right) const  { deviceSU3 ret(*this); ret.Mul(right); return ret; }
    __device__ __inline__ deviceSU3 Mulc(const cuComplex& right) const { deviceSU3 ret(*this); ret.Mul(right); return ret; }
    __device__ __inline__ deviceSU3 Scalec(const cuComplex& right) const  { deviceSU3 ret(*this); ret.Scale(right); return ret; }

    #pragma endregion operators

    #pragma region useful functions

    /**
    * ret = ||U||, where U is a d=3 matrix
    * we do NOT need ||SU3|| because ||SU3||=1
    */
    __device__ static __inline__ cuComplex _deviceD3MatrixDeterminent(const cuComplex* u)
    {
        return cuCsubf(
            cuCaddf(
                cuCaddf(
                    cuCmulf(cuCmulf(u[0], u[4]), u[8]),
                    cuCmulf(cuCmulf(u[1], u[5]), u[6])
                ),
                cuCmulf(cuCmulf(u[3], u[7]), u[2])
            ),
            cuCaddf(
                cuCaddf(
                    cuCmulf(cuCmulf(u[2], u[4]), u[6]),
                    cuCmulf(cuCmulf(u[1], u[3]), u[8])
                ),
                cuCmulf(cuCmulf(u[0], u[5]), u[7])
            )
        );
    }

    /**
    * Re[Tr[U]]
    */
    __device__ __inline__ FLOAT ReTr()
    {
        return m_me[0].x + m_me[4].x + m_me[8].x;
    }

    /**
    * res = Conjugate[Transpose[U]]
    */
    __device__ __inline__ void Dagger()
    {
        cuComplex res[9];
        res[0] = cuConjf(m_me[0]);
        res[1] = cuConjf(m_me[3]);
        res[2] = cuConjf(m_me[6]);

        res[3] = cuConjf(m_me[1]);
        res[4] = cuConjf(m_me[4]);
        res[5] = cuConjf(m_me[7]);

        res[6] = cuConjf(m_me[2]);
        res[7] = cuConjf(m_me[5]);
        res[8] = cuConjf(m_me[8]);

        memcpy(m_me, res, sizeof(cuComplex) * 9);
    }

    /**
    * A=(U-U^dagger)/2
    * res = A - tr(A)/3
    */
    __device__ __inline__ void TrTa()
    {
        //new [0], [4], [8] = 2img I / 2
        //we DO NOT calculate it here
        //m_me[0] = make_cuComplex(0.0f, cuCimagf(m_me[0]));
        //m_me[4] = make_cuComplex(0.0f, cuCimagf(m_me[4]));
        //m_me[8] = make_cuComplex(0.0f, cuCimagf(m_me[8]));

        //new [1] = [1] - conj([3])
        //new [3] = [3] - conj([1]) = -conj(new [1])
        cuComplex new1 = cuCsubf(m_me[1], cuConjf(m_me[3]));
        cuComplex new2 = cuCsubf(m_me[2], cuConjf(m_me[6]));
        cuComplex new5 = cuCsubf(m_me[5], cuConjf(m_me[7]));
        m_me[1] = make_cuComplex(0.5f * cuCrealf(new1), 0.5f * cuCimagf(new1));
        m_me[3] = make_cuComplex(-cuCrealf(m_me[1]), cuCimagf(m_me[1]));
        m_me[2] = make_cuComplex(0.5f * cuCrealf(new2), 0.5f * cuCimagf(new2));
        m_me[6] = make_cuComplex(-cuCrealf(m_me[2]), cuCimagf(m_me[2]));
        m_me[5] = make_cuComplex(0.5f * cuCrealf(new5), 0.5f * cuCimagf(new5));
        m_me[7] = make_cuComplex(-cuCrealf(m_me[5]), cuCimagf(m_me[5]));

        FLOAT fImgTr = (m_me[0].y + m_me[4].y + m_me[8].y) / 3.0f;
        m_me[0] = make_cuComplex(0.0f, m_me[0].y - fImgTr);
        m_me[4] = make_cuComplex(0.0f, m_me[4].y - fImgTr);
        m_me[8] = make_cuComplex(0.0f, m_me[8].y - fImgTr);
    }

    /**
    * make any matrix to SU3
    */
    __device__ __inline__ void Norm()
    {
        FLOAT fDiv1 = 1.0f / sqrtf(cuCabsSqf(m_me[0]) + cuCabsSqf(m_me[1]) + cuCabsSqf(m_me[2]));
        m_me[0] = cuCmulf(m_me[0], fDiv1);
        m_me[1] = cuCmulf(m_me[1], fDiv1);
        m_me[2] = cuCmulf(m_me[2], fDiv1);

        //it is the name in Bridge++
        cuComplex sp1 = cuConjf(
            cuCaddf(cuCaddf(
                cuCmulf(m_me[0], cuConjf(m_me[3]))
                , cuCmulf(m_me[1], cuConjf(m_me[4])))
                , cuCmulf(m_me[2], cuConjf(m_me[5])))
        );

        m_me[3] = cuCsubf(m_me[3], cuCmulf(sp1, m_me[0]));
        m_me[4] = cuCsubf(m_me[4], cuCmulf(sp1, m_me[1]));
        m_me[5] = cuCsubf(m_me[5], cuCmulf(sp1, m_me[2]));
        FLOAT fDiv2 = 1.0f / sqrtf(cuCabsSqf(m_me[3]) + cuCabsSqf(m_me[4]) + cuCabsSqf(m_me[5]));
        m_me[3] = cuCmulf(m_me[3], fDiv2);
        m_me[4] = cuCmulf(m_me[4], fDiv2);
        m_me[5] = cuCmulf(m_me[5], fDiv2);

        m_me[6] = cuConjf(cuCsubf(cuCmulf(m_me[1], m_me[5]), cuCmulf(m_me[2], m_me[4])));
        m_me[7] = cuConjf(cuCsubf(cuCmulf(m_me[2], m_me[3]), cuCmulf(m_me[0], m_me[5])));
        m_me[8] = cuConjf(cuCsubf(cuCmulf(m_me[0], m_me[4]), cuCmulf(m_me[1], m_me[3])));
    }

    /**
    * U' = exp(aU) = (1 + a U + a^2 U^2/2 +  ... + a^N U^N/N!)
    *    = 1 + a U (1 + a U /2 (1 + a U/3 ...))
    */
    __device__ __inline__ deviceSU3 Exp(const cuComplex& a, UINT uiPrecision) const
    {
        deviceSU3 identity = makeSU3Id();
        deviceSU3 tmp(*this);

        /**
        * tmp = U
        * loop1: tmp = 1+aU/N
        * loop2: tmp = 1+(aU/(N-1))(1+aU/N)
        * loop3: tmp = 1+(aU/(N-2))(1+aU/(N-1)(1+aU/N))
        * ...
        * loopN: tmp = 1+aU (1 + aU/2 (1+...))
        */
        for (int i = 0; i < uiPrecision; ++i)
        {
            FLOAT exp_factor = 1.0f / (uiPrecision - i);
            cuComplex alpha = cuCmulf(a, exp_factor);
            tmp.Mul(alpha);
            tmp = identity.Addc(tmp);
        }
        tmp.Norm();
        return tmp;
    }

    #pragma endregion useful functions

    cuComplex m_me[9];
};


__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
