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

__BEGIN_NAMESPACE

struct CLGAPI deviceSU3
{
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
        return cuCrealf(cuCmulf(cuCmulf(m_me[0], m_me[4]), m_me[8]));
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
    * make any matrix to SU3
    */
    void __device__ __inline__ Norm()
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

    #pragma endregion useful functions

    cuComplex m_me[9];
};


__END_NAMESPACE

#endif //#ifndef _SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
