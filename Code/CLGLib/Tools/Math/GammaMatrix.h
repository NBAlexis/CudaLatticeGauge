//=============================================================================
// FILENAME : GammaMatrix.h
// 
// DESCRIPTION:
//
// If define as 0, it is gamma matrix implementation in Bridge++, each cost 64 bytes for float and 128 for double
//
// If define as 1, it is:
// The Gamma matrix is a USHORT divded into 8 parts (2 bits for each), we use UINT to align it to 4 bytes
//
// The data is 0000000000000000 | a_1 | a_2 | a_3 | a_4 | b_1 | b_2 | b_3 | b_4
// Where a_i and b_i = (0, 1, 2, 3)
//
// The a_i element of row(i) is none-zero, the value is b_i element of Z4
// 
// The single float version, 0 is 8% faster
// Not tested, but when double float, 1 might be better
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================

#ifndef _GAMMAMATRIX_H_
#define _GAMMAMATRIX_H_

//define as 0 is 8% faster for single float
#define __GAMMA_IMPLEMENTATION_SWITCH 0

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#if __GAMMA_IMPLEMENTATION_SWITCH

    struct alignas(4) gammaMatrix
    {
        //Initialize is not allowed for constant variables
        __device__ gammaMatrix()
        {

        }
        __device__ gammaMatrix(UINT uiValue) : m_uiValue(uiValue) { }
        __device__ gammaMatrix(const gammaMatrix& other) : m_uiValue(other.m_uiValue) { }


        __device__ void Reset()
        {
            m_uiValue = 0;
        }

        /**
        * b = 0, 1, 2, 3
        * for 1, i, -1, -i
        */
        __device__ __inline__ static deviceSU3Vector MultZ4(const deviceSU3Vector& v, UINT uiB)
        {
            // real : 0,1,2,3  -> 1, 0, -1, 0
            // imag : 0,1,2,3  -> 0, 1, 0, -1
            return v.MulCompC(_make_cuComplex(3 == uiB ? F(0.0) : (F(1.0) - uiB), 0 == uiB ? F(0.0) : (F(2.0) - uiB)));
        }

        __device__ __inline__ void Set(UINT x, UINT y, UINT b)
        {
            assert(x < 4);
            assert(y < 4);
            assert(b < 4);

            UINT uiShift = x << 1;
            UINT uiIndex = y << (uiShift + 8);
            UINT uiV = b << uiShift;
            m_uiValue = m_uiValue | uiIndex | uiV;
        }

        /**
        * ret = (Gamma . Spinor)
        */
        __device__ __inline__ deviceWilsonVectorSU3 MulWilsonC(const deviceWilsonVectorSU3& other) const
        {
            deviceWilsonVectorSU3 ret;
            UINT uiIndexOfRow1 = (m_uiValue >> 8) & 3;
            UINT uiIndexOfRow2 = (m_uiValue >> 10) & 3;
            UINT uiIndexOfRow3 = (m_uiValue >> 12) & 3;
            UINT uiIndexOfRow4 = (m_uiValue >> 14) & 3;
            UINT uiZ4OfRow1 = (m_uiValue/*>>0*/) & 3;
            UINT uiZ4OfRow2 = (m_uiValue >> 2) & 3;
            UINT uiZ4OfRow3 = (m_uiValue >> 4) & 3;
            UINT uiZ4OfRow4 = (m_uiValue >> 6) & 3;

            ret.m_d[0] = MultZ4(other.m_d[uiIndexOfRow1], uiZ4OfRow1);
            ret.m_d[1] = MultZ4(other.m_d[uiIndexOfRow2], uiZ4OfRow2);
            ret.m_d[2] = MultZ4(other.m_d[uiIndexOfRow3], uiZ4OfRow3);
            ret.m_d[3] = MultZ4(other.m_d[uiIndexOfRow4], uiZ4OfRow4);

            return ret;
        }

        /*
        * We are on device. Do not use thie function...
        * this is on device, so use print
        */
        __device__ __inline__ void Print() const
        {
            INT reals[4];
            INT imags[4];
            reals[0] = 1;
            reals[1] = 0;
            reals[2] = -1;
            reals[3] = 0;
            imags[0] = 0;
            imags[1] = 1;
            imags[2] = 0;
            imags[3] = -1;
            for (int row = 0; row < 4; ++row)
            {
                UINT uiShift = (row << 1);
                UINT uiRow = (m_uiValue >> (8 + uiShift)) & 3;
                UINT uiB = (m_uiValue >> uiShift) & 3;

                if (uiRow == 0)
                {
                    printf("(%2d,%2d) 0     0     0\n",
                        reals[uiB], imags[uiB]);
                }
                else if (uiRow == 1)
                {
                    printf("0     (%2d,%2d) 0     0\n",
                        reals[uiB], imags[uiB]);
                }
                else if (uiRow == 2)
                {
                    printf("0     0     (%2d,%2d) 0     0\n",
                        reals[uiB], imags[uiB]);
                }
                else if (uiRow == 3)
                {
                    printf("0     0     0     (%2d,%2d)\n",
                        reals[uiB], imags[uiB]);
                }
            }
        }

        friend class gammaMatrixSet;

    private:

        /**
        * For those have (g1 * g2)__{ab} = g1_{ac} g2_{cb}
        */
        __device__ __inline__ gammaMatrix _mult(const gammaMatrix& other) const
        {
            UINT uiNew = 0;

            for (int row = 0; row < 4; ++row)
            {
                UINT uiShiftMe = (row << 1);
                UINT uiRowMe = (m_uiValue >> (8 + uiShiftMe)) & 3;
                UINT uiBMe = (m_uiValue >> uiShiftMe) & 3;

                UINT uiShiftOther = uiRowMe << 1;
                UINT uiRowOther = (other.m_uiValue >> (8 + uiShiftOther)) & 3;
                UINT uiBOther = (other.m_uiValue >> uiShiftOther) & 3;

                uiBOther = (uiBOther + uiBMe) & 3;

                UINT uiIndexNew = uiRowOther << (uiShiftMe + 8);
                UINT uiBNew = uiBOther << uiShiftMe;
                uiNew = uiNew | uiIndexNew | uiBNew;
            }

            return gammaMatrix(uiNew);
        }

        /**
        * For those have (g1 * g2)__{ab} = i g1_{ac} g2_{cb}
        */
        __device__ __inline__ gammaMatrix _mult_i(const gammaMatrix& other) const
        {
            UINT uiNew = 0;

            for (int row = 0; row < 4; ++row)
            {
                UINT uiShiftMe = (row << 1);
                UINT uiRowMe = (m_uiValue >> (8 + uiShiftMe)) & 3;
                UINT uiBMe = (m_uiValue >> uiShiftMe) & 3;

                UINT uiShiftOther = uiRowMe << 1;
                UINT uiRowOther = (other.m_uiValue >> (8 + uiShiftOther)) & 3;
                UINT uiBOther = (other.m_uiValue >> uiShiftOther) & 3;

                uiBOther = (uiBOther + uiBMe + 1) & 3;

                UINT uiIndexNew = uiRowOther << (uiShiftMe + 8);
                UINT uiBNew = uiBOther << uiShiftMe;
                uiNew = uiNew | uiIndexNew | uiBNew;
            }

            return gammaMatrix(uiNew);
        }

        UINT m_uiValue;
    };

#else

//#if _CLG_DOUBLEFLOAT
//#define __GAMMA_ALIGN 128
//#else
//#define __GAMMA_ALIGN 64
//#endif

    struct alignas(16) gammaMatrix
    {
        //Initialize is not allowed for constant variables
        __device__ gammaMatrix() 
        { 

        }

        __device__ void Reset() { }

        __device__ __inline__ void Set(BYTE x, BYTE y, BYTE b)
        {
            assert(x < 4);
            assert(y < 4);
            assert(b < 4);

            m_uiIndex[x] = y;
            m_byZ4[x] = b;
            //m_cV[x] = _make_cuComplex(3 == b ? F(0.0) : (F(1.0) - b), 0 == b ? F(0.0) : (F(2.0) - b));
        }

        /**
        * ret = (Gamma . Spinor)
        */
        __device__ __inline__ deviceWilsonVectorSU3 MulWilsonC(const deviceWilsonVectorSU3& other) const
        {
            deviceWilsonVectorSU3 ret;

            //ret.m_d[0] = other.m_d[m_uiIndex[0]].MulCompC(m_cV[0]);
            //ret.m_d[1] = other.m_d[m_uiIndex[1]].MulCompC(m_cV[1]);
            //ret.m_d[2] = other.m_d[m_uiIndex[2]].MulCompC(m_cV[2]);
            //ret.m_d[3] = other.m_d[m_uiIndex[3]].MulCompC(m_cV[3]);
            ret.m_d[0] = other.m_d[m_uiIndex[0]].MulZ4C(m_byZ4[0]);
            ret.m_d[1] = other.m_d[m_uiIndex[1]].MulZ4C(m_byZ4[1]);
            ret.m_d[2] = other.m_d[m_uiIndex[2]].MulZ4C(m_byZ4[2]);
            ret.m_d[3] = other.m_d[m_uiIndex[3]].MulZ4C(m_byZ4[3]);

            return ret;
        }

        /*
        * We are on device. Do not use thie function...
        * this is on device, so use print
        */
        __device__ __inline__ void Print() const
        {
            #pragma unroll
            for (BYTE row = 0; row < 4; ++row)
            {
                const CLGComplex cv = _make_cuComplex(3 == m_byZ4[row] ? F(0.0) : (F(1.0) - m_byZ4[row]), 0 == m_byZ4[row] ? F(0.0) : (F(2.0) - m_byZ4[row]));
                const BYTE byNoneZero = m_uiIndex[row];
                if (byNoneZero == 0)
                {
                    printf("(%2f,%2f) 0     0     0\n",
                        cv.x, cv.y);
                }
                else if (byNoneZero == 1)
                {
                    printf("0     (%2f,%2f) 0     0\n",
                        cv.x, cv.y);
                }
                else if (byNoneZero == 2)
                {
                    printf("0     0     (%2f,%2f) 0\n",
                        cv.x, cv.y);
                }
                else if (byNoneZero == 3)
                {
                    printf("0     0     0     (%2f,%2f)\n",
                        cv.x, cv.y);
                }
            }
        }

        friend class gammaMatrixSet;

    private:

        /**
        * For those have (g1 * g2)__{ab} = g1_{ac} g2_{cb}
        */
        __device__ __inline__ gammaMatrix _mult(const gammaMatrix& other) const
        {
            gammaMatrix ret;

            #pragma unroll
            for (int row = 0; row < 4; ++row)
            {
                ret.m_uiIndex[row] = other.m_uiIndex[m_uiIndex[row]];
                ret.m_byZ4[row] = (m_byZ4[row] + other.m_byZ4[m_uiIndex[row]]) & 0x03;
                //ret.m_cV[row] = _cuCmulf(m_cV[row], other.m_cV[m_uiIndex[row]]);
            }

            return ret;
        }

        /**
        * For those have (g1 * g2)__{ab} = i g1_{ac} g2_{cb}
        */
        __device__ __inline__ gammaMatrix _mult_i(const gammaMatrix& other) const
        {
            gammaMatrix ret;

            #pragma unroll
            for (int row = 0; row < 4; ++row)
            {
                ret.m_uiIndex[row] = other.m_uiIndex[m_uiIndex[row]];
                ret.m_byZ4[row] = (m_byZ4[row] + other.m_byZ4[m_uiIndex[row]] + 1) & 0x03;
                //ret.m_cV[row] = _cuCmulf(_make_cuComplex(F(0.0), F(1.0)), _cuCmulf(m_cV[row], other.m_cV[m_uiIndex[row]]));
            }

            return ret;
        }

        /**
        * Note that, except for Unitary and Gamma5, SIGMA12 and SIGMA43 matrix
        *
        * When m_uiIndex[0] = 1, we must have m_uiIndex[1] = 0.
        * at the same time, m_uiIndex[2] = 3, m_uiIndex[3] = 2.
        *
        * When m_uiIndex[0] = 2, we must have m_uiIndex[2] = 0.
        * at the same time, m_uiIndex[1] = 3, m_uiIndex[3] = 1.
        *
        * So, we only need to know m_uiIndex[0], i, we known every m_uiIndex[]
        *
        * for m_uiIndex[0] = 0, m_uiIndex[1] = 1, m_uiIndex[2] = 2, m_uiIndex[3] = 3
        * for m_uiIndex[0] = m and i = n, m_uiIndex[0] = m, m_uiIndex[m] = 0, m_uiIndex[i] = n(k!=0,m,n) m_uiIndex[n] = i
        *
        */
        __device__ __inline__ void CalculateNextSymmetricIndex()
        {
            if (0 == m_uiIndex[0])
            {
                m_byNextSymmetricIndex = 0;
                //m_uiIndex[0] = 0, m_uiIndex[1] = 1, m_uiIndex[2] = 2, m_uiIndex[3] = 3
            }
            else if (1 == m_uiIndex[0])
            {
                m_byNextSymmetricIndex = 2;
                //m_uiIndex[0] = 1, m_uiIndex[1] = 0, m_uiIndex[2] = 3, m_uiIndex[3] = 2
            }
            //two cases
            //m_uiIndex[0] = 2, m_uiIndex[1] = 3, m_uiIndex[2] = 0, m_uiIndex[3] = 1
            //or
            //m_uiIndex[0] = 3, m_uiIndex[1] = 2, m_uiIndex[3] = 0, m_uiIndex[2] = 1
            //Both can be set as m_byNextSymmetricIndex = 1, but we also can distinguish them
            else if (2 == m_uiIndex[0])
            {
                m_byNextSymmetricIndex = 3; //can be both 1 or 3
                //m_uiIndex[0] = 2, m_uiIndex[1] = 3, m_uiIndex[2] = 0, m_uiIndex[3] = 1
            }
            else
            {
                m_byNextSymmetricIndex = 1; //can be both 1 or 2
                //m_uiIndex[0] = 3, m_uiIndex[1] = 2, m_uiIndex[3] = 0, m_uiIndex[2] = 1
            }
        }

    public:

        //CLGComplex m_cV[4];
        //For Z4, we store b
        //i^b, b=0,1,2,3
        //(i^b)* = i ^ (((b & 1) << 1) + b) & 3
        //(i^b1 * i^b2) = i ^ (b1 + b2) & 3
        //CLGComplex = _make_cuComplex(3 == b ? F(0.0) : (F(1.0) - b), 0 == b ? F(0.0) : (F(2.0) - b));
        BYTE m_byZ4[4];
        BYTE m_uiIndex[4];
        BYTE m_byNextSymmetricIndex;
    };

#endif

#if defined(__cplusplus)
}
#endif /* __cplusplus */

enum EGammaMatrixSet
{
    EGMS_Dirac,
    EGMS_Chiral,
};

class gammaMatrixSet 
{
public:
    __device__ static void CreateGammaMatrix(EGammaMatrixSet eSet, gammaMatrix* gmarray)
    {
        //Gamma matrix is not initialized
        #pragma unroll
        for (UINT i = 0; i < EGM_MAX; ++i)
        {
            gmarray[i].Reset();
        }
        gmarray[UNITY].Set(0, 0, 0);
        gmarray[UNITY].Set(1, 1, 0);
        gmarray[UNITY].Set(2, 2, 0);
        gmarray[UNITY].Set(3, 3, 0);

        if (EGMS_Dirac == eSet)
        {
            gmarray[GAMMA1].Set(0, 3, 3);
            gmarray[GAMMA1].Set(1, 2, 3);
            gmarray[GAMMA1].Set(2, 1, 1);
            gmarray[GAMMA1].Set(3, 0, 1);

            gmarray[GAMMA2].Set(0, 3, 2);
            gmarray[GAMMA2].Set(1, 2, 0);
            gmarray[GAMMA2].Set(2, 1, 0);
            gmarray[GAMMA2].Set(3, 0, 2);

            gmarray[GAMMA3].Set(0, 2, 3);
            gmarray[GAMMA3].Set(1, 3, 1);
            gmarray[GAMMA3].Set(2, 0, 1);
            gmarray[GAMMA3].Set(3, 1, 3);

            gmarray[GAMMA4].Set(0, 0, 0);
            gmarray[GAMMA4].Set(1, 1, 0);
            gmarray[GAMMA4].Set(2, 2, 2);
            gmarray[GAMMA4].Set(3, 3, 2);

            gmarray[GAMMA5].Set(0, 2, 0);
            gmarray[GAMMA5].Set(1, 3, 0);
            gmarray[GAMMA5].Set(2, 0, 0);
            gmarray[GAMMA5].Set(3, 1, 0);
        }
        else
        {
            gmarray[GAMMA1].Set(0, 3, 3);
            gmarray[GAMMA1].Set(1, 2, 3);
            gmarray[GAMMA1].Set(2, 1, 1);
            gmarray[GAMMA1].Set(3, 0, 1);

            gmarray[GAMMA2].Set(0, 3, 2);
            gmarray[GAMMA2].Set(1, 2, 0);
            gmarray[GAMMA2].Set(2, 1, 0);
            gmarray[GAMMA2].Set(3, 0, 2);

            gmarray[GAMMA3].Set(0, 2, 3);
            gmarray[GAMMA3].Set(1, 3, 1);
            gmarray[GAMMA3].Set(2, 0, 1);
            gmarray[GAMMA3].Set(3, 1, 3);

            //this is from bridge++, however, the Ref.~book2010, shows it is incorrect.
            //gmarray[GAMMA4].Set(0, 2, 2);
            //gmarray[GAMMA4].Set(1, 3, 2);
            //gmarray[GAMMA4].Set(2, 0, 2);
            //gmarray[GAMMA4].Set(3, 1, 2);
            gmarray[GAMMA4].Set(0, 2, 0);
            gmarray[GAMMA4].Set(1, 3, 0);
            gmarray[GAMMA4].Set(2, 0, 0);
            gmarray[GAMMA4].Set(3, 1, 0);

            gmarray[GAMMA5].Set(0, 0, 0);
            gmarray[GAMMA5].Set(1, 1, 0);
            gmarray[GAMMA5].Set(2, 2, 2);
            gmarray[GAMMA5].Set(3, 3, 2);
        }

        gmarray[GAMMA51] = gmarray[GAMMA5]._mult(gmarray[GAMMA1]);
        gmarray[GAMMA52] = gmarray[GAMMA5]._mult(gmarray[GAMMA2]);
        gmarray[GAMMA53] = gmarray[GAMMA5]._mult(gmarray[GAMMA3]);
        gmarray[GAMMA54] = gmarray[GAMMA5]._mult(gmarray[GAMMA4]);

        gmarray[GAMMA15] = gmarray[GAMMA1]._mult(gmarray[GAMMA5]);
        gmarray[GAMMA25] = gmarray[GAMMA2]._mult(gmarray[GAMMA5]);
        gmarray[GAMMA35] = gmarray[GAMMA3]._mult(gmarray[GAMMA5]);
        gmarray[GAMMA45] = gmarray[GAMMA4]._mult(gmarray[GAMMA5]);

        gmarray[SIGMA12] = gmarray[GAMMA2]._mult_i(gmarray[GAMMA1]);
        gmarray[SIGMA23] = gmarray[GAMMA3]._mult_i(gmarray[GAMMA2]);
        gmarray[SIGMA31] = gmarray[GAMMA1]._mult_i(gmarray[GAMMA3]);

        gmarray[SIGMA41] = gmarray[GAMMA1]._mult_i(gmarray[GAMMA4]);
        gmarray[SIGMA42] = gmarray[GAMMA2]._mult_i(gmarray[GAMMA4]);
        gmarray[SIGMA43] = gmarray[GAMMA3]._mult_i(gmarray[GAMMA4]);

        gmarray[CHARGECONJG] = gmarray[GAMMA4]._mult(gmarray[GAMMA2]);

        #pragma unroll
        for (UINT i = 0; i < EGM_MAX; ++i)
        {
            gmarray[i].CalculateNextSymmetricIndex();
        }
    }
};


__END_NAMESPACE

#endif //#ifndef _GAMMAMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================
