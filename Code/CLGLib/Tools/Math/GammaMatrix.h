//=============================================================================
// FILENAME : GammaMatrix.h
// 
// DESCRIPTION:
// This is a gamma matrix implementation from Bridge++
//
// m_uiIndex[row]: index of element in raw which is non-zero (is m_me[row])
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================

#ifndef _GAMMAMATRIX_H_
#define _GAMMAMATRIX_H_

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    struct alignas(64) gammaMatrix
    {
    public:
        __device__ gammaMatrix()
        {
            for (INT i = 0; i < 4; ++i)
            {
                m_uiIndex[i] = 0;
                m_me[i] = make_cuComplexI(0, 0);
            }
        }

        __device__ __inline__ void Set(UINT x, UINT y, INT r, INT i)
        {
            assert(x < 4);
            assert(y < 4);
            m_uiIndex[x] = y;
            m_me[x] = make_cuComplexI(r, i);
        }

        /**
        * ret = (Gamma . Spinor)
        */
        __device__ __inline__ deviceWilsonVectorSU3 MulWilsonC(const deviceWilsonVectorSU3& other) const
        {
            deviceWilsonVectorSU3 ret = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
            ret.m_d[0] = other.m_d[m_uiIndex[0]];
            ret.m_d[0].MulCompI(m_me[0]);

            ret.m_d[1] = other.m_d[m_uiIndex[1]];
            ret.m_d[1].MulCompI(m_me[1]);

            ret.m_d[2] = other.m_d[m_uiIndex[2]];
            ret.m_d[2].MulCompI(m_me[2]);

            ret.m_d[3] = other.m_d[m_uiIndex[3]];
            ret.m_d[3].MulCompI(m_me[3]);

            return ret;
        }

        UINT m_uiIndex[4];
        cuComplexI m_me[4];

        /*
        * We are on device. Do not use thie function...
        * this is on device, so use print
        */
        __device__ __inline__ void Print()
        {
            for (int row = 0; row < 4; ++row)
            {
                if (m_uiIndex[row] == 0)
                {
                    printf("(%2d,%2d) 0     0     0\n",
                        m_me[row].x, m_me[row].y);
                }
                else if (m_uiIndex[row] == 1)
                {
                    printf("0     (%2d,%2d) 0     0\n",
                        m_me[row].x, m_me[row].y);
                }
                else if (m_uiIndex[row] == 2)
                {
                    printf("0     0     (%2d,%2d) 0     0\n",
                        m_me[row].x, m_me[row].y);
                }
                else if (m_uiIndex[row] == 3)
                {
                    printf("0     0     0     (%2d,%2d)\n",
                        m_me[row].x, m_me[row].y);
                }
            }
        }

        gammaMatrix __device__ __inline__ Mult(INT n) const
        {
            gammaMatrix ret;

            for (int row = 0; row < 4; ++row)
            {
                ret.m_uiIndex[row] = m_uiIndex[row];
                ret.m_me[row] = cuCmulI(m_me[row], make_cuComplexI(n, 0));
            }

            return ret;
        }

        friend class gammaMatrixSet;

    private:

        /**
        * For those have (g1 * g2)__{ab} = g1_{ac} g2_{cb}
        */
        gammaMatrix __device__ __inline__ _mult(const gammaMatrix& other) const
        {
            gammaMatrix ret;

            for (int row = 0; row < 4; ++row)
            {
                ret.m_uiIndex[row] = other.m_uiIndex[m_uiIndex[row]];
                ret.m_me[row] = cuCmulI(m_me[row], other.m_me[m_uiIndex[row]]);
            }

            return ret;
        }

        /**
        * For those have (g1 * g2)__{ab} = i g1_{ac} g2_{cb}
        */
        gammaMatrix __device__ __inline__ _mult_i(const gammaMatrix& other) const
        {
            gammaMatrix ret;

            for (int row = 0; row < 4; ++row)
            {
                ret.m_uiIndex[row] = other.m_uiIndex[m_uiIndex[row]];
                ret.m_me[row] = cuCmulI(make_cuComplexI(0, 1), cuCmulI(m_me[row], other.m_me[m_uiIndex[row]]));
            }

            return ret;
        }
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

enum EGammaMatrix
{
    ZERO,
    UNITY,
    GAMMA1,
    GAMMA2,
    GAMMA3,
    GAMMA4,
    GAMMA5,
    GAMMA51,
    GAMMA52,
    GAMMA53,
    GAMMA54,
    GAMMA15,
    GAMMA25,
    GAMMA35,
    GAMMA45,
    SIGMA12,
    SIGMA23,
    SIGMA31,
    SIGMA41,
    SIGMA42,
    SIGMA43,
    CHARGECONJG,
    EGM_MAX,
};

enum EGammaMatrixSet
{
    EGMS_Dirac,
    EGMS_Chiral,
};

class gammaMatrixSet
{
public:
    __device__ gammaMatrixSet(EGammaMatrixSet eSet)
    {
        m_gm[UNITY].Set(0, 0, 1, 0);
        m_gm[UNITY].Set(1, 1, 1, 0);
        m_gm[UNITY].Set(2, 2, 1, 0);
        m_gm[UNITY].Set(3, 3, 1, 0);

        if (EGMS_Dirac == eSet)
        {
            m_gm[GAMMA1].Set(0, 3, 0, -1);
            m_gm[GAMMA1].Set(1, 2, 0, -1);
            m_gm[GAMMA1].Set(2, 1, 0, 1);
            m_gm[GAMMA1].Set(3, 0, 0, 1);

            m_gm[GAMMA2].Set(0, 3, -1, 0);
            m_gm[GAMMA2].Set(1, 2, 1, 0);
            m_gm[GAMMA2].Set(2, 1, 1, 0);
            m_gm[GAMMA2].Set(3, 0, -1, 0);

            m_gm[GAMMA3].Set(0, 2, 0, -1);
            m_gm[GAMMA3].Set(1, 3, 0, 1);
            m_gm[GAMMA3].Set(2, 0, 0, 1);
            m_gm[GAMMA3].Set(3, 1, 0, -1);

            m_gm[GAMMA4].Set(0, 0, 1, 0);
            m_gm[GAMMA4].Set(1, 1, 1, 0);
            m_gm[GAMMA4].Set(2, 2, -1, 0);
            m_gm[GAMMA4].Set(3, 3, -1, 0);

            m_gm[GAMMA5].Set(0, 2, 1, 0);
            m_gm[GAMMA5].Set(1, 3, 1, 0);
            m_gm[GAMMA5].Set(2, 0, 1, 0);
            m_gm[GAMMA5].Set(3, 1, 1, 0);
        }
        else 
        {
            m_gm[GAMMA1].Set(0, 3, 0, -1);
            m_gm[GAMMA1].Set(1, 2, 0, -1);
            m_gm[GAMMA1].Set(2, 1, 0, 1);
            m_gm[GAMMA1].Set(3, 0, 0, 1);

            m_gm[GAMMA2].Set(0, 3, -1, 0);
            m_gm[GAMMA2].Set(1, 2, 1, 0);
            m_gm[GAMMA2].Set(2, 1, 1, 0);
            m_gm[GAMMA2].Set(3, 0, -1, 0);

            m_gm[GAMMA3].Set(0, 2, 0, -1);
            m_gm[GAMMA3].Set(1, 3, 0, 1);
            m_gm[GAMMA3].Set(2, 0, 0, 1);
            m_gm[GAMMA3].Set(3, 1, 0, -1);

            m_gm[GAMMA4].Set(0, 2, -1, 0);
            m_gm[GAMMA4].Set(1, 3, -1, 0);
            m_gm[GAMMA4].Set(2, 0, -1, 0);
            m_gm[GAMMA4].Set(3, 1, -1, 0);

            m_gm[GAMMA5].Set(0, 0, 1, 0);
            m_gm[GAMMA5].Set(1, 1, 1, 0);
            m_gm[GAMMA5].Set(2, 2, -1, 0);
            m_gm[GAMMA5].Set(3, 3, -1, 0);
        }

        m_gm[GAMMA51] = m_gm[GAMMA5]._mult(m_gm[GAMMA1]);
        m_gm[GAMMA52] = m_gm[GAMMA5]._mult(m_gm[GAMMA2]);
        m_gm[GAMMA53] = m_gm[GAMMA5]._mult(m_gm[GAMMA3]);
        m_gm[GAMMA54] = m_gm[GAMMA5]._mult(m_gm[GAMMA4]);

        m_gm[GAMMA15] = m_gm[GAMMA1]._mult(m_gm[GAMMA5]);
        m_gm[GAMMA25] = m_gm[GAMMA2]._mult(m_gm[GAMMA5]);
        m_gm[GAMMA35] = m_gm[GAMMA3]._mult(m_gm[GAMMA5]);
        m_gm[GAMMA45] = m_gm[GAMMA4]._mult(m_gm[GAMMA5]);

        m_gm[SIGMA12] = m_gm[GAMMA2]._mult_i(m_gm[GAMMA1]);
        m_gm[SIGMA23] = m_gm[GAMMA3]._mult_i(m_gm[GAMMA2]);
        m_gm[SIGMA31] = m_gm[GAMMA1]._mult_i(m_gm[GAMMA3]);

        m_gm[SIGMA41] = m_gm[GAMMA1]._mult_i(m_gm[GAMMA4]);
        m_gm[SIGMA42] = m_gm[GAMMA2]._mult_i(m_gm[GAMMA4]);
        m_gm[SIGMA43] = m_gm[GAMMA3]._mult_i(m_gm[GAMMA4]);

        m_gm[CHARGECONJG] = m_gm[GAMMA4]._mult(m_gm[GAMMA2]);
    }

    __device__ void Print()
    {
        for (INT i = 0; i < EGM_MAX; ++i)
        {
            m_gm[i].Print();
        }
    }

    gammaMatrix m_gm[EGM_MAX];
};

__END_NAMESPACE

#endif //#ifndef _GAMMAMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================
