//=============================================================================
// FILENAME : CRationalApproximation.h
// 
// DESCRIPTION:
// We do not calculate it directly, but just read from file
//
//
// REVISION:
//  [05/06/2020 nbale]
//=============================================================================

#ifndef _CRATIONALAPPROXIMATION_H_
#define _CRATIONALAPPROXIMATION_H_

__BEGIN_NAMESPACE

/**
 * In Mathematica 12.0
 *
 * Call: (see OtherProgram/RemesExample.nb )
 *
 * << FunctionApproximations`;
 * res1 = Apart[MiniMaxApproximation[1/Sqrt[x], {x, {0.003, 1}, 3, 3}][[2]][[1]]]
 * Join[{Part[res1, 1]}, Table[Numerator[Part[res1, n]], {n, 2, 4}], Table[Denominator[Part[res1, n]] /. x -> 0, {n, 2, 4}]]
 *
 * "1/Sqrt[x]" is the function to be approximated
 * "{0.003, 1}" is the region
 * "3,3" is the degree
 * The final result in the above case is
 * {0.39046, 0.0511094, 0.140829, 0.596485, 0.00127792, 0.0286165, 0.4106}
 *
 * means {c, a1, a2, a3, b1, b2, b3} such that
 *
 * f(x) = c + sum _i ai / (x + bi)
 *
 * "3" is uiDegree and {c, a1, a2, a3, b1, b2, b3} is the parameters
 */
class CLGAPI CRatinalApproximation
{
public:

    CRatinalApproximation(UINT uiDegree, TArray<Real> parameters)
        : m_uiDegree(uiDegree)
        , m_fC(parameters[0])
    {
        for (UINT i = 0; i < uiDegree; ++i)
        {
            m_lstA.AddItem(parameters[1 + i]);
            m_lstB.AddItem(parameters[1 + uiDegree + i]);
        }
    }

    ~CRatinalApproximation()
    {

    }

    /**
     * Test Function
     */
    Real fx(Real x)
    {
        Real fRet = m_fC;
        for (UINT i = 0; i < m_uiDegree; ++i)
        {
            fRet += m_lstA[i] / (x + m_lstB[i]);
        }
        return fRet;
    }

    UINT m_uiDegree;
    Real m_fC;
    TArray<Real> m_lstA;
    TArray<Real> m_lstB;
};

__END_NAMESPACE

#endif //#ifndef _CRATIONALAPPROXIMATION_H_

//=============================================================================
// END OF FILE
//=============================================================================
