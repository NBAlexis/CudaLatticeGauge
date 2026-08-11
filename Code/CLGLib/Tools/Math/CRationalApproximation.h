//=============================================================================
// FILENAME : CRationalApproximation.h
// 
// DESCRIPTION:
// We do not calculate it directly, but just read from file
//
//
// REVISION:
//  [mm/dd/yy]
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
    CRatinalApproximation() 
        : m_uiDegree(0)
        , m_fC(F(0.0))
        , m_pDeviceData(NULL)
    {
        
    }
    CRatinalApproximation(const TArray<Real>& parameters);
    ~CRatinalApproximation();

    void Initial(const TArray<Real>& parameters);

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

    inline UBOOL operator==(const CRatinalApproximation& Other) const
    {
        if (m_uiDegree != Other.m_uiDegree)
        {
            return FALSE;
        }
        if (appAbs(m_fC - Other.m_fC) > _CLG_FLT_MIN_)
        {
            return FALSE;
        }
        if (m_lstA.Num() != m_lstB.Num())
        {
            return FALSE;
        }
        for (INT i = 0; i < m_lstA.Num(); ++i)
        {
            if (appAbs(m_lstA[i] - Other.m_lstA[i]) > _CLG_FLT_MIN_)
            {
                return FALSE;
            }
            if (appAbs(m_lstB[i] - Other.m_lstB[i]) > _CLG_FLT_MIN_)
            {
                return FALSE;
            }
        }
        return TRUE;
    }

    UINT m_uiDegree;
    Real m_fC;
    TArray<Real> m_lstA;
    TArray<Real> m_lstB;
    Real* m_pDeviceData;
};

/**
* avoid copy of rational approximation parameters when copy the fields
*/
class CRatinalApproximationSet
{
public:
    void Quit()
    {
        for (INT i = 0; i < m_pRASet.Num(); ++i)
        {
            appSafeDelete(m_pRASet[i]);
        }
        m_pRASet.RemoveAll();
    }

    INT Add(const TArray<Real>& ra)
    {
        INT ret = m_pRASet.Num();
        m_pRASet.AddItem(new CRatinalApproximation(ra));
        return ret;
    }

    TArray<CRatinalApproximation*> m_pRASet;
};

extern CLGAPI CRatinalApproximationSet GRASet;

__END_NAMESPACE

#endif //#ifndef _CRATIONALAPPROXIMATION_H_

//=============================================================================
// END OF FILE
//=============================================================================
