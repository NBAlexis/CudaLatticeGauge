//=============================================================================
// FILENAME : CIntegratorMultiLevelOmelyan.cpp
// 
// DESCRIPTION:
// This is the Approximate force gradient integrator for HMC
//
// REVISION:
//  [08/18/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorMultiLevelOmelyan)

void CIntegratorMultiLevelOmelyan::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CMultiLevelNestedIntegrator::Initial(pOwner, pLattice, params);
    m_f2Lambda = OmelyanLambda2;
    if (!params.FetchValueReal(_T("Omelyan2Lambda"), m_f2Lambda))
    {
        m_f2Lambda = OmelyanLambda2;
    }
}

void CIntegratorMultiLevelOmelyan::Evaluate()
{
    NestedEvaluate(0, m_fTotalStepLength, TRUE, TRUE);
    FinishEvaluate();
}

/**
 * For last nest, it is exp(h/6 S) exp(h/2 T) exp(2/3 h S + 1/72 h^3 C) exp(h/2 T) exp(h/6 S)
 * For other nest, it is exp(h/6 S) NextNest exp(2/3 h S + 1/72 h^3 C) NextNest exp(h/6 S)
 */
void CIntegratorMultiLevelOmelyan::NestedEvaluate(INT iLevel, Real fNestedStepLength, UBOOL bFirst, UBOOL bLast)
{
    //It is h/2m for nested, so times 0.5
    const UINT uiStepAll = (0 == iLevel) ? m_uiStepCount : m_uiNestedStep[iLevel - 1];
    fNestedStepLength = fNestedStepLength / uiStepAll;
    const Real fHalfEstep = F(0.5) * fNestedStepLength;
    appDetailed(_T("  Omelyan level %d step 0, -------- length : %f \n"), iLevel, fNestedStepLength);

    UpdateP(m_f2Lambda * fHalfEstep, iLevel,
        bFirst ? ESP_StartTrajectory : ESP_InTrajectory,
        FALSE,
        TRUE);

    if (m_bDebugForce)
    {
        appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, CalcForce());
    }

    for (UINT uiStep = 1; uiStep < uiStepAll + 1; ++uiStep)
    {
        // middle step, exp(h/2 T) or NextNest
        if (iLevel != m_uiNestedStep.Num())
        {
            //NextNest
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(iLevel + 1, fHalfEstep, bFirst, FALSE);
            }
            else
            {
                NestedEvaluate(iLevel + 1, fHalfEstep, bFirst, FALSE);
            }
        }
        else
        {
            //exp(h/2 T)
            //appDetailed("  Omelyan nested level %d step %d U1\n", iLevel, uiStep);
            UpdateU(fHalfEstep);
        }

        appDetailed("  Omelyan nested level %d step %d P Middle\n", iLevel, uiStep);
        UpdateP(fNestedStepLength * (F(1.0) - m_f2Lambda), iLevel, 
            ESP_InTrajectory, FALSE, TRUE);

        if (m_bDebugForce)
        {
            appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, CalcForce());
        }

        if (iLevel != m_uiNestedStep.Num())
        {
            //NextNest
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(iLevel + 1, fHalfEstep, FALSE, uiStep == uiStepAll);
            }
            else
            {
                NestedEvaluate(iLevel + 1, fHalfEstep, FALSE, uiStep == uiStepAll);
            }
        }
        else
        {
            //exp(h/2 T)
            //appDetailed("  Omelyan nested level %d step %d U2\n", iLevel, uiStep);
            UpdateU(fHalfEstep);
        }

        if (uiStep < uiStepAll)
        {
            appDetailed("  Omelyan nested level %d step %d P\n", iLevel, uiStep);
            UpdateP(fNestedStepLength * m_f2Lambda, iLevel, ESP_InTrajectory, FALSE, TRUE);
        }
        else
        {
            appDetailed("  Omelyan nested level %d last step %d\n", iLevel, uiStep);
            UpdateP(m_f2Lambda * fHalfEstep, iLevel, bLast ? ESP_EndTrajectory : ESP_InTrajectory, bLast, TRUE);
        }
        if (m_bDebugForce)
        {
            appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, CalcForce());
        }
    }
}

CCString CIntegratorMultiLevelOmelyan::GetInfos(const CCString& sTab) const
{
    CCString sRet = sTab + _T("Name : Nested Omelyan\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + sTab + _T("##Omelyan2Lambda = 2 x lambda\n");
    sRet = sRet + sTab + _T("Omelyan2Lambda : ") + appToString(m_f2Lambda) + _T("\n");
    sRet = sRet + GetNestedInfo(sTab);
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================