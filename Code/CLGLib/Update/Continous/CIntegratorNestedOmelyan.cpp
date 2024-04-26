//=============================================================================
// FILENAME : CIntegratorForceGradient.cpp
// 
// DESCRIPTION:
// This is the Approximate force gradient integrator for HMC
//
// REVISION:
//  [03/05/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorNestedOmelyan)

void CIntegratorNestedOmelyan::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CNestedIntegrator::Initial(pOwner, pLattice, params);
    m_f2Lambda = OmelyanLambda2;
    if (!params.FetchValueReal(_T("Omelyan2Lambda"), m_f2Lambda))
    {
        m_f2Lambda = OmelyanLambda2;
    }
    m_fNestedStepLength = F(0.5) * m_fNestedStepLength;
}

void CIntegratorNestedOmelyan::Evaluate()
{
    if (m_lstActions.Num() < 2)
    {
        appCrucial(_T("Nested Updator only work with actions more than 2!"));
        _FAIL_EXIT;
    }

    const Real fHalfEstep = F(0.5) * m_fEStep;
    appDetailed("  Omelyan sub step 0\n");
    UpdatePF(m_f2Lambda * fHalfEstep, ESP_StartTrajectory);

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        if (m_bInnerLeapFrog)
        {
            NestedEvaluateLeapfrog(FALSE);
        }
        else
        {
            NestedEvaluate(FALSE);
        }
        
        UpdatePF(m_fEStep * (F(1.0) - m_f2Lambda), ESP_InTrajectory);

        if (uiStep < m_uiStepCount)
        {
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(FALSE);
            }
            else
            {
                NestedEvaluate(FALSE);
            }
            appDetailed("  Omelyan sub step %d\n", uiStep);
            UpdatePF(m_fEStep * m_f2Lambda, ESP_InTrajectory);
        }
        else
        {
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(TRUE);
            }
            else
            {
                NestedEvaluate(TRUE);
            }
            appDetailed("  Omelyan last step %d\n", uiStep);
            UpdatePF(m_f2Lambda * fHalfEstep, ESP_EndTrajectory);
        }
    }

    FinishEvaluate();
}

void CIntegratorNestedOmelyan::NestedEvaluate(UBOOL bLast)
{
    const Real fHalfEstep = F(0.5) * m_fNestedStepLength;
    appDetailed("  Omelyan nested sub step 0\n");
    UpdatePG(m_f2Lambda * fHalfEstep, FALSE);

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(fHalfEstep);
        UpdatePG(m_fNestedStepLength * (F(1.0) - m_f2Lambda), FALSE);
        UpdateU(fHalfEstep);

        if (uiStep < m_uiNestedStep)
        {
            appDetailed("  Omelyan nested sub step %d\n", uiStep);
            UpdatePG(m_fNestedStepLength * m_f2Lambda, FALSE);
        }
        else
        {
            appDetailed("  Omelyan nested last step %d\n", uiStep);
            UpdatePG(m_f2Lambda * fHalfEstep, bLast);
        }
    }
}

CCString CIntegratorNestedOmelyan::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Nested Omelyan\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appAnyToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appAnyToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appAnyToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + sTab + _T("##Omelyan2Lambda = 2 x lambda\n");
    sRet = sRet + sTab + _T("Omelyan2Lambda : ") + appAnyToString(m_f2Lambda) + _T("\n");
    sRet = sRet + GetNestedInfo(sTab);
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================