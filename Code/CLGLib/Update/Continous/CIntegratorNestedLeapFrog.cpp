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

__CLGIMPLEMENT_CLASS(CIntegratorNestedLeapFrog)

void CIntegratorNestedLeapFrog::Evaluate()
{
    if (m_lstActions.Num() < 2)
    {
        appCrucial(_T("Nested Updator only work with actions more than 2!"));
        _FAIL_EXIT;
    }
    Real fHalfPstep = F(0.5) * m_fEStep;
    UpdatePF(fHalfPstep, ESP_StartTrajectory);
    appDetailed("  leap frog sub step 0\n");

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        if (uiStep < m_uiStepCount)
        {
            NestedEvaluate(FALSE);
            UpdatePF(m_fEStep, ESP_InTrajectory);
            appDetailed("  leap frog sub step %d\n", uiStep);
        }
        else
        {
            NestedEvaluate(TRUE);
            UpdatePF(fHalfPstep, ESP_EndTrajectory);
            appDetailed("  leap frog last step %d\n", uiStep);
        }
    }

    FinishEvaluate();
}

void CIntegratorNestedLeapFrog::NestedEvaluate(UBOOL bLast)
{
    Real fHalfPstep = F(0.5) * m_fNestedStepLength;
    UpdatePG(fHalfPstep, FALSE);
    appDetailed("  leap frog nested sub step 0\n");

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(m_fNestedStepLength);

        if (uiStep < m_uiNestedStep)
        {
            UpdatePG(m_fNestedStepLength, FALSE);
            appDetailed("  leap frog nested sub step %d\n", uiStep);
        }
        else
        {
            UpdatePG(fHalfPstep, bLast);
            appDetailed("  leap frog nested last step %d\n", uiStep);
        }
    }
}

CCString CIntegratorNestedLeapFrog::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Nested leap frog\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appFloatToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appIntToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appFloatToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + GetNestedInfo(sTab);
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================