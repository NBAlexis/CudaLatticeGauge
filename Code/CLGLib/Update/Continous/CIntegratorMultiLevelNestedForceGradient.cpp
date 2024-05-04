//=============================================================================
// FILENAME : CIntegratorMultiLevelNestedForceGradient.cpp
// 
// DESCRIPTION:
// This is the Approximate force gradient integrator for HMC
//
// REVISION:
//  [08/17/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorMultiLevelNestedForceGradient)

void CIntegratorMultiLevelNestedForceGradient::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CMultiLevelNestedIntegrator::Initial(pOwner, pLattice, params);
    
    CreateBackupFields();
}

void CIntegratorMultiLevelNestedForceGradient::Evaluate()
{
    NestedEvaluate(0, m_fTotalStepLength, TRUE, TRUE);
    FinishEvaluate();
}

/**
 * For last nest, it is exp(h/6 S) exp(h/2 T) exp(2/3 h S + 1/72 h^3 C) exp(h/2 T) exp(h/6 S)
 * For other nest, it is exp(h/6 S) NextNest exp(2/3 h S + 1/72 h^3 C) NextNest exp(h/6 S)
 */
void CIntegratorMultiLevelNestedForceGradient::NestedEvaluate(INT iLevel, Real fNestedStep, UBOOL bFirst, UBOOL bLast)
{
    //It is h/2m for nested, so times 0.5
    const UINT uiStepAll = (0 == iLevel) ? m_uiStepCount : m_uiNestedStep[iLevel - 1];
    const Real fNestedStepLength = fNestedStep / uiStepAll;
    const Real f1Over2EStep = fNestedStepLength * F(0.5);
    const Real f1Over6Estep = fNestedStepLength * OneOver6;
    const Real f1Over3Estep = f1Over6Estep * F(2.0);
    const Real f2Over3Estep = f1Over3Estep * F(2.0);
    const Real f1Over24EstepSq = fNestedStepLength * fNestedStepLength * OneOver24;

    appDetailed(_T("  Force Gradient nested level %d step 0, length=%f\n"), iLevel, fNestedStepLength);

    //exp(h/6 S), only iLevel == 0 will refresh force
    UpdateP(f1Over6Estep, iLevel, 
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
            appDetailed("  Force Gradient nested level %d step %d\n", iLevel, uiStep);
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(iLevel + 1, f1Over2EStep, bFirst, FALSE);
            }
            else
            {
                NestedEvaluate(iLevel + 1, f1Over2EStep, bFirst, FALSE);
            }
        }
        else
        {
            //exp(h/2 T)
            appDetailed("  Force Gradient nested level %d step %d\n", iLevel, uiStep);
            UpdateU(f1Over2EStep);
        }

        //exp(2/3 h S + 1/72 h^3 C)
        //Use U to calculate Force1
        //Use Force1 to update U'
        //Use U' to update P
        //Recover U' to U

        //If not update momentum, the step length is unused. So set to 1.0
        UpdateP(F(1.0), iLevel,
            ESP_InTrajectory,
            FALSE,
            FALSE);

        PreserveFields();
        AddForceToFieldDirectly(f1Over24EstepSq);

        UpdateP(f2Over3Estep, iLevel, ESP_InTrajectory, FALSE, TRUE);
        if (m_bDebugForce)
        {
            appGeneral(_T(" ------ Force (%d) = %f\n"), iLevel, CalcForce());
        }
        //restore U
        RecoverFields();

        if (iLevel != m_uiNestedStep.Num())
        {
            //NextNest
            appDetailed("  Force Gradient nested level %d step %d\n", iLevel, uiStep);
            if (m_bInnerLeapFrog)
            {
                NestedEvaluateLeapfrog(iLevel + 1, f1Over2EStep, FALSE, uiStep == uiStepAll);
            }
            else
            {
                NestedEvaluate(iLevel + 1, f1Over2EStep, FALSE, uiStep == uiStepAll);
            }
        }
        else
        {
            //exp(h/2 T)
            appDetailed("  Force Gradient nested level %d step %d\n", iLevel, uiStep);
            UpdateU(f1Over2EStep);
        }

        if (uiStep < uiStepAll)
        {
            appDetailed("  Force Gradient nested level %d step %d\n", iLevel, uiStep);
            UpdateP(f1Over3Estep, iLevel, ESP_InTrajectory, FALSE, TRUE);
        }
        else
        {
            appDetailed("  Force Gradient nested last step %d\n", uiStep);
            UpdateP(f1Over6Estep, iLevel, bLast ? ESP_EndTrajectory : ESP_InTrajectory, bLast, TRUE);
        }
        if (m_bDebugForce)
        {
            appGeneral(_T(" ------ Force (%d) = %f\n"), iLevel, CalcForce());
        }
    }
}

CCString CIntegratorMultiLevelNestedForceGradient::GetInfos(const CCString& sTab) const
{
    CCString sRet = sTab + _T("Name : Nested Force Gradient\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + GetNestedInfo(sTab);
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================