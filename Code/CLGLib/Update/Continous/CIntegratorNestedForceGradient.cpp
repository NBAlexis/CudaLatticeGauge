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

__CLGIMPLEMENT_CLASS(CIntegratorNestedForceGradient)


void CIntegratorNestedForceGradient::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CNestedIntegrator::Initial(pOwner, pLattice, params);

    CreateBackupFields();
    m_fNestedStepLength = F(0.5) * m_fNestedStepLength;
}

void CIntegratorNestedForceGradient::Evaluate()
{
    if (m_lstActions.Num() < 2)
    {
        appCrucial(_T("Nested Updator only work with actions more than 2!"));
        _FAIL_EXIT;
    }

    //Real f1Over2EStep = m_fEStep * F(0.5);
    const Real f1Over6Estep = m_fEStep * OneOver6;
    const Real f1Over3Estep = f1Over6Estep * F(2.0);
    const Real f2Over3Estep = f1Over3Estep * F(2.0);
    const Real f1Over24EstepSq = m_fEStep * m_fEStep * OneOver24;

    appDetailed("  Force Gradient sub step 0\n");
    UpdatePF(f1Over6Estep, ESP_StartTrajectory);

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

        // middle step
        ZeroForce();
        checkCudaErrors(cudaDeviceSynchronize());
        for (INT i = 0; i < m_lstActions.Num(); ++i)
        {
            //this is accumulate
            if (m_lstActions[i]->IsFermion())
            {
                m_lstActions[i]->CalculateForce(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pForceField.GetData(), m_pBosonForceFields.GetData(), NULL, ESP_InTrajectory);
            }

            checkCudaErrors(cudaDeviceSynchronize());
        }

        PreserveFields();
        AddForceToFieldDirectly(f1Over24EstepSq);

        UpdatePF(f2Over3Estep, ESP_InTrajectory);

        //restore U
        RecoverFields();
        

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
            
            appDetailed("  Force Gradient sub step %d\n", uiStep);
            UpdatePF(f1Over3Estep, ESP_InTrajectory);
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
            
            appDetailed("  Force Gradient last step %d\n", uiStep);
            UpdatePF(f1Over6Estep, ESP_EndTrajectory);
        }
    }

    FinishEvaluate();
}

void CIntegratorNestedForceGradient::NestedEvaluate(UBOOL bLast)
{
    const Real f1Over2EStep = m_fNestedStepLength * F(0.5);
    const Real f1Over6Estep = m_fNestedStepLength * OneOver6;
    const Real f1Over3Estep = f1Over6Estep * F(2.0);
    const Real f2Over3Estep = f1Over3Estep * F(2.0);
    const Real f1Over24EstepSq = m_fNestedStepLength * m_fNestedStepLength * OneOver24;

    appDetailed("  Force Gradient nested sub step 0\n");
    UpdatePG(f1Over6Estep, FALSE);

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(f1Over2EStep);

        // middle step
        ZeroForce();
        checkCudaErrors(cudaDeviceSynchronize());

        for (INT i = 0; i < m_lstActions.Num(); ++i)
        {
            if (!m_lstActions[i]->IsFermion())
            {
                m_lstActions[i]->CalculateForce(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pForceField.GetData(), m_pBosonForceFields.GetData(), NULL, ESP_Once);
            }
        }

        PreserveFields();
        AddForceToFieldDirectly(f1Over24EstepSq);

        UpdatePG(f2Over3Estep, FALSE);

        //restore U
        RecoverFields();
        UpdateU(f1Over2EStep);

        if (uiStep < m_uiNestedStep)
        {
            appDetailed("  Force Gradient nested sub step %d\n", uiStep);
            UpdatePG(f1Over3Estep, FALSE);
        }
        else
        {
            appDetailed("  Force Gradient nested last step %d\n", uiStep);
            UpdatePG(f1Over6Estep, bLast);
        }
    }
}

CCString CIntegratorNestedForceGradient::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Nested Force Gradient\n");
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