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

__CLGIMPLEMENT_CLASS(CIntegratorForceGradient)

CIntegratorForceGradient::~CIntegratorForceGradient()
{
    appSafeDelete(m_pUPrime);
}

void CIntegratorForceGradient::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);
    
    m_pUPrime = dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField->GetCopy());
}

void CIntegratorForceGradient::Evaluate()
{
    const Real f1Over2EStep = m_fEStep * F(0.5);
    const Real f1Over6Estep = m_fEStep * OneOver6;
    const Real f1Over3Estep = f1Over6Estep * F(2.0);
    const Real f2Over3Estep = f1Over3Estep * F(2.0);
    const Real f1Over24EstepSq = m_fEStep * m_fEStep * OneOver24;

    appDetailed("  Force Gradient sub step 0\n");
    UpdateP(f1Over6Estep, FALSE, ESP_StartTrajectory);

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        UpdateU(f1Over2EStep);

        // middle step
        m_pForceField->Zero();
        checkCudaErrors(cudaDeviceSynchronize());
        for (INT i = 0; i < m_lstActions.Num(); ++i)
        {
            //this is accumulate
            m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL, ESP_InTrajectory);
            checkCudaErrors(cudaDeviceSynchronize());
        }

        m_pGaugeField->CopyTo(m_pUPrime);
        m_pForceField->ExpMult(f1Over24EstepSq, m_pGaugeField);

        UpdateP(f2Over3Estep, FALSE, ESP_InTrajectory);

        //restore U
        m_pUPrime->CopyTo(m_pGaugeField);
        UpdateU(f1Over2EStep);

        if (uiStep < m_uiStepCount)
        {
            appDetailed("  Force Gradient sub step %d\n", uiStep);
            UpdateP(f1Over3Estep, FALSE, ESP_InTrajectory);
        }
        else
        {
            appDetailed("  Force Gradient last step %d\n", uiStep);
            UpdateP(f1Over6Estep, TRUE, ESP_EndTrajectory);
        }
    }

    FinishEvaluate();
}

CCString CIntegratorForceGradient::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Force Gradient\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appFloatToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appIntToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appFloatToString(m_fEStep * m_uiStepCount) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================