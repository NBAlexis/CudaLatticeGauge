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

CIntegratorNestedForceGradient::~CIntegratorNestedForceGradient()
{
    appSafeDelete(m_pUPrime);
}

void CIntegratorNestedForceGradient::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CNestedIntegrator::Initial(pOwner, pLattice, params);
    m_fNestedStepLength = F(0.5) * m_fNestedStepLength;
    m_pUPrime = dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField->GetCopy());
}

void CIntegratorNestedForceGradient::Evaluate()
{
    if (m_lstActions.Num() < 2)
    {
        appCrucial(_T("Nested Updator only work with actions more than 2!"));
        _FAIL_EXIT;
    }

    //Real f1Over2EStep = m_fEStep * F(0.5);
    Real f1Over6Estep = m_fEStep * OneOver6;
    Real f1Over3Estep = f1Over6Estep * F(2.0);
    Real f2Over3Estep = f1Over3Estep * F(2.0);
    Real f1Over24EstepSq = m_fEStep * m_fEStep * OneOver24;

    appDetailed("  Force Gradient sub step 0\n");
    UpdatePF(f1Over6Estep);

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        NestedEvaluate(FALSE);

        // middle step
        m_pForceField->Zero();
        checkCudaErrors(cudaDeviceSynchronize());
        for (INT i = 1; i < m_lstActions.Num(); ++i)
        {
            //this is accumulate
            m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL);
            checkCudaErrors(cudaDeviceSynchronize());
        }

        m_pGaugeField->CopyTo(m_pUPrime);
        m_pForceField->ExpMult(f1Over24EstepSq, m_pGaugeField);

        UpdatePF(f2Over3Estep);

        //restore U
        m_pUPrime->CopyTo(m_pGaugeField);
        

        if (uiStep < m_uiStepCount)
        {
            NestedEvaluate(FALSE);
            appDetailed("  Force Gradient sub step %d\n", uiStep);
            UpdatePF(f1Over3Estep);
        }
        else
        {
            NestedEvaluate(TRUE);
            appDetailed("  Force Gradient last step %d\n", uiStep);
            UpdatePF(f1Over6Estep);
        }
    }

    FinishEvaluate();
}

void CIntegratorNestedForceGradient::NestedEvaluate(UBOOL bLast)
{
    Real f1Over2EStep = m_fNestedStepLength * F(0.5);
    Real f1Over6Estep = m_fNestedStepLength * OneOver6;
    Real f1Over3Estep = f1Over6Estep * F(2.0);
    Real f2Over3Estep = f1Over3Estep * F(2.0);
    Real f1Over24EstepSq = m_fNestedStepLength * m_fNestedStepLength * OneOver24;

    appDetailed("  Force Gradient nested sub step 0\n");
    UpdatePG(f1Over6Estep, FALSE);

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(f1Over2EStep);

        // middle step
        m_pForceField->Zero();
        checkCudaErrors(cudaDeviceSynchronize());
        m_lstActions[0]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL);

        m_pGaugeField->CopyTo(m_pUPrime);
        m_pForceField->ExpMult(f1Over24EstepSq, m_pGaugeField);

        UpdatePG(f2Over3Estep, FALSE);

        //restore U
        m_pUPrime->CopyTo(m_pGaugeField);
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