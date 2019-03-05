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

#define OneOver6 (F(0.16666666666666666666666666666667))
#define OneOver24 (F(0.04166666666666666666666666666667))

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
    Real f1Over2EStep = m_fEStep * F(0.5);
    Real f1Over6Estep = m_fEStep * OneOver6;
    Real f1Over3Estep = f1Over6Estep * F(2.0);
    Real f2Over3Estep = f1Over3Estep * F(2.0);
    Real f1Over24EstepSq = m_fEStep * m_fEStep * OneOver24;

    appDetailed("  Force Gradient sub step 0\n");
    UINT uiForceStep = 0;
    UpdateP(f1Over6Estep, FALSE, uiForceStep);
    ++uiForceStep;

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        UpdateU(f1Over2EStep);

        // middle step
        m_pForceField->Zero();
        checkCudaErrors(cudaDeviceSynchronize());
        for (INT i = 0; i < m_lstActions.Num(); ++i)
        {
            //this is accumulate
            m_lstActions[i]->CalculateForceOnGauge(uiForceStep, m_pGaugeField, m_pForceField, NULL);
            checkCudaErrors(cudaDeviceSynchronize());
        }
        uiForceStep++;

        m_pGaugeField->CopyTo(m_pUPrime);
        m_pForceField->ExpMult(f1Over24EstepSq, m_pGaugeField);

        UpdateP(f2Over3Estep, FALSE, uiForceStep);
        ++uiForceStep;

        //restore U
        m_pUPrime->CopyTo(m_pGaugeField);
        UpdateU(f1Over2EStep);

        if (uiStep < m_uiStepCount)
        {
            appDetailed("  Force Gradient sub step %d\n", uiStep);
            UpdateP(f1Over3Estep, FALSE, uiForceStep);
            ++uiForceStep;
        }
        else
        {
            appDetailed("  Force Gradient last step %d\n", uiStep);
            UpdateP(f1Over6Estep, TRUE, uiForceStep);
            ++uiForceStep;
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