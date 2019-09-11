//=============================================================================
// FILENAME : CIntegrator.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CIntegrator::~CIntegrator()
{
    appSafeDelete(m_pGaugeField);
    appSafeDelete(m_pForceField);
    appSafeDelete(m_pMomentumField);
    appSafeDelete(m_pStapleField);
}

/**
* Create the fields here
*/
void CIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    m_pOwner = pOwner;
    m_pLattice = pLattice;
    m_lstActions = pLattice->m_pActionList;
    m_bStapleCached = FALSE;

    INT iStepCount = 50;
    params.FetchValueINT(_T("IntegratorStep"), iStepCount);
    Real fStepLength = F(1.0);
    params.FetchValueReal(_T("IntegratorStepLength"), fStepLength);
    m_uiStepCount = (UINT)iStepCount;
    m_fEStep = fStepLength / m_uiStepCount;

    m_pGaugeField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pGaugeField->m_pOwner = pLattice;
    m_pGaugeField->InitialField(EFIT_Zero);

    m_pForceField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pForceField->m_pOwner = pLattice;
    m_pForceField->InitialField(EFIT_Zero);

    m_pMomentumField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pMomentumField->m_pOwner = pLattice;
    m_pMomentumField->InitialField(EFIT_Zero);

    if (CCommonData::m_bStoreStaple)
    {
        m_pStapleField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
        m_pStapleField->m_pOwner = pLattice;
        m_pStapleField->InitialField(EFIT_Zero);
    }
    else
    {
        m_pStapleField = NULL;
    }
}

void CIntegrator::Prepare(UBOOL bLastAccepted, UINT uiStep)
{
    //we may not accept the evaluation, so we need to copy it first
    if (!bLastAccepted || 0 == uiStep)
    {
        m_pLattice->m_pGaugeField->CopyTo(m_pGaugeField);
        m_bStapleCached = FALSE;
        checkCudaErrors(cudaDeviceSynchronize());
    }

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->PrepareForHMC(m_pGaugeField, uiStep);
    }

    //generate a random momentum field to start
    m_pMomentumField->MakeRandomGenerator();
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_pGaugeField->CopyTo(m_pLattice->m_pGaugeField);
    }
    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->OnFinishTrajectory(bAccepted);
    }
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::UpdateU(Real fStep) const
{
    //U(k) = exp (i e P) U(k-1)
    m_pMomentumField->ExpMult(fStep, m_pGaugeField);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::UpdateP(Real fStep, UBOOL bCacheStaple, ESolverPhase ePhase)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, (0 == i && bCacheStaple) ? m_pStapleField : NULL, ePhase);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //P = P + e F
    m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    m_pMomentumField->Axpy(fStep, m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::FinishEvaluate() const
{
    m_pGaugeField->ElementNormalize();
}

Real CIntegrator::GetEnergy(UBOOL bBeforeEvolution) const
{
    Real retv = m_pMomentumField->CalculateKinematicEnergy();

#if _CLG_DEBUG
    CCString sLog = _T("");
    sLog.Format(_T("kin:%f, "), retv);
#endif

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
#if _CLG_DEBUG
        Real fActionEnergy = m_bStapleCached ? m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField);
        CCString sThisActionInfo = _T("");
        sThisActionInfo.Format(_T(" Action%d:%f, "), i + 1, fActionEnergy);
        sLog += sThisActionInfo;
        retv += fActionEnergy;
#else
        retv += m_bStapleCached ? m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField);
#endif
    }
#if _CLG_DEBUG
    appDetailed(_T("H (%s) = %s \n"), bBeforeEvolution ? "before" : "after" , sLog.c_str());
#endif
    return retv;
}

void CNestedIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);

    INT iNestedStepCount = 3;
    params.FetchValueINT(_T("NestedStep"), iNestedStepCount);
    if (iNestedStepCount < 1)
    {
        appCrucial(_T("NestedStep must >= 1, but set to be %d!\n"), iNestedStepCount)
    }
    m_uiNestedStep = static_cast<UINT>(iNestedStepCount);
    m_fNestedStepLength = m_fEStep / m_uiNestedStep;
}

CCString CNestedIntegrator::GetNestedInfo(const CCString & sTab) const
{
    return sTab + _T("Nested : ") + appIntToString(static_cast<INT>(m_uiNestedStep)) + _T("\n");
}

void CNestedIntegrator::UpdatePF(Real fStep, ESolverPhase ePhase)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    for (INT i = 1; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL, ePhase);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //P = P + e F
    m_pMomentumField->Axpy(fStep, m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CNestedIntegrator::UpdatePG(Real fStep, UBOOL bCacheStaple)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    m_lstActions[0]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, bCacheStaple ? m_pStapleField : NULL, ESP_Once);
    checkCudaErrors(cudaDeviceSynchronize());

    //P = P + e F
    m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    m_pMomentumField->Axpy(fStep, m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());
}

__CLGIMPLEMENT_CLASS(CIntegratorLeapFrog)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================