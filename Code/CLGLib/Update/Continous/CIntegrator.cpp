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
    m_pGaugeField->CachePlaqutteIndexes();

    m_pForceField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pForceField->m_pOwner = pLattice;
    m_pForceField->InitialField(EFIT_Zero);

    m_pMomentumField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pMomentumField->m_pOwner = pLattice;
    m_pMomentumField->InitialField(EFIT_Zero);

    m_pStapleField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pStapleField->m_pOwner = pLattice;
    m_pStapleField->InitialField(EFIT_Zero);
}

void CIntegrator::Prepare(UBOOL bLastAccepted)
{
    //we may not accept the evaluation, so we need to copy it first
    if (!bLastAccepted)
    {
        m_pLattice->m_pGaugeField->CopyTo(m_pGaugeField);
        m_bStapleCached = FALSE;
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //generate a random momentum field to start
    m_pMomentumField->MakeRandomGenerator();
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::Accept()
{
    m_pGaugeField->CopyTo(m_pLattice->m_pGaugeField);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::UpdateU(Real fStep)
{
    //U(k) = exp (i e P) U(k-1)
    m_pMomentumField->ExpMult(_make_cuComplex(F(0.0), fStep), m_pGaugeField);
    checkCudaErrors(cudaDeviceSynchronize());
    //m_pGaugeField->DebugPrintMe();
}

void CIntegrator::UpdateP(Real fStep, UBOOL bCacheStaple)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, (0 == i && bCacheStaple) ? m_pStapleField : NULL);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //P = P + e F
    m_bStapleCached = bCacheStaple;
    m_pMomentumField->Axpy(_make_cuComplex(F(0.0), fStep), m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());
}

Real CIntegrator::GetEnergy() const
{
#if _CLG_DEBUG
    CCString sLog = _T("");
#endif

    //get P energy
    Real retv = m_pMomentumField->CalculateKinematicEnergy();

#if _CLG_DEBUG
    sLog.Format(_T("Kin:%f"), retv);
#endif

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
#if _CLG_DEBUG
        Real fActionEnergy = m_bStapleCached ? m_lstActions[i]->Energy(m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(m_pGaugeField);
        CCString sThisActionInfo = _T("");
        sThisActionInfo.Format(_T(" Action%d:%f, "), i + 1, fActionEnergy);
        sLog += sThisActionInfo;
        retv += fActionEnergy;
#else
        retv += m_bStapleCached ? m_lstActions[i]->Energy(m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(m_pGaugeField);
#endif
    }
#if _CLG_DEBUG
    appDetailed(_T("H:%s \n"), sLog.c_str());
#endif
    return retv;
}

__CLGIMPLEMENT_CLASS(CIntegratorLeapFrog)

__CLGIMPLEMENT_CLASS(CIntegratorOmelyan)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================