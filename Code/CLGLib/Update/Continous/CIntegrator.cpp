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

/**
* Create the fields here
*/
void CIntegrator::Initial(const TArray<class CAction*>& actionList, UINT uiStepCount, UINT uiExpPrecision)
{

}

void CIntegrator::Prepare()
{
    //we may not accept the evaluation, so we need to copy it first
    m_pLattice->m_pGaugeField->CopyTo(m_pGaugeField);
    //generate a random momentum field to start
    m_pMomentumField->MakeRandomGenerator();
}

void CIntegrator::UpdateU(FLOAT fStep)
{
    m_pMomentumField->ExpMult(make_cuComplex(fStep, 0.0f), m_uiExpPrecision, m_pGaugeField);
}

void CIntegrator::UpdateP(FLOAT fStep)
{
    if (!m_bPCached)
    {
        // recalc force
        m_pForceField->Zero();

        for (INT i = 0; i < m_lstActions.Num(); ++i)
        {
            //this is accumulate
            m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, m_pStapleField);
        }

        m_bPCached = TRUE;
    }

    m_pMomentumField->Axpy(fStep, m_pForceField);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================