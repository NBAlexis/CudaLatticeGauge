//=============================================================================
// FILENAME : CIntegratorOmelyan.cpp
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [02/12/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorOmelyan)

void CIntegratorOmelyan::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);
    m_f2Lambda = OmelyanLambda2;
    if (!params.FetchValueReal(_T("Omelyan2Lambda"), m_f2Lambda))
    {
        m_f2Lambda = OmelyanLambda2;
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================