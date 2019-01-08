//=============================================================================
// FILENAME : CActionGaugePlaquette.cu
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquette)

CActionGaugePlaquette::CActionGaugePlaquette()
    : CAction()
{
}

void CActionGaugePlaquette::Initial(class CLatticeData* pOwner, const CParameters& param)
{
    m_pOwner = pOwner;
    Real fBeta = 0.1f;
    param.FetchValueReal(_T("Beta"), fBeta);
    if (NULL != pOwner->m_pGaugeField && EFT_GaugeSU3 == pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / 3;
    }
    m_cMinusBetaOverN = _make_cuComplex(-fBeta, 0);
}

void CActionGaugePlaquette::CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const
{
    pGauge->CalculateForceAndStaple(pGauge, pStaple, m_cMinusBetaOverN);
    checkCudaErrors(cudaDeviceSynchronize());
}

/**
* The implementation depends on the type of gauge field
*/
Real CActionGaugePlaquette::Energy(const class CFieldGauge* pGauge) const
{
    return pGauge->CalculatePlaqutteEnergy(m_cMinusBetaOverN);
    checkCudaErrors(cudaDeviceSynchronize());
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================