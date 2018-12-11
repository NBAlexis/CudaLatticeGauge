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

__global__ void _kernelInitialActionPlaquette(CLatticeData * pOwner, 
    CDeviceLattice *& pDevice)
{
    pDevice = pOwner->GetDeviceInstance();
}

CActionGaugePlaquette::CActionGaugePlaquette()
    : CAction()
{
    m_pOwner = CLatticeData::GetInstance();

    _kernelInitialActionPlaquette << <1, 1 >> > (m_pOwner, m_pDeviceLattce);
}

void CActionGaugePlaquette::CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const
{
    preparethread;
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_cMinusBetaOverN);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================