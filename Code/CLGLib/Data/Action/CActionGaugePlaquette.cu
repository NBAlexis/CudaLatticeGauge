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

CActionGaugePlaquette::CActionGaugePlaquette()
    : CAction()
{
}

void CActionGaugePlaquette::CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const
{
    pGauge->CalculateForceAndStaple(pGauge, pStaple, m_cMinusBetaOverN);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================