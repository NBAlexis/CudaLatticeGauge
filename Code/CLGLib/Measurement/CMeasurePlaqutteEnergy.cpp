//=============================================================================
// FILENAME : CMeasurePlaqutteEnergy.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePlaqutteEnergy)

void CMeasurePlaqutteEnergy::OnConfigurationAcceptedSingleField(const CFieldGauge* pAcceptGauge, const CFieldGauge* pCorrespondingStaple)
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE plaqutteEneregy = 0.0;
#else
    Real plaqutteEneregy = F(0.0);
#endif
    if (NULL == pCorrespondingStaple || !CCommonData::m_bStoreStaple)
    {
#if !_CLG_DOUBLEFLOAT
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergy(1.0 / GetDefaultMatrixN());
#else
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergy(F(1.0) / GetDefaultMatrixN());
#endif
    }
    else
    {
#if !_CLG_DOUBLEFLOAT
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergyUsingStable(1.0 / GetDefaultMatrixN(), pCorrespondingStaple);
#else
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergyUsingStable(F(1.0) / GetDefaultMatrixN(), pCorrespondingStaple);
#endif
    }
    plaqutteEneregy = plaqutteEneregy / _HC_PlaqutteCount;
#if !_CLG_DOUBLEFLOAT
    const Real plaqEnergy = F(1.0) - static_cast<Real>(plaqutteEneregy);
#else
    const Real plaqEnergy = F(1.0) - plaqutteEneregy;
#endif
    UpdateRealResult(plaqEnergy);
    appParanoiac(_T(" === Plaqutte Energy Measured === energy = %f\n"), plaqEnergy);
}

void CMeasurePlaqutteEnergy::Report()
{
    Average();
    appGeneral(_T(" === Plaqutte Energy Averaged === energy = %f\n\n"), GetAverageRealRes());
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================