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

void CMeasurePlaqutteEnergy::OnConfigurationAccepted(const CFieldGauge* pAcceptGauge, const CFieldGauge* pCorrespondingStaple)
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE plaqutteEneregy = 0.0;
#else
    Real plaqutteEneregy = F(0.0);
#endif
    if (NULL == pCorrespondingStaple || !CCommonData::m_bStoreStaple)
    {
#if !_CLG_DOUBLEFLOAT
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergy(1.0 / _HC_SUN);
#else
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergy(F(1.0) / _HC_SUN);
#endif
    }
    else
    {
#if !_CLG_DOUBLEFLOAT
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergyUsingStable(1.0 / _HC_SUN, pCorrespondingStaple);
#else
        plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergyUsingStable(F(1.0) / _HC_SUN, pCorrespondingStaple);
#endif
    }
    plaqutteEneregy = plaqutteEneregy / _HC_PlaqutteCount;
#if !_CLG_DOUBLEFLOAT
    const Real plaqEnergy = F(1.0) - static_cast<Real>(plaqutteEneregy);
#else
    const Real plaqEnergy = F(1.0) - plaqutteEneregy;
#endif
    m_lstData.AddItem(plaqEnergy);
    m_fLastRealResult = plaqEnergy;
    appParanoiac(_T(" === Plaqutte Energy Measured === energy = %f\n"), plaqEnergy);
}

void CMeasurePlaqutteEnergy::Average(UINT )
{
    Real fAdd = F(0.0);
    for (INT i = 0; i < m_lstData.Num(); ++i)
    {
        fAdd += m_lstData[i];
    }
    m_fLastRealResult = fAdd / m_lstData.Num();
    appParanoiac(_T(" === Plaqutte Energy Averaged (%d measures) === energy = %f\n"), m_lstData.Num(), m_fLastRealResult);
}

void CMeasurePlaqutteEnergy::Report()
{
    appGeneral(_T(" === Plaqutte Energy Averaged === energy = %f\n\n"), m_fLastRealResult);
}

void CMeasurePlaqutteEnergy::Reset()
{
    m_lstData.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================