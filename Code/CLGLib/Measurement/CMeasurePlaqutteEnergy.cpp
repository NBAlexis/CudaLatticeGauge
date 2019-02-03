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

void CMeasurePlaqutteEnergy::OnConfigurationAccepted()
{
    Real plaqutteEneregy = m_pLatticeData->m_pGaugeField->CalculatePlaqutteEnergy(F(1.0) / _HC_SUN);
    plaqutteEneregy = plaqutteEneregy / _HC_PlaqutteCount;
    Real plaqEnergy = F(1.0) - plaqutteEneregy;
    m_lstData.AddItem(plaqEnergy);
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