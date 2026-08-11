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
#include "CMeasurePlaqutteEnergy.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePlaqutteEnergy)

void CMeasurePlaqutteEnergy::OnConfigurationAcceptedSingleField(const CFieldGauge* pAcceptGauge, const CFieldGauge* pCorrespondingStaple)
{
    DOUBLE plaqutteEneregy = 0.0;
    plaqutteEneregy = pAcceptGauge->CalculatePlaqutteEnergy(1.0 / pAcceptGauge->MatrixN());
    //P4-3.1: CalculatePlaqutteEnergy returns the LOCAL partial sum on multi-GPU;
    //reduce to the global value before normalizing / storing.
    GlobalSumReal(plaqutteEneregy);
    //plaqutteEneregy = (*pAcceptGauge)::
    //P4-3.2: after the global reduction the energy is the GLOBAL sum, so the
    //normalization count must also be the GLOBAL plaquette count (the local
    //_HC_PlaqutteCount is half of it on -n2 -> the result would be x2 off,
    //verified: plaq -0.999 on -n2 vs 3.07e-04 on -n1 before this fix).
    plaqutteEneregy = plaqutteEneregy / GlobalPlaqutteCount();
    const Real plaqEnergy = F(1.0) - static_cast<Real>(plaqutteEneregy);

    DOUBLE u0Raw = pAcceptGauge->CalculatePlaqutteEnergyOriginal(1.0 / pAcceptGauge->MatrixN());
    GlobalSumReal(u0Raw);
    appGeneral(_T("plaq = %f\n"), 1.0 - u0Raw / GlobalPlaqutteCount());
    DOUBLE u0 = pow(1.0 - u0Raw / GlobalPlaqutteCount(), 0.25);

    DOUBLE v0 = 0.0;
    if (m_bV0)
    {
        CGaugeSmearing* smearing = appGetGaugeSmearing(pAcceptGauge->m_byFieldId);
        if (NULL != smearing && smearing->CalledWhenUpdate())
        {
            smearing->GaugeSmearingC(pAcceptGauge);
            v0 = smearing->GetEffectiveGauge()->CalculatePlaqutteEnergyOriginal(1.0 / pAcceptGauge->MatrixN());
            GlobalSumReal(v0);
            v0 = pow(1.0 - v0 / GlobalPlaqutteCount(), 0.25);
            m_lstV0.AddItem(v0);
        }
    }

    UpdateRealResult(plaqEnergy);
    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("Plaquette"), plaqEnergy);
        m_pOwner->AddOneConfigurationResult(this, _T("U0"), u0);
        if (m_bV0)
        {
            m_pOwner->AddOneConfigurationResult(this, _T("V0"), v0);
        }
    }

    if (m_bShowResult)
    {
        if (m_bV0)
        {
            appGeneral(_T(" === Average gauge action = %f, u0 = %f, v0 = %f\n"), plaqEnergy, u0, v0);
        }
        else
        {
            appGeneral(_T(" === Average gauge action = %f, u0 = %f\n"), plaqEnergy, u0);
        }
    }
}

void CMeasurePlaqutteEnergy::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    INT iValue = 0;
    if (param.FetchValueINT(_T("v0"), iValue))
    {
        m_bV0 = (0 != iValue);
    }
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