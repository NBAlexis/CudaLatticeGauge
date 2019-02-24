//=============================================================================
// FILENAME : TestGaugeSmearing.cpp
// 
// DESCRIPTION:
//
//     
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestGaugeSmearing(CParameters&)
{
    Real fBetaOverN = F(2.5);
    CFieldGauge* pGauge = appGetLattice()->m_pGaugeField;
    CFieldFermionWilsonSquareSU3* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));

    CFieldGauge* pForce = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    CFieldGauge* pStaple = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());

    //calculate staple for smearing
    pGauge->CalculateOnlyStaple(pStaple);

    //Calculate energy before smearing
    Real fGaugeEnergy = pGauge->CalculatePlaqutteEnergyUsingStable(fBetaOverN, pStaple);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    assert(NULL != pPooled);
    pFermion->CopyTo(pPooled);
    pPooled->InverseD(pGauge);
    Real fFermionEnergy = pPooled->Dot(pPooled).x;

    //Do the smearing
    appGetGaugeSmearing()->GaugeSmearing(pGauge, pStaple);

    Real fGaugeEnergy2 = pGauge->CalculatePlaqutteEnergy(fBetaOverN);
    pFermion->CopyTo(pPooled);
    pPooled->InverseD(pGauge);
    Real fFermionEnergy2 = pPooled->Dot(pPooled).x;

    appGeneral(_T("E: G:%f(%f) F:%f(%f), Total = %f(%f) diff=%f"), 
        fGaugeEnergy, 
        fGaugeEnergy2, 
        fFermionEnergy, 
        fFermionEnergy2, 
        fGaugeEnergy + fFermionEnergy,
        fGaugeEnergy2 + fFermionEnergy2,
        fGaugeEnergy + fFermionEnergy - (fGaugeEnergy2 + fFermionEnergy2));

    return 0;
}

//930 - 990 ms
__REGIST_TEST(TestGaugeSmearing, Misc, TestGaugeSmearingAPEProj);

__REGIST_TEST(TestGaugeSmearing, Misc, TestGaugeSmearingAPEStout);


//=============================================================================
// END OF FILE
//=============================================================================
