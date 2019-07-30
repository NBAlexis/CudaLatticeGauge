//=============================================================================
// FILENAME : TestZ3Symmetry.cpp
// 
// DESCRIPTION:
//
//     Test the Z3 Symmetry
//
// REVISION:
//  [06/23/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestZ3Symmetry(CParameters&)
{
    //appGetLattice()->m_pUpdator->Update(10, FALSE);
    //appGetLattice()->m_pUpdator->Update(20, TRUE);
    //appGetLattice()->m_pMeasurements->Report();
    CFieldFermionWilsonSquareSU3* pCpy = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2)->GetCopy());

    SFermionSource sourceData;
    sourceData.m_eSourceType = EFS_Point;
    sourceData.m_sSourcePoint.x = 0;
    sourceData.m_sSourcePoint.y = 0;
    sourceData.m_sSourcePoint.z = 0;
    sourceData.m_sSourcePoint.w = 0;

    sourceData.m_byColorIndex = 1;
    sourceData.m_bySpinIndex = 1;
    pCpy->InitialAsSource(sourceData);

    pCpy->InverseD(appGetLattice()->m_pGaugeField);

    return 0;
}


__REGIST_TEST(TestZ3Symmetry, Misc, TestZ3Symmetry);



//=============================================================================
// END OF FILE
//=============================================================================
