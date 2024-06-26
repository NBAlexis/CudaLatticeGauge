//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggered.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [09/28/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE


__CLGIMPLEMENT_CLASS(CMeasureConnectedSusceptibilityKS)

CMeasureConnectedSusceptibilityKS::~CMeasureConnectedSusceptibilityKS()
{

}

void CMeasureConnectedSusceptibilityKS::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const CFieldGauge* const* pStapleField)
{
    m_pSourceZero = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId()));
    CFieldFermion* pSourceZeroCopy = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId()));
    SFermionSource sour;
    sour.m_byColorIndex = 0;
    sour.m_eSourceType = EFS_Point;
    sour.m_sSourcePoint = _HC_Center;//SSmallInt4(0, 0, 0, 0);
    //appGeneral(_T("point1 %d, point2 %d\n"), _hostGetSiteIndex(sour.m_sSourcePoint), _hostGetSiteIndex(CCommonData::m_sCenter));
    m_pSourceZero->InitialAsSource(sour);
    m_pSourceZero->FixBoundary();
    m_pSourceZero->CopyTo(pSourceZeroCopy);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
#if !_CLG_DOUBLEFLOAT
    const cuDoubleComplex color1 = pSourceZeroCopy->Dot(m_pSourceZero);
#else
    const CLGComplex color1 = pSourceZeroCopy->Dot(m_pSourceZero);
#endif
    sour.m_byColorIndex = 1;
    m_pSourceZero->InitialAsSource(sour);
    m_pSourceZero->CopyTo(pSourceZeroCopy);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
#if !_CLG_DOUBLEFLOAT
    const cuDoubleComplex color2 = pSourceZeroCopy->Dot(m_pSourceZero);
#else
    const CLGComplex color2 = pSourceZeroCopy->Dot(m_pSourceZero);
#endif

    sour.m_byColorIndex = 2;
    m_pSourceZero->InitialAsSource(sour);
    m_pSourceZero->CopyTo(pSourceZeroCopy);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
    m_pSourceZero->InverseD(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson);
#if !_CLG_DOUBLEFLOAT
    const cuDoubleComplex color3 = pSourceZeroCopy->Dot(m_pSourceZero);
#else
    const CLGComplex color3 = pSourceZeroCopy->Dot(m_pSourceZero);
#endif

#if !_CLG_DOUBLEFLOAT
    UpdateComplexResult(_cToFloat(cuCadd(cuCadd(color1, color2), color3)));
#else
    UpdateComplexResult(_cuCaddf(_cuCaddf(color1, color2), color3));
#endif
    pSourceZeroCopy->Return();

    if (m_bShowResult)
    {
        appGeneral(_T("Connected Chiral susp:"));
        LogGeneralComplex(GetLastCmpRes(), FALSE);
        appGeneral(_T("\n"));
    }
    m_pSourceZero->Return();
}

void CMeasureConnectedSusceptibilityKS::Report()
{
    appGeneral(_T(" =======================================================\n"));
    appGeneral(_T(" ================ Connected Susceptibility =============\n"));
    appGeneral(_T(" =======================================================\n\n"));
    CLGComplex average = _zeroc;
    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        LogGeneralComplex(CmpResAtI(i), FALSE);
    }
    appGeneral(_T("}\n\n"));

    average.x = average.x / m_uiConfigurationCount;
    average.y = average.y / m_uiConfigurationCount;

    appGeneral(_T(" ============average: "));
    LogGeneralComplex(average, FALSE);
    appGeneral(_T("========= \n"));
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================