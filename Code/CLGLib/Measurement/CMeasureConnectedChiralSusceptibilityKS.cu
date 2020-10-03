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

void CMeasureConnectedSusceptibilityKS::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasureConnectedSusceptibilityKS::OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField)
{
    m_pSourceZero = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    SFermionSource sour;
    sour.m_byColorIndex = 0;
    sour.m_eSourceType = EFS_Point;
    sour.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
    m_pSourceZero->InitialAsSource(sour);
    m_pSourceZero->InverseD(pGaugeField);
    const CLGComplex color1 = m_pSourceZero->Dot(m_pSourceZero);

    //sour.m_byColorIndex = 1;
    //m_pSourceZero->InitialAsSource(sour);
    //m_pSourceZero->InverseD(pGaugeField);
    //const CLGComplex color2 = m_pSourceZero->Dot(m_pSourceZero);

    //sour.m_byColorIndex = 2;
    //m_pSourceZero->InitialAsSource(sour);
    //m_pSourceZero->InverseD(pGaugeField);
    //const CLGComplex color3 = m_pSourceZero->Dot(m_pSourceZero);
    //m_lstResults.AddItem(_cuCaddf(_cuCaddf(color1, color2), color3));
    m_lstResults.AddItem(color1);

    if (m_bShowResult)
    {
        appGeneral(_T("Connected Chiral susp:"));
        LogGeneralComplex(m_lstResults[m_lstResults.Num() - 1], FALSE);
        appGeneral(_T("\n"));
    }
    m_pSourceZero->Return();
    ++m_uiConfigurationCount;
}

void CMeasureConnectedSusceptibilityKS::Average(UINT)
{

}

void CMeasureConnectedSusceptibilityKS::Report()
{
    appGeneral(_T(" =======================================================\n"));
    appGeneral(_T(" ================ Connected Susceptibility =============\n"));
    appGeneral(_T(" =======================================================\n\n"));
    CLGComplex average = _zeroc;
    appGeneral(_T("{"));
    for (INT i = 0; i < m_lstResults.Num(); ++i)
    {
        LogGeneralComplex(m_lstResults[m_lstResults.Num() - 1], FALSE);
    }
    appGeneral(_T("}\n\n"));

    average.x = average.x / m_lstResults.Num();
    average.y = average.y / m_lstResults.Num();

    appGeneral(_T(" ============average: "));
    LogGeneralComplex(average, FALSE);
    appGeneral(_T("========= \n"));
}

void CMeasureConnectedSusceptibilityKS::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstResults.RemoveAll();
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================