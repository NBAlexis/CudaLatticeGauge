//=============================================================================
// FILENAME : CMeasureChiralCondensate.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [06/13/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChiralCondensate)

CMeasureChiralCondensate::~CMeasureChiralCondensate()
{

}

void CMeasureChiralCondensate::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    Reset();

    INT iValue = 2;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("FieldCount"), iValue);
    m_uiFieldCount = static_cast<UINT>(iValue);
}

void CMeasureChiralCondensate::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    CLGComplex res = _make_cuComplex(F(0.0), F(0.0));

    CFieldFermion* pF1 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    CFieldFermion* pF2 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    UINT uiVolume = appGetLattice()->m_pIndexCache->m_uiSiteNumber[m_byFieldId];

    for (UINT i = 0; i < m_uiFieldCount; ++i)
    {
        pF1->InitialField(EFIT_RandomZ4);
        pF1->FixBoundary();
        //pF1->DebugPrintMe();
        //pGauge->DebugPrintMe();
        pF1->CopyTo(pF2);
        pF1->InverseD(pGauge);
        CLGComplex thisRes = pF2->Dot(pF1);
        res.x = res.x + thisRes.x / uiVolume;
        res.y = res.y + thisRes.y / uiVolume;
    }

    pF1->Return();
    pF2->Return();
    res.x = res.x / m_uiFieldCount;
    res.y = res.y / m_uiFieldCount;
    appDetailed(_T("\nChiral Condensate = %2.12f + %2.12f\n"), res.x, res.y);
    ++m_uiConfigurationCount;
    m_lstCondensate.AddItem(res);
}

void CMeasureChiralCondensate::Average(UINT )
{
    //nothing to do
}

void CMeasureChiralCondensate::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstCondensate.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Chiral Condensate (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- each configuration ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        tmpChargeSum.x += m_lstCondensate[i].x;
        tmpChargeSum.y += m_lstCondensate[i].y;
        LogGeneralComplex(m_lstCondensate[i]);
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"), 
        tmpChargeSum.x / m_uiConfigurationCount,
        tmpChargeSum.y / m_uiConfigurationCount);

    appGeneral(_T("==========================================================================\n\n"));
    appSetLogDate(TRUE);
}

void CMeasureChiralCondensate::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstCondensate.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================