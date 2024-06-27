//=============================================================================
// FILENAME : CMeasureBosonCond.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [28/06/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureBosonCond)


void CMeasureBosonCond::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple)
{
    if (m_lstBosonFieldIds.Num() < 1)
    {
        appCrucial(_T("CMeasureBosonCond, BosonFields not properly set\n"));
        return;
    }
    INT ibosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, pAcceptBoson, m_lstBosonFieldIds[0]);
    if (ibosonidx < 0 || ibosonidx >= bosonNum)
    {
        appCrucial(_T("CMeasureBosonCond, BosonFields not properly set\n"));
        return;
    }

    const CFieldBoson* bosonfield = pAcceptBoson[ibosonidx];
#if !_CLG_DOUBLEFLOAT
    UpdateComplexResult(_cToFloat(bosonfield->Dot(bosonfield)));
#else
    UpdateComplexResult(bosonfield->Dot(bosonfield));
#endif
}

void CMeasureBosonCond::Report()
{
    Average();

    appPushLogDate(FALSE);

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Boson condensation (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        LogGeneralComplex(CmpResAtI(i));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average = "));
    ReportAverageComplexRes();

    appGeneral(_T(" --------\n==========================================================================\n"));
    appPopLogDate();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================