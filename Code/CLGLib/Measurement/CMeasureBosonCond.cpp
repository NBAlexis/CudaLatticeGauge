//=============================================================================
// FILENAME : CMeasureBosonCond.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [28/06/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CMeasureBosonCond.h"

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
    const CLGComplex condsq = _cToRealC(bosonfield->Dot(bosonfield));
    UpdateRealResult(_cuCabsf(condsq) / _HC_Volume / bosonfield->VectorN());

    m_lstElement.AddItem(bosonfield->Sum());
}

void CMeasureBosonCond::Reset()
{
    CMeasure::Reset();
    m_lstElement.RemoveAll();
    m_lstAverageElement.RemoveAll();
}

void CMeasureBosonCond::Average()
{
    CMeasure::Average();
    assert(m_lstElement.Num() == static_cast<INT>(m_uiConfigurationCount));
    assert(m_uiConfigurationCount > 0);

    m_lstAverageElement.RemoveAll();
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        if (0 == i)
        {
            for (INT j = 0; j < m_lstElement[0].Num(); ++j)
            {
                m_lstAverageElement.AddItem(m_lstElement[0][j]);
            }
        }
        else
        {
            for (INT j = 0; j < m_lstElement[i].Num(); ++j)
            {
                m_lstAverageElement[j] += m_lstElement[i][j];
            }
        }
    }
    for (INT j = 0; j < m_lstAverageElement.Num(); ++j)
    {
        m_lstAverageElement[j] = m_lstAverageElement[j] / m_uiConfigurationCount;
    }
}

void CMeasureBosonCond::Report()
{
    Average();

    appPushLogDate(FALSE);

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Boson condensation(Phi) (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("{\n"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        appGeneral(_T("{"));
        for (INT j = 0; j < m_lstElement[i].Num(); ++j)
        {
            appGeneral(_T("%f%s"), m_lstElement[i][j], (j != (m_lstElement[i].Num() - 1)) ? _T(", ") : _T(""));
        }
        appGeneral(_T("}"));
        if ((i + 1) == m_uiConfigurationCount)
        {
            appGeneral(_T("\n}"));
        }
        else
        {
            appGeneral(_T(",\n"));
        }
    }

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Boson condensation(Phi^2) (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        appGeneral(_T("%f, "), RealResAtI(i));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average(Phi) = "));
    for (INT j = 0; j < m_lstAverageElement.Num(); ++j)
    {
        appGeneral(_T("%f, "), m_lstAverageElement[j]);
    }

    appGeneral(_T("\n ----------- average(Phi^2) = "));
    appGeneral(_T("%f, "), GetAverageRealRes());

    appGeneral(_T(" --------\n==========================================================================\n"));
    appPopLogDate();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================