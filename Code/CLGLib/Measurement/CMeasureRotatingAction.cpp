//=============================================================================
// FILENAME : CMeasureRotatingAction.cpp
//
// DESCRIPTION:
// This is the class for rotating gauge action energy measurement
//
// REVISION:
//  [05/22/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CMeasureRotatingAction.h"
#include "Data/Action/CActionGaugePlaquetteRotatingT.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureRotatingAction)

void CMeasureRotatingAction::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = m_iActionIndex;
    param.FetchValueINT(_T("ActionIndex"), iValue);
    m_iActionIndex = iValue;
}

void CMeasureRotatingAction::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple)
{
    DOUBLE rotatingEnergy = 0.0;
    DOUBLE s0 = 0.0;
    DOUBLE s1 = 0.0;
    DOUBLE s2 = 0.0;

    CActionGaugePlaquetteRotating* pAction = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(static_cast<BYTE>(m_iActionIndex)));
    if (NULL != pAction)
    {
        rotatingEnergy = pAction->Energy(FALSE, gaugeNum, bosonNum, tensor2Num, pAcceptGauge, pAcceptBoson, tensor2Fields, NULL);
        s0 = pAction->GetS0Energy();
        s1 = pAction->GetS1Energy();
        s2 = pAction->GetS2Energy();
    }
    else
    {
        appCrucial(_T("CMeasureRotatingAction: ActionIndex %d not found or is not a rotating action\n"), m_iActionIndex);
    }

#if !_CLG_DOUBLEFLOAT
    UpdateRealResult(static_cast<Real>(rotatingEnergy));
#else
    UpdateRealResult(rotatingEnergy);
#endif
    m_lstS0.AddItem(s0);
    m_lstS1.AddItem(s1);
    m_lstS2.AddItem(s2);

    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("Energy"), GetLastRealRes());
        m_pOwner->AddOneConfigurationResult(this, _T("S0"), static_cast<DOUBLE>(s0));
        m_pOwner->AddOneConfigurationResult(this, _T("S1"), static_cast<DOUBLE>(s1));
        m_pOwner->AddOneConfigurationResult(this, _T("S2"), static_cast<DOUBLE>(s2));
    }

    appParanoiac(_T(" === Rotating Action Energy Measured === total = %f, S0 = %f, S1 = %f, S2 = %f\n"), rotatingEnergy, s0, s1, s2);
}

void CMeasureRotatingAction::Reset()
{
    CMeasure::Reset();
    m_lstS0.RemoveAll();
    m_lstS1.RemoveAll();
    m_lstS2.RemoveAll();
}

void CMeasureRotatingAction::Report()
{
    Average();
    appGeneral(_T(" === Rotating Action Averaged === total = %f\n\n"), GetAverageRealRes());

    if (m_lstS0.Num() > 0)
    {
        DOUBLE fAvgS0 = 0.0;
        for (INT i = 0; i < m_lstS0.Num(); ++i)
        {
            fAvgS0 += m_lstS0[i];
        }
        fAvgS0 = fAvgS0 / m_lstS0.Num();
        appGeneral(_T(" === S0 Averaged === %f\n"), fAvgS0);
    }

    if (m_lstS1.Num() > 0)
    {
        DOUBLE fAvgS1 = 0.0;
        for (INT i = 0; i < m_lstS1.Num(); ++i)
        {
            fAvgS1 += m_lstS1[i];
        }
        fAvgS1 = fAvgS1 / m_lstS1.Num();
        appGeneral(_T(" === S1 Averaged === %f\n"), fAvgS1);
    }

    if (m_lstS2.Num() > 0)
    {
        DOUBLE fAvgS2 = 0.0;
        for (INT i = 0; i < m_lstS2.Num(); ++i)
        {
            fAvgS2 += m_lstS2[i];
        }
        fAvgS2 = fAvgS2 / m_lstS2.Num();
        appGeneral(_T(" === S2 Averaged === %f\n\n"), fAvgS2);
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
