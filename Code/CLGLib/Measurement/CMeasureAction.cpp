//=============================================================================
// FILENAME : CMeasureAction.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [10/06/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CMeasureAction.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAction)

void CMeasureAction::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("FermiomFieldCount"), iValue);
    m_iFermionFieldCount = iValue > 1 ? iValue : 1;
}

void CMeasureAction::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple)
{
    DOUBLE plaqutteEneregy = 0.0;

    for (INT i = 0; i < appGetLattice()->m_pActionList.Num(); ++i)
    {
        if (appGetLattice()->m_pActionList[i]->IsFermion())
        {
            DOUBLE fToBeAdd = 0.0;
            for (UINT j = 0; j < m_iFermionFieldCount; ++j)
            {
                appGetLattice()->m_pActionList[i]->PrepareForHMC(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, 0);
                fToBeAdd += appGetLattice()->m_pActionList[i]->Energy(FALSE, gaugeNum, bosonNum, tensor2Num, pAcceptGauge, pAcceptBoson, tensor2Fields, NULL);
            }
            plaqutteEneregy += (fToBeAdd / m_iFermionFieldCount);
        }
        else
        {
            plaqutteEneregy += appGetLattice()->m_pActionList[i]->Energy(FALSE, gaugeNum, bosonNum, tensor2Num, pAcceptGauge, pAcceptBoson, tensor2Fields, NULL);
        }
    }

#if !_CLG_DOUBLEFLOAT
    UpdateRealResult(static_cast<Real>(plaqutteEneregy));
#else
    UpdateRealResult(plaqutteEneregy);
#endif
    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("Action"), GetLastRealRes());
    }
    appParanoiac(_T(" === Action Energy Measured === energy = %f\n"), plaqutteEneregy);
}

void CMeasureAction::Report()
{
    Average();
    appGeneral(_T(" === Action Averaged === energy = %f\n\n"), GetAverageRealRes());
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================