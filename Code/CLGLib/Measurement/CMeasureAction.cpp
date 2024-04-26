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

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAction)

void CMeasureAction::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("FermiomFieldCount"), iValue);
    m_iFermionFieldCount = iValue > 1 ? iValue : 1;
}

void CMeasureAction::OnConfigurationAccepted(const CFieldGauge* pAcceptGauge, const CFieldGauge* )
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE plaqutteEneregy = 0.0;
#else
    Real plaqutteEneregy = F(0.0);
#endif

    for (INT i = 0; i < appGetLattice()->m_pActionList.Num(); ++i)
    {
        if (appGetLattice()->m_pActionList[i]->IsFermion())
        {
            //appGeneral(_T("field count:%d\n"), m_iFermionFieldCount);
#if !_CLG_DOUBLEFLOAT
            DOUBLE fToBeAdd = 0.0;
#else
            Real fToBeAdd = F(0.0);
#endif
            for (UINT j = 0; j < m_iFermionFieldCount; ++j)
            {
                appGetLattice()->m_pActionList[i]->PrepareForHMCSingleField(pAcceptGauge, 0);
                fToBeAdd += appGetLattice()->m_pActionList[i]->EnergySingleField(FALSE, pAcceptGauge, NULL);
            }
            plaqutteEneregy += (fToBeAdd / m_iFermionFieldCount);
        }
        else
        {
            plaqutteEneregy += appGetLattice()->m_pActionList[i]->EnergySingleField(FALSE, pAcceptGauge, NULL);
        }
        
    }

#if !_CLG_DOUBLEFLOAT
    UpdateRealResult(static_cast<Real>(plaqutteEneregy));
#else
    UpdateRealResult(plaqutteEneregy);
#endif
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