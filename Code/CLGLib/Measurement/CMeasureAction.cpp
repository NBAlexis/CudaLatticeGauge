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
                appGetLattice()->m_pActionList[i]->PrepareForHMC(pAcceptGauge, 0);
                fToBeAdd += appGetLattice()->m_pActionList[i]->Energy(FALSE, pAcceptGauge);
            }
            plaqutteEneregy += (fToBeAdd / m_iFermionFieldCount);
        }
        else
        {
            plaqutteEneregy += appGetLattice()->m_pActionList[i]->Energy(FALSE, pAcceptGauge);
        }
        
    }

#if !_CLG_DOUBLEFLOAT
    m_lstData.AddItem(static_cast<Real>(plaqutteEneregy));
#else
    m_lstData.AddItem(plaqutteEneregy);
#endif
    appParanoiac(_T(" === Action Energy Measured === energy = %f\n"), plaqutteEneregy);
}

void CMeasureAction::Average(UINT )
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE fAdd = 0.0;
    for (INT i = 0; i < m_lstData.Num(); ++i)
    {
        fAdd += m_lstData[i];
    }
    m_fLastRealResult = static_cast<Real>(fAdd / m_lstData.Num());
    appParanoiac(_T(" === Action Averaged (%d measures) === energy = %f\n"), m_lstData.Num(), m_fLastRealResult);
#else
    Real fAdd = F(0.0);
    for (INT i = 0; i < m_lstData.Num(); ++i)
    {
        fAdd += m_lstData[i];
    }
    m_fLastRealResult = fAdd / m_lstData.Num();
    appParanoiac(_T(" === Action Averaged (%d measures) === energy = %f\n"), m_lstData.Num(), m_fLastRealResult);
#endif
}

void CMeasureAction::Report()
{
    appGeneral(_T(" === Action Averaged === energy = %f\n\n"), m_fLastRealResult);
}

void CMeasureAction::Reset()
{
    m_lstData.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================