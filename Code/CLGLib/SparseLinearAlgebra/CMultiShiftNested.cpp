//=============================================================================
// FILENAME : CMultiShiftNested.cpp
// 
// DESCRIPTION:
// This is the class for GMRES Solver
//
// REVISION:
//  [15/06/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CMultiShiftNested.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMultiShiftNested)

CMultiShiftNested::CMultiShiftNested()
    : CMultiShiftSolver()
    , m_pNestedSolver(NULL)
{

}

CMultiShiftNested::~CMultiShiftNested()
{
    appSafeDelete(m_pNestedSolver);
}

void CMultiShiftNested::Configurate(const CParameters& param)
{
    //CMultiShiftSolver::Configurate(param);

    CParameters nestedSolver;
    if (!param.FetchParameterValue(_T("NestedSolver"), nestedSolver))
    {
        appCrucial(_T("Nested Solver not specified!"));
        _FAIL_EXIT;
    }

    CCString sSolverName = _T("CSLASolverBiCGStab");
    nestedSolver.FetchStringValue(_T("SolverName"), sSolverName);
    INT byFieldId = 2;
    nestedSolver.FetchValueINT(_T("SolverForFieldId"), byFieldId);
    CField* pField = appGetLattice()->GetFieldById(static_cast<BYTE>(byFieldId));
    if (NULL == pField)
    {
        appCrucial(_T("Solver must be created for a specified field!\n"));
    }

    CBase* pSolver = appCreate(sSolverName);
    m_pNestedSolver = dynamic_cast<CSLASolver*>(pSolver);
    if (NULL == m_pNestedSolver)
    {
        appCrucial(_T("Create Fermion Solver %s failed!\n"), sSolverName.c_str());
        _FAIL_EXIT;
    }
    m_pNestedSolver->Configurate(nestedSolver);
    m_pNestedSolver->AllocateBuffers(pField);

    appGeneral(_T("Create sparse linear algebra solver: %s \n"), sSolverName.c_str());

}

void CMultiShiftNested::AllocateBuffers(const CField*)
{

}

UBOOL CMultiShiftNested::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    EFieldOperator uiM, ESolverPhase , const CField* )
{
    appPushLogDate(FALSE);
#if _CLG_DEBUG
    if (cn.Num() > 2)
    {
        appGeneral(_T("Better to use a REAL multi-shift solver!\n"));
    }
#endif

    EFieldOperator toSolve = EFO_F_D_WithMass;
    switch (uiM)
    {
    case EFO_F_D:
        toSolve = EFO_F_D_WithMass;
        break;
    case EFO_F_Ddagger:
        toSolve = EFO_F_Ddagger_WithMass;
        break;
    case EFO_F_DDdagger:
        toSolve = EFO_F_DDdagger_WithMass;
        break;
    default:
        appCrucial(_T("Operator not supported in nested shift solver!\n"));
        break;
    }
    for (INT i = 0; i < cn.Num(); ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);

        if (NULL == dynamic_cast<CFieldFermionKS*>(pFieldX[i]))
        {
            appCrucial(_T("CMultiShiftNested only support Staggered fermion\n"));
            continue;
        }

        CFieldFermionKS* pFieldKS = dynamic_cast<CFieldFermionKS*>(pFieldX[i]);

        if (NULL == pFieldKS)
        {
            appCrucial(_T("CMultiShiftNested only support Staggered fermion\n"));
            continue;
        }

        if (appAbs(cn[i].y) > F(0.0000001))
        {
            appCrucial(_T("CMultiShiftNested only support real number shift!\n"));
        }
        Real fShiftedMass;
        if (EFO_F_DDdagger_WithMass == toSolve)
        {
            fShiftedMass = pFieldKS->GetMass() * pFieldKS->GetMass() + cn[i].x;
            if (fShiftedMass < F(0.0000001))
            {
                appCrucial(_T("CMultiShiftNested only support real number shift and mass+shift > 0!\n"));
                continue;
            }
            fShiftedMass = _hostsqrt(fShiftedMass);
        }
        else
        {
            fShiftedMass = pFieldKS->GetMass() + cn[i].x;
            if (fShiftedMass < F(0.0000001))
            {
                appCrucial(_T("CMultiShiftNested only support real number shift and mass+shift > 0!\n"));
                continue;
            }
        }

        CCommonData::m_fShiftedMass = fShiftedMass;
        m_pNestedSolver->Solve(pFieldKS, pFieldB, gaugeNum, bosonNum, gaugeFields, bosonFields, toSolve, ESP_Once, NULL);
    }
    appPopLogDate();
    return TRUE;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================