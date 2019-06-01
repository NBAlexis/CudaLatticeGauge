//=============================================================================
// FILENAME : CLatticeData.cpp
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

Real CLGAPI CCommonData::m_fBeta = F(0.0);
Real CLGAPI CCommonData::m_fKai = F(0.0);
Real CLGAPI CCommonData::m_fOmega = F(0.0);
SSmallInt4 CLGAPI CCommonData::m_sCenter = SSmallInt4(0,0,0,0);
UBOOL CLGAPI CCommonData::m_bStoreStaple = TRUE;
UBOOL CLGAPI CCommonData::m_bStoreLastSolution = TRUE;

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_pGaugeField(NULL)
    , m_pUpdator(NULL)

    , m_pDeviceRandom(NULL)
    , m_pIndex(NULL)

    , m_pFermionSolver(NULL)
    , m_pMeasurements(NULL)
    , m_pFieldCache(NULL)
    , m_pIndexCache(NULL)
    , m_pGaugeSmearing(NULL)

    , m_uiRandomType(0)
    , m_uiRandomSeed(0)
{
    m_pFieldCache = new CFieldCache();    
}

CLatticeData::~CLatticeData()
{
    if (NULL != m_pUpdator)
    {
        appSafeDelete(m_pUpdator);
    }

    for (INT i = 0; i < m_pActionList.Num(); ++i)
    {
        appSafeDelete(m_pActionList[i]);
    }

    if (NULL != m_pIndex)
    {
        appSafeDelete(m_pIndex);
    }
    if (NULL != m_pDeviceRandom)
    {
        checkCudaErrors(cudaFree(m_pDeviceRandom));
        m_pDeviceRandom = NULL;
    }

    for (INT i = 0; i < m_pOtherFields.Num(); ++i)
    {
        appSafeDelete(m_pOtherFields[i]);
    }
    m_pOtherFields.RemoveAll();

    for (INT i = 0; i < m_pAllBoundaryFields.Num(); ++i)
    {
        appSafeDelete(m_pAllBoundaryFields[i]);
    }
    m_pAllBoundaryFields.RemoveAll();

    for (INT i = 0; i < m_pFieldPools.Num(); ++i)
    {
        appSafeDelete(m_pFieldPools[i]);
    }

    appSafeDelete(m_pFieldCache);
    appSafeDelete(m_pGaugeField);
    appSafeDelete(m_pIndexCache);
    appSafeDelete(m_pRandom);
    appSafeDelete(m_pMeasurements);
    appSafeDelete(m_pFermionSolver);
    appSafeDelete(m_pGaugeSmearing);
}

void CLatticeData::CreateFermionSolver(const CCString& sSolver, const CParameters& param, const CField* pFermionField)
{
    CBase* pSolver = appCreate(sSolver);
    m_pFermionSolver = dynamic_cast<CSLASolver*>(pSolver);
    if (NULL == m_pFermionSolver)
    {
        appCrucial(_T("Create Fermion Solver %s failed!\n"), sSolver.c_str());
        _FAIL_EXIT;
    }
    m_pFermionSolver->Configurate(param);
    m_pFermionSolver->AllocateBuffers(pFermionField);

    appGeneral(_T("Create sparse linear algebra solver: %s \n"), sSolver.c_str());
}

void CLatticeData::CreateFieldPool(BYTE byFieldId, UINT uiCount)
{
    if (m_pFieldPoolMap.Exist(byFieldId))
    {
        appGeneral(_T("Create field pool, but field id already exist!\n"));
        return;
    }
    CField* pField = GetFieldById(byFieldId);
    if (NULL == pField)
    {
        appCrucial(_T("Create field pool, but field cannot be found!\n"));
        return;
    }

    CFieldPool * pFieldPool = new CFieldPool(pField, uiCount);
    m_pFieldPools.AddItem(pFieldPool);
    m_pFieldPoolMap.SetAt(byFieldId, pFieldPool);
}

CField* CLatticeData::GetPooledFieldById(BYTE byId)
{
    if (!m_pFieldPoolMap.Exist(byId))
    {
        appCrucial(_T("Get Pooled field failed!\n"));
        return NULL;
    }

    return m_pFieldPoolMap[byId]->GetOne();
}

void CLatticeData::OnUpdatorConfigurationAccepted(const CFieldGauge* pAcceptGauge, const CFieldGauge* pCorrespondingStaple)
{
    //accept gauge already copy to m_pGaugeField in Updator, maybe change this behavour in the furture.
    if (NULL != m_pMeasurements)
    {
        m_pMeasurements->OnConfigurationAccepted(pAcceptGauge, pCorrespondingStaple);
    }
}

void CLatticeData::OnUpdatorFinished(UBOOL bMeasured, UBOOL bReport)
{
    if (NULL != m_pMeasurements && bMeasured)
    {
        m_pMeasurements->OnUpdateFinished(bReport);
    }
}

void CLatticeData::FixAllFieldBoundary()
{
    if (m_pIndex->NeedToFixBoundary())
    {
        for (UINT i = 0; i < kMaxFieldCount; ++i)
        {
            CField* pField = GetFieldById(static_cast<BYTE>(i));
            //only need to fix gauge, because fermion field do not change
            if (NULL != pField) 
            {
                pField->FixBoundary();
            }
        }
    }
}

void CLatticeData::SetFieldBoundaryCondition(BYTE byFieldId, const SBoundCondition& bc)
{
    m_pIndex->m_pBoundaryCondition->SetFieldSpecificBc(byFieldId, bc);
}

CCString CLatticeData::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    CCString sInfos;
    sInfos.Format(_T("LatticeSize : [%d, %d, %d, %d]\n"), _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt);
    sRet = sTab + sInfos;
    sInfos.Format(_T("Random : %s\n"), __ENUM_TO_STRING(ERandom, static_cast<ERandom>(m_uiRandomType)).c_str());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("sizeof(Real) : %d\n"), sizeof(Real));
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("RandomSeed : %d\n"), m_uiRandomSeed);
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("IndexType : %s\n"), NULL == m_pIndex ? _T("None") : m_pIndex->GetClass()->GetName());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("BoundaryCondition : %s\n"), (NULL == m_pIndex || NULL == m_pIndex->m_pBoundaryCondition) ? _T("None") : m_pIndex->m_pBoundaryCondition->GetClass()->GetName());
    sRet = sRet + sTab + sInfos;
    if (NULL != m_pFermionSolver)
    {
        sRet = sRet + sTab + _T("Solver : \n");
        sRet = sRet + m_pFermionSolver->GetInfos(sTab + _T("    "));
    }
    if (NULL != m_pUpdator)
    {
        sRet = sRet + sTab + _T("Updator : \n");
        sRet = sRet + m_pUpdator->GetInfos(sTab + _T("    "));
    }

    if (NULL != m_pGaugeField)
    {
        sRet = sRet + sTab + _T("GaugeField : \n");
        sRet = sRet + m_pGaugeField->GetInfos(sTab + _T("    "));
    }
    sInfos.Format(_T("OtherFieldCount : %d\n"), m_pOtherFields.Num());
    sRet = sRet + sTab + sInfos;
    for (INT i = 0; i < m_pOtherFields.Num(); ++i)
    {
        sRet = sRet + sTab + _T("OtherField") + appIntToString(i) + _T(" : \n");
        sRet = sRet + m_pOtherFields[i]->GetInfos(sTab + _T("    "));
    }
    sInfos.Format(_T("ActionCount : %d\n"), m_pActionList.Num());
    sRet = sRet + sTab + sInfos;
    for (INT i = 0; i < m_pActionList.Num(); ++i)
    {
        sRet = sRet + sTab + _T("OtherField") + appIntToString(i) + _T(" : \n");
        sRet = sRet + m_pActionList[i]->GetInfos(sTab + _T("    "));
    }
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

