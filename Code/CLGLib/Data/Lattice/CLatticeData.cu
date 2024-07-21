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

#if !_CLG_DOUBLEFLOAT
DOUBLE CLGAPI CCommonData::m_fBeta = 0.0;
#else
Real CLGAPI CCommonData::m_fBeta = F(0.0);
#endif
Real CLGAPI CCommonData::m_fKai = F(0.0);
Real CLGAPI CCommonData::m_fG = F(0.0);
Real CLGAPI CCommonData::m_fShiftedMass = F(0.0);
//SSmallInt4 CLGAPI CCommonData::m_sCenter = SSmallInt4(0,0,0,0);
UBOOL CLGAPI CCommonData::m_bStoreStaple = TRUE;
UBOOL CLGAPI CCommonData::m_bStoreLastSolution = TRUE;
UBOOL CLGAPI CCommonData::m_bStochasticGaussian = FALSE;
UINT CLGAPI CCommonData::m_uiMaxThreadPerBlock = 0;

Real CLGAPI CCommonData::m_fBz = F(0.0);
Real CLGAPI CCommonData::m_fEz = F(0.0);

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_uiRandomType(0)
    , m_uiRandomSeed(0)

    , m_pAphys(NULL)
    , m_pUpure(NULL)
    , m_pFieldCache(NULL)

    , m_pUpdator(NULL)
    , m_pMeasurements(NULL)
    , m_pDeviceRandom(NULL)
    , m_pIndex(NULL)
    , m_pIndexCache(NULL)

    , m_pGaugeSmearing(NULL)
    , m_pGaugeFixing(NULL)
{
    m_pFieldCache = new CFieldCache();
    memset(m_pFermionSolver, 0, sizeof(CSLASolver*) * kMaxFieldCount);
    memset(m_pFermionMultiShiftSolver, 0, sizeof(CMultiShiftSolver*) * kMaxFieldCount);
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
    appSafeDelete(m_pAphys);
    appSafeDelete(m_pUpure);
    appSafeDelete(m_pIndexCache);
    appSafeDelete(m_pRandom);
    appSafeDelete(m_pMeasurements);
    for (BYTE byField = 0; byField < kMaxFieldCount; ++byField)
    {
        appSafeDelete(m_pFermionSolver[byField]);
        appSafeDelete(m_pFermionMultiShiftSolver[byField]);
    }
    appSafeDelete(m_pGaugeSmearing);
    appSafeDelete(m_pGaugeFixing);
}

void CLatticeData::CreateFermionSolver(const CCString& sSolver, const CParameters& param, const CField* pFermionField, BYTE byFieldId)
{
    CBase* pSolver = appCreate(sSolver);
    m_pFermionSolver[byFieldId] = dynamic_cast<CSLASolver*>(pSolver);
    if (NULL == m_pFermionSolver[byFieldId])
    {
        appCrucial(_T("Create Fermion Solver %s failed!\n"), sSolver.c_str());
        _FAIL_EXIT;
    }
    m_pFermionSolver[byFieldId]->Configurate(param);
    m_pFermionSolver[byFieldId]->AllocateBuffers(pFermionField);

    appGeneral(_T("Create sparse linear algebra solver: %s \n"), sSolver.c_str());
}

void CLatticeData::CreateMultiShiftSolver(const CCString& sSolver, const CParameters& param, const CField* pFermionField, BYTE byFieldId)
{
    CBase* pSolver = appCreate(sSolver);
    m_pFermionMultiShiftSolver[byFieldId] = dynamic_cast<CMultiShiftSolver*>(pSolver);
    if (NULL == m_pFermionMultiShiftSolver[byFieldId])
    {
        appCrucial(_T("Create Fermion Solver %s failed!\n"), sSolver.c_str());
        _FAIL_EXIT;
    }
    m_pFermionMultiShiftSolver[byFieldId]->Configurate(param);
    m_pFermionMultiShiftSolver[byFieldId]->AllocateBuffers(pFermionField);

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

    appGeneral(_T("Create field pool, with field id: %d and count: %d\n"), byFieldId, uiCount);
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

CField* CLatticeData::GetPooledCopy(const CField* pField)
{
    CField* pooled = GetPooledFieldById(pField->m_byFieldId);
    if (NULL != pooled)
    {
        pField->CopyTo(pooled);
    }
    return pooled;
}

void CLatticeData::ReCopyPooled() const
{
    for (INT i = 0; i < m_pFieldPools.Num(); ++i)
    {
        m_pFieldPools[i]->ReCopyAll();
    }
}

void CLatticeData::ReCopyPooled(BYTE byId) const
{
    if (!m_pFieldPoolMap.Exist(byId))
    {
        appCrucial(_T("Get Pooled field failed!\n"));
        return;
    }
    const CFieldPool* pool = m_pFieldPoolMap.GetAt(byId);
    pool->ReCopyAll();
}

void CLatticeData::OnUpdatorConfigurationAccepted(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pAcceptGauge, const CFieldBoson* const* pAcceptBoson, const CFieldGauge* const* pCorrespondingStaple) const
{
    //accept gauge already copy to m_pGaugeField in Updator, maybe change this behavour in the furture.
    if (NULL != m_pMeasurements)
    {
        m_pMeasurements->OnConfigurationAccepted(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, pCorrespondingStaple);
    }
}

void CLatticeData::OnUpdatorFinished(UBOOL bMeasured, UBOOL bReport) const
{
    if (NULL != m_pMeasurements && bMeasured)
    {
        m_pMeasurements->OnUpdateFinished(bReport);
    }
}

void CLatticeData::FixAllFieldBoundary() const
{
    if (NULL == m_pIndex)
    {
        return;
    }

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

void CLatticeData::SetFieldBoundaryCondition(BYTE byFieldId, const SBoundCondition& bc) const
{
    m_pIndex->m_pBoundaryCondition->SetFieldSpecificBc(byFieldId, bc);
}

void CLatticeData::SetAPhys(const CFieldGauge* pUatCoulomb)
{
    if (NULL == m_pAphys)
    {
        m_pAphys = dynamic_cast<CFieldGauge*>(pUatCoulomb->GetCopy());
    }
    else
    {
        pUatCoulomb->CopyTo(m_pAphys);
    }
    m_pAphys->TransformToIA();
}

void CLatticeData::SetAPure(const CFieldGauge* pUnow)
{
    if (NULL == m_pAphys)
    {
        appCrucial(_T("CLatticeData: A phys is undefined, need to define A phys first!\n"));
        return;
    }
    if (NULL == m_pUpure)
    {
        m_pUpure = dynamic_cast<CFieldGauge*>(pUnow->GetCopy());
    }
    else
    {
        pUnow->CopyTo(m_pUpure);
    }
    m_pUpure->TransformToIA();
    m_pUpure->AxpyMinus(m_pAphys);
    m_pUpure->TransformToU();
}

UINT CLatticeData::GetDefaultSUN() const
{
    if (m_pGaugeField.Num() > 0 && NULL != m_pGaugeField[0])
    {
        return m_pGaugeField[0]->MatrixN();
    }
    return 3;
}

CCString CLatticeData::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    CCString sInfos;
    BYTE realByte[8];
    Real testReal = F(-1.2345);
    memset(realByte, 0, 8);
    memcpy(realByte, &testReal, sizeof(Real));
    CCString sRealByte;
    for (UINT i = 0; i < 8; ++i)
    {
        sRealByte += appToString(realByte[i]) + _T(", ");
    }
    sInfos.Format(_T("Version : %d\n"), appVersion());
    sInfos.Format(_T("LatticeSize : [%d, %d, %d, %d]\n"), _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt);
    sRet = sTab + sInfos;
    sInfos.Format(_T("Center : [%d, %d, %d, %d]\n"), _HC_Centerx, _HC_Centery, _HC_Centerz, _HC_Centert);
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("Random : %s\n"), __ENUM_TO_STRING(ERandom, static_cast<ERandom>(m_uiRandomType)).c_str());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("sizeof(Real) : %d and -1.2345 is %s\n"), sizeof(Real), sRealByte.c_str());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("RandomSeed : %d\n"), m_uiRandomSeed);
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("IndexType : %s\n"), NULL == m_pIndex ? _T("None") : m_pIndex->GetClass()->GetName());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("BoundaryCondition : %s\n"), (NULL == m_pIndex || NULL == m_pIndex->m_pBoundaryCondition) ? _T("None") : m_pIndex->m_pBoundaryCondition->GetClass()->GetName());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("A field Definition (is Log(U) or U.TA()) : %d\n"), _HC_ALog);
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("SUN : %d\n"), GetDefaultSUN());
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("Beta : %f\n"), CCommonData::m_fBeta);
    sRet = sRet + sTab + sInfos;
    sInfos.Format(_T("H diff : %f\n"), m_pUpdator->GetLastHDiff());
    sRet = sRet + sTab + sInfos;

    for (INT i = 0; i < kMaxFieldCount; ++i)
    {
        if (NULL != m_pFermionSolver[i])
        {
            sRet = sRet + sTab + _T("Solver : \n");
            sRet = sRet + m_pFermionSolver[i]->GetInfos(sTab + _T("    "));
        }

        if (NULL != m_pFermionMultiShiftSolver[i])
        {
            sRet = sRet + sTab + _T("Multi Shift Solver : \n");
            sRet = sRet + m_pFermionMultiShiftSolver[i]->GetInfos(sTab + _T("    "));
        }
    }

    if (NULL != m_pUpdator)
    {
        sRet = sRet + sTab + _T("Updator : \n");
        sRet = sRet + m_pUpdator->GetInfos(sTab + _T("    "));
    }

    sInfos.Format(_T("All Field Count : %d\n"), m_pOtherFields.Num());
    sRet = sRet + sTab + sInfos;
    for (INT i = 0; i < m_pOtherFields.Num(); ++i)
    {
        sRet = sRet + sTab + _T("Field") + appToString(i) + _T(" : \n");
        sRet = sRet + m_pOtherFields[i]->GetInfos(sTab + _T("    "));
    }
    sInfos.Format(_T("ActionCount : %d\n"), m_pActionList.Num());
    sRet = sRet + sTab + sInfos;
    for (INT i = 0; i < m_pActionList.Num(); ++i)
    {
        sRet = sRet + sTab + _T("Action") + appToString(i) + _T(" : \n");
        sRet = sRet + m_pActionList[i]->GetInfos(sTab + _T("    "));
    }

    if (NULL == m_pMeasurements)
    {
        sRet = sRet + sTab + _T("MeasureCount : 0\n");
    }
    else
    {
        sInfos.Format(_T("MeasureCount : %d\n"), m_pMeasurements->m_lstAllMeasures.Num());
        sRet = sRet + sTab + sInfos;
        for (INT i = 0; i < m_pMeasurements->m_lstAllMeasures.Num(); ++i)
        {
            sRet = sRet + sTab + _T("Measure") + appToString(i) + _T(" : \n");
            sRet = sRet + m_pMeasurements->m_lstAllMeasures[i]->GetInfos(sTab + _T("    "));
        }
    }

    if (NULL != m_pGaugeSmearing)
    {
        sRet = sRet + sTab + _T("Gauge Smearing : \n");
        sRet = sRet + m_pGaugeSmearing->GetInfos(sTab + _T("    "));
    }

    if (NULL != m_pGaugeFixing)
    {
        sRet = sRet + sTab + _T("Gauge Fixing : \n");
        sRet = sRet + m_pGaugeFixing->GetInfos(sTab + _T("    "));
    }
    return sRet;
}

INT CLatticeData::GetGaugeFieldIndexById(INT num, const CFieldGauge* const* gaugeFields, BYTE byFieldId)
{
    for (INT i = 0; i < num; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == byFieldId)
        {
            return i;
        }
    }
    return -1;
}

INT CLatticeData::GetBosonFieldIndexById(INT num, const CFieldBoson* const* bosonFields, BYTE byFieldId)
{
    for (INT i = 0; i < num; ++i)
    {
        if (NULL != bosonFields[i] && bosonFields[i]->m_byFieldId == byFieldId)
        {
            return i;
        }
    }
    return -1;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

