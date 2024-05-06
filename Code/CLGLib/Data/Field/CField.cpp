//=============================================================================
// FILENAME : CField.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CField::CField() 
    : CBase()
    , m_pOwner(NULL)
    , m_fLength(F(1.0))
    , m_pPool(NULL)
{
    
}

void CField::Return()
{
    assert(NULL != m_pPool);
    m_pPool->Return(this);
}

CCString CField::SaveToFile(const CCString& fileName, EFieldFileType eType) const
{
    switch (eType)
    {
    case EFFT_CLGBin:
        {
            UINT uiSize = 0;
            BYTE* byToSave = CopyDataOut(uiSize);
            appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
            CCString MD5 = CLGMD5Hash(byToSave, uiSize);
            free(byToSave);
            return MD5;
        }
    case EFFT_CLGBinCompressed:
        {
            return SaveToCompressedFile(fileName);
        }
    case EFFT_CLGBinFloat:
        {
            UINT uiSize = 0;
            BYTE* byToSave = CopyDataOutFloat(uiSize);
            appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
            CCString MD5 = CLGMD5Hash(byToSave, uiSize);
            free(byToSave);
            return MD5;
        }
    case EFFT_CLGBinDouble:
        {
            UINT uiSize = 0;
            BYTE* byToSave = CopyDataOutDouble(uiSize);
            appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
            CCString MD5 = CLGMD5Hash(byToSave, uiSize);
            free(byToSave);
            return MD5;
        }
    default:
        break;
    }

    appCrucial(_T("Save for this type not implemented: CFieldGaugeSU3 : %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
    return _T("Not supported");
}

CCString CField::GetInfos(const CCString& tab) const
{
    CCString sRet = CBase::GetInfos(tab);
    sRet = sRet + tab + _T("FieldId : ") + appToString(m_byFieldId) + _T("\n");
    sRet = sRet + tab + _T("GaugeFields : ") + appToString(m_byGaugeFieldIds) + _T("\n");
    sRet = sRet + tab + _T("BosonFields : ") + appToString(m_byBosonFieldIds) + _T("\n");
    return sRet;
}

CFieldFermion::CFieldFermion()
: CField()
{
    m_uiLinkeCount = _HC_Volume * _HC_Dir;
    m_uiSiteCount = _HC_Volume;
}

CFieldMatrixOperation* CFieldMatrixOperation::Create(EFieldType ef)
{
    if (ef == EFT_FermionWilsonSquareSU3)
    {
        return new CFieldMatrixOperationWilsonSquareSU3();
    }

    appCrucial(_T("Matrix operation for field type %s not implemented!\n"), __ENUM_TO_STRING(EFieldType, ef).c_str());
    return NULL;
}

UBOOL CFieldFermion::RationalApproximation(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const class CRatinalApproximation* pRational)
{
    if (NULL == pRational)
    {
        return FALSE;
    }
    if (pRational->m_uiDegree < 1)
    {
        return FALSE;
    }
    CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
    if (NULL == solver)
    {
        return FALSE;
    }
    TArray<CField*> solutions;
    TArray<CLGComplex> shifts;
    for (UINT i = 0; i < pRational->m_uiDegree; ++i)
    {
        CField* pPooled = appGetLattice()->GetPooledFieldById(m_byFieldId);
        solutions.AddItem(pPooled);
        shifts.AddItem(_make_cuComplex(pRational->m_lstB[i], F(0.0)));
    }

    solver->Solve(solutions, shifts, this, gaugeNum, bosonNum, gaugeFields, pBoson, op);

    ScalarMultply(pRational->m_fC);

    for (UINT i = 0; i < pRational->m_uiDegree; ++i)
    {
        Axpy(pRational->m_lstA[i], solutions[i]);
        solutions[i]->Return();
    }
    return TRUE;
}

void CFieldFermionKS::InitialOtherParameters(CParameters& params)
{
    CFieldFermion::InitialOtherParameters(params);

    params.FetchValueReal(_T("Mass"), m_f2am);
    if (m_f2am < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionKS: Mass is nearly 0!\n"));
    }

    INT iEachEta = 0;
    params.FetchValueINT(_T("EachSiteEta"), iEachEta);
    m_bEachSiteEta = (0 != iEachEta);

    TArray<Real> coeffs;
    params.FetchValueArrayReal(_T("MD"), coeffs);
    m_rMD.Initial(coeffs);

    params.FetchValueArrayReal(_T("MC"), coeffs);
    m_rMC.Initial(coeffs);

    //params.FetchValueArrayReal(_T("EN"), coeffs);
    //m_rEN.Initial(coeffs);

    if (NULL != m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(m_pMDNumerator));
    }
    checkCudaErrors(cudaMalloc((void**)&m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));
    Real* hostNumerator = (Real*)appAlloca(sizeof(Real) * m_rMD.m_uiDegree);
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        hostNumerator[i] = m_rMD.m_lstA[i];
    }
    checkCudaErrors(cudaMemcpy(m_pMDNumerator, hostNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyHostToDevice));

}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================