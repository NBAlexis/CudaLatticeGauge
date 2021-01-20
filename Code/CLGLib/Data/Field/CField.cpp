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
    }

    appCrucial(_T("Save for this type not implemented: CFieldGaugeSU3 : %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
    return _T("Not supported");
}

CFieldFermion::CFieldFermion()
: CField()
, m_byEvenFieldId(-1)
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

UBOOL CFieldFermion::RationalApproximation(EFieldOperator op, const CField* pGauge, const class CRatinalApproximation* pRational)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
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

    solver->Solve(solutions, shifts, this, pFieldSU3, op);

    ScalarMultply(pRational->m_fC);

    for (UINT i = 0; i < pRational->m_uiDegree; ++i)
    {
        Axpy(pRational->m_lstA[i], solutions[i]);
        solutions[i]->Return();
    }
    return TRUE;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================