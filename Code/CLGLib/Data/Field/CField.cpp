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
    , m_pPool(NULL)
    , m_fLength(F(1.0))
{
    
}

void CField::Return()
{
    assert(NULL != m_pPool);
    m_pPool->Return(this);
}

CFieldFermion::CFieldFermion() : CField()
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

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================