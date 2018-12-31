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

CField::CField() : CBase(), m_pOwner(NULL)
{
    m_uiThreadCount = _HC_DecompLx * _HC_DecompLy * _HC_DecompLz;
}

CFieldGauge::CFieldGauge() : CField()
{
    m_uiLinkeCount = _HC_Volumn * _HC_Dir;
}

CFieldFermion::CFieldFermion() : CField()
{
    m_uiLinkeCount = _HC_Volumn * _HC_Dir;
    m_uiSiteCount = _HC_Volumn;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================