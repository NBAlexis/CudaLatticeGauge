//=============================================================================
// FILENAME : CFieldGauge.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/03/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CFieldGauge::CFieldGauge()
    : CField()
    , m_uiLinkeCount(_HC_Volume * _HC_Dir)
{

}

CFieldGauge::~CFieldGauge()
{

}


void CFieldGauge::CopyTo(CField* U) const
{
    CFieldGauge* pFieldGauge = dynamic_cast<CFieldGauge*>(U);
    if (NULL == pFieldGauge)
    {
        return;
    }

    CField::CopyTo(U);

    pFieldGauge->m_uiLinkeCount = m_uiLinkeCount;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================