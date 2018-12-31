//=============================================================================
// FILENAME : CFieldFermion.h
// 
// DESCRIPTION:
// This is the class for all fermion fields
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDFERMION_H_
#define _CFIELDFERMION_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldFermion : public CField
{
public:
    CFieldFermion();

    virtual void PrepareForHMC(const CFieldGauge* pGauge) = 0;
    virtual void CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce) = 0;

protected:

    UINT m_uiLinkeCount;
    UINT m_uiSiteCount;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================