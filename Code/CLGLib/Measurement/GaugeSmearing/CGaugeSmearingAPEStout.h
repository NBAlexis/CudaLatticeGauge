//=============================================================================
// FILENAME : CGaugeSmearingAPEStout.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================

#ifndef _CGAUGESMEARINGAPESTOUT_H_
#define _CGAUGESMEARINGAPESTOUT_H_

__BEGIN_NAMESPACE

class CLGAPI CGaugeSmearingAPEStout : public CGaugeSmearing
{
    __CLGDECLARE_CLASS(CGaugeSmearingAPEStout)
public:

    CGaugeSmearingAPEStout() : CGaugeSmearing(), m_fRhoS(F(0.1)), m_fRhoT(F(0.0)), m_bHasT(FALSE), m_uiIterate(1) {  }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);
    virtual void GaugeSmearing(class CFieldGauge* pGauge, CFieldGauge* pStaple);
    virtual CCString GetInfos(const CCString& sTab) const;

    UINT m_uiIterate;
    Real m_fRhoS;
    Real m_fRhoT;
    UBOOL m_bHasT;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGAPESTOUT_H_

//=============================================================================
// END OF FILE
//=============================================================================