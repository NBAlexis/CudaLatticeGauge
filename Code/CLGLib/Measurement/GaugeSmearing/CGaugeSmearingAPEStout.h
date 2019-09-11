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

__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingAPEStout)

class CLGAPI CGaugeSmearingAPEStout : public CGaugeSmearing
{
    __CLGDECLARE_CLASS(CGaugeSmearingAPEStout)
public:

    CGaugeSmearingAPEStout() 
        : CGaugeSmearing()
        , m_fRhoS(F(0.1))
        , m_fRhoT(F(0.0))
        , m_bHasT(FALSE)
        , m_uiIterate(1) {  }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeSmearing(class CFieldGauge* pGauge, CFieldGauge* pStaple) override;
    CCString GetInfos(const CCString& sTab) const override;

    Real m_fRhoS;
    Real m_fRhoT;
    UBOOL m_bHasT;
    UINT m_uiIterate;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGAPESTOUT_H_

//=============================================================================
// END OF FILE
//=============================================================================