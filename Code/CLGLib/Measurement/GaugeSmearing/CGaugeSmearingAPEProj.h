//=============================================================================
// FILENAME : CGaugeSmearingAPEProj.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================

#ifndef _CGAUGESMEARINGAPEPROJ_H_
#define _CGAUGESMEARINGAPEPROJ_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingAPEProj)

class CLGAPI CGaugeSmearingAPEProj : public CGaugeSmearing
{
    __CLGDECLARE_CLASS(CGaugeSmearingAPEProj)
public:

    CGaugeSmearingAPEProj() 
        : CGaugeSmearing()
        , m_fAlpha(F(0.1))
        , m_uiIterate(1)
        , m_byProjIterate(6) {  }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);
    virtual void GaugeSmearing(class CFieldGauge* pGauge, CFieldGauge* pStaple);
    virtual CCString GetInfos(const CCString& sTab) const;

    Real m_fAlpha;
    UINT m_uiIterate;
    BYTE m_byProjIterate;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGAPEPROJ_H_

//=============================================================================
// END OF FILE
//=============================================================================