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
        , m_fAlphaLeft(F(1.0))
        , m_fAlphaRight(F(0.1))
        , m_bCMProj(FALSE)
        , m_byProjIterate(6)
    {  }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeSmearing(class CFieldGauge* pGauge, CFieldGauge* pStaple) override;
    CCString GetInfos(const CCString& sTab) const override;

    Real m_fAlphaLeft;
    Real m_fAlphaRight;
    UBOOL m_bCMProj;
    BYTE m_byProjIterate;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGAPEPROJ_H_

//=============================================================================
// END OF FILE
//=============================================================================