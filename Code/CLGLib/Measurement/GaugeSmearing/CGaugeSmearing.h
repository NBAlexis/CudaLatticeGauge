//=============================================================================
// FILENAME : CGaugeSmearing.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================

#ifndef _CGAUGESMEARING_H_
#define _CGAUGESMEARING_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EGaugeSmearingType,
    EGST_APEStout,

    EGST_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CGaugeSmearing : public CBase
{
public:

    CGaugeSmearing()
    : m_pOwner(NULL)
    , m_pTmpStaple(NULL)
    , m_bHasT(TRUE)
    , m_uiIterate(10)
    {
        
    }

    ~CGaugeSmearing()
    {
        appSafeDelete(m_pTmpStaple);
    }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);
    virtual void GaugeSmearing(class CFieldGauge* pGauge, class CFieldGauge* pStaple) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;

    //The fat link with T direction is just staple, this calculate fat link without T direction
    virtual void CalculateSpatialFatLink(const class CFieldGauge* pGauge, class CFieldGauge* pFatlink) const;

    class CLatticeData* m_pOwner;
    class CFieldGauge* m_pTmpStaple;
    UBOOL m_bHasT;
    UINT m_uiIterate;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARING_H_

//=============================================================================
// END OF FILE
//=============================================================================