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

    CGaugeSmearing() : m_pOwner(NULL) {  }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;
    virtual void GaugeSmearing(class CFieldGauge* pGauge, CFieldGauge* pStaple) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;

    class CLatticeData* m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARING_H_

//=============================================================================
// END OF FILE
//=============================================================================