//=============================================================================
// FILENAME : CGaugeFixing.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/18/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXING_H_
#define _CGAUGEFIXING_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EGaugeFixingType,
    EGFT_LandauCornell,
    EGFT_LandauCornellFFT,
    EGFT_CoulombCornellFFT,

    EGFT_ForceDWORD = 0x7fffffff,
    )


class CLGAPI CGaugeFixing : public CBase
{
public:

    CGaugeFixing()
    : m_pOwner(NULL)
    , m_fAccuracy(F(0.00000000001))
    , m_iIterate(0)
    , m_iMaxIterate(1000000)
    , m_iShowErrorStep(1000)
    {
        
    }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;

    /**
     * pResGauge should be a copy of the gauge, and will be changed.
     */
    virtual void GaugeFixing(CFieldGauge* pResGauge) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;
    virtual Real CheckRes(const CFieldGauge* pGauge) = 0;

    class CLatticeData* m_pOwner;
    Real m_fAccuracy;
    UINT m_iIterate;
    UINT m_iMaxIterate;
    UINT m_iShowErrorStep;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXING_H_

//=============================================================================
// END OF FILE
//=============================================================================