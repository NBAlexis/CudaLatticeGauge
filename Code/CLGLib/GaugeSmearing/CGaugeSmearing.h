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
    , m_uiIterate(1)
    , m_bCalledWhenUpdate(FALSE)
    , m_pEffecitveGauge(NULL)
    {
        
    }

    ~CGaugeSmearing()
    {
        //if (NULL != m_pEffecitveGauge)
        //{
        //    m_pEffecitveGauge->Return();
        //}
        appSafeDelete(m_pTmpStaple);
    }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);

    //This is used in measurement (and update via GaugeSmearingC)
    //For some gauge smearing to work, pOrignal must be filled, while for others, no need to do so
    //pGauge = pOrignal is garentted and is free to change
    virtual void GaugeSmearing(class CFieldGauge* pGauge, const CFieldGauge* pOrignal, class CFieldGauge* pStaple, UBOOL bProject = TRUE) = 0;

    //This is used for update
    virtual void GaugeSmearingC(const class CFieldGauge* pGauge)
    {
        if (NULL == m_pEffecitveGauge)
        {
            m_pEffecitveGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        }
        else
        {
            pGauge->CopyTo(m_pEffecitveGauge);
        }
        GaugeSmearing(m_pEffecitveGauge, pGauge, NULL);
    }

    CCString GetInfos(const CCString& sTab) const override;

    virtual const CFieldGauge* GetEffectiveGaugeLevel1() const
    {
        return NULL;
    }

    virtual const CFieldGauge* GetNaikLink() const
    {
        return NULL;
    }

    virtual const CFieldGauge* GetEffectiveGauge() const
    {
        return m_pEffecitveGauge;
    }

    UBOOL CalledWhenUpdate() const 
    {
        return m_bCalledWhenUpdate;
    }

    //The fat link with T direction is just staple, this calculate fat link without T direction
    virtual void CalculateSpatialFatLink(const class CFieldGauge* pGauge, class CFieldGauge* pFatlink) const;
    virtual void DerivateOnU(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const
    {
        appCrucial(_T("This gauge smearing does not support calculation of force!\n"));
    }

    class CLatticeData* m_pOwner;
    class CFieldGauge* m_pTmpStaple;
    UBOOL m_bHasT;
    UINT m_uiIterate;
    BYTE m_byFieldId;
    UBOOL m_bCalledWhenUpdate;
    class CFieldGauge* m_pEffecitveGauge;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARING_H_

//=============================================================================
// END OF FILE
//=============================================================================