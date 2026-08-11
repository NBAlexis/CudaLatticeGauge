//=============================================================================
// FILENAME : CStapleCache.h
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [05/17/2025 nbale]
//=============================================================================
#include "GaugeSmearing/CGaugeSmearing.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"

#ifndef _CSTAPLECACHE_H_
#define _CSTAPLECACHE_H_

__BEGIN_NAMESPACE

inline class CGaugeSmearing* appGetGaugeSmearing(BYTE byFieldId);

class CLGAPI CStapleCache : public CBase
{
    //For REM (rotation + EM), the cached rotation links must be built with the
    //U(1) phase of a given charge, so the cache needs to know the charge and
    //the U(1) field id. KSTR does not need this (charge = 0, no phase).
public:

    CStapleCache()
        : m_byFieldId(0)
        , m_bCacheStaple(FALSE)
        , m_bCachePlaqutte(FALSE)
        , m_bCacheFmunu(TRUE)
        , m_bCacheRotationLink(FALSE)
        , m_fCacheRotationLinkCharge(F(0.0))
        , m_byCacheRotationLinkPhaseFieldId(0)

        , m_bFermionUpdate(TRUE)
        , m_bGaugeUpdate(FALSE)
        , m_bUseEffectiveGauge(TRUE)
    {

    }

    virtual void InitialBuffers(BYTE byFieldId) = 0;
    virtual void Cache(const CFieldGauge* pGauge, ECacheCall eCC) = 0;
    virtual const void* const* GetStaples() const = 0;
    virtual const void* const* GetPlaquttes() const = 0;
    virtual const void* GetFmunu() const = 0;
    virtual const void* GetRotationBuffer() const = 0;
    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);

    UBOOL FermionUpdate() const
    {
        return m_bFermionUpdate;
    }
    UBOOL GaugeUpdate() const
    {
        return m_bGaugeUpdate;
    }

    BYTE m_byFieldId;
    UBOOL m_bCacheStaple;
    UBOOL m_bCachePlaqutte;
    UBOOL m_bCacheFmunu;
    UBOOL m_bCacheRotationLink;
    Real m_fCacheRotationLinkCharge;
    BYTE m_byCacheRotationLinkPhaseFieldId;
    UBOOL m_bFermionUpdate;
    UBOOL m_bGaugeUpdate;

    UBOOL m_bUseEffectiveGauge;
};

template<class FieldType>
class __DLL_EXPORT CStapleCacheT : public CStapleCache
{
public:

    CStapleCacheT()
        : CStapleCache()
        , m_pDeviceStaplePtr(NULL)
        , m_pDevicePlaqPtr(NULL)
        , m_pDeviceFmunuPtr(NULL)
        , m_pDeviceGaugeRotationLinksPtr(NULL)
    {

    }

    //to remove dependency on cuda, move the checkCudaErrors to .cpp
    ~CStapleCacheT();
    void InitialBuffers(BYTE byFieldId) override;

    void Cache(const CFieldGauge* pGauge, ECacheCall eCC) override
    {
        //============ before smearing, if m_bCacheOrignalGauge = TRUE, cache orignal gauge
        if (ECC_BeforeGaugeUpdate == eCC)
        {
            //This may be helpful in calculation of gauge force
            //I don't see the motivation for cache original gauge, so not implemented

            return;
        }

        //============ before smearing, if m_bCacheOrignalGauge = TRUE, cache orignal gauge
        if (ECC_BeforeFermionUpdateBeforeSmearing == eCC || ECC_BeforeAllUpdateBeforeSmearing == eCC)
        {
            //This may be helpful to calculate gauge smearing
            //I don't see the motivation for cache original gauge, so not implemented

            return;
        }

        //============ after smearing, if m_bCacheEffectiveGauge = TRUE and NULL != effective, cache effective gauge, if m_bCacheEffectiveGauge = TRUE and NULL == effective, cache effective gauge
        if (ECC_BeforeFermionUpdateAfterSmearing == eCC || ECC_BeforeAllUpdateAfterSmearing == eCC)
        {
            if ((NULL == m_pDeviceStaplePtr && m_bCacheStaple)
             || (NULL == m_pDevicePlaqPtr && m_bCachePlaqutte)
             || (NULL == m_pDeviceFmunuPtr && m_bCacheFmunu)
             || (NULL == m_pDeviceGaugeRotationLinksPtr && m_bCacheRotationLink)
                )
            {
                InitialBuffers(pGauge->m_byFieldId);
            }
            const FieldType* pGaugeT = (m_bUseEffectiveGauge && NULL != appGetGaugeSmearing(pGauge->m_byFieldId)) ? dynamic_cast<const FieldType*>(appGetGaugeSmearing(pGauge->m_byFieldId)->GetEffectiveGauge()) : dynamic_cast<const FieldType*>(pGauge);
            if (m_bUseEffectiveGauge && NULL != appGetGaugeSmearing(pGauge->m_byFieldId))
            {
                appParanoiac(_T("CStapleCacheT::Cache using effective gauge\n"));
            }

            if (m_bCacheStaple)
            {
                appParanoiac(_T("CStapleCacheT::Cache staple\n"));
                pGaugeT->CalculateAllStaples(m_pDeviceStaplePtr);
            }
            if (m_bCachePlaqutte)
            {
                appParanoiac(_T("CStapleCacheT::Cache plaquttes\n"));
                pGaugeT->CalculateAllPlaquttes(m_pDevicePlaqPtr);
            }
            if (m_bCacheFmunu)
            {
                appParanoiac(_T("CStapleCacheT::Cache Fmunu\n"));
                pGaugeT->CalculateFmunu(m_pDeviceFmunuPtr);
            }
            if (m_bCacheRotationLink)
            {
                appParanoiac(_T("CStapleCacheT::Cache rotation link\n"));
                if (0 != m_byCacheRotationLinkPhaseFieldId)
                {
                    //REM (rotation + EM): the rotation links carry the U(1) phase
                    //of the given charge, so generate the cache with phase.
                    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byCacheRotationLinkPhaseFieldId));
                    if (NULL == pU1)
                    {
                        appCrucial(_T("CStapleCacheT::Cache rotation link: phase field not found!\n"));
                        return;
                    }
                    pGaugeT->CacheRotationKSBuffer(m_pDeviceGaugeRotationLinksPtr, TRUE, m_fCacheRotationLinkCharge, (const Real*)pU1->m_pDeviceData);
                }
                else
                {
                    pGaugeT->CacheRotationKSBuffer(m_pDeviceGaugeRotationLinksPtr, FALSE, F(0.0), NULL);
                }
            }
        }
    }

    const void* const * GetStaples() const override
    {
        return (const void* const*)m_pDeviceStaplePtr;
    }

    const void* const* GetPlaquttes() const override
    {
        return (const void* const*)m_pDevicePlaqPtr;
    }

    const void* GetFmunu() const override
    {
        return (const void*)m_pDeviceFmunuPtr;
    }

    const void* GetRotationBuffer() const override
    {
        return (const void*)m_pDeviceGaugeRotationLinksPtr;
    }

    TArray<FieldType*> m_pPooledFields;
    typename FieldType::_Gauge** m_pDeviceStaplePtr;
    typename FieldType::_Gauge** m_pDevicePlaqPtr;

    /**
    * 6 buffers:
    * x-y
    * x-z
    * x-t
    * y-z
    * y-t
    * z-t
    * 
    * Fmunu = (Pmunu - Pmunu^+) / 8i, (arXiv:1311.6312)
    * This buffer stores (Pmunu - Pmunu^+)
    *
    * Pmunu  = + (mu, nu) + (nu, -mu) - (-nu, -mu) - (mu, -nu) 10.1016/0550-3213(85)90002-1
    * Pmunu+ = + (nu, mu) + (-mu, nu) - (-mu, -nu) - (-nu, mu)
    * Pmunu - Pmunu+  = Qmunu - Qmunu+
    * Qmunu = + (mu, nu) - (-mu, nu) + (-mu, -nu) - (mu, -nu)
    */
    typename FieldType::_Gauge* m_pDeviceFmunuPtr;

    typename FieldType::_Gauge* m_pDeviceGaugeRotationLinksPtr;

    //Improve-1 (multi-GPU-improve1.md 3.4.1/3.4.2, I6): owner-side halo
    //registration. The staple/plaquette device POINTER ARRAYS get buffer-set
    //handles whose members are the pooled fields their entries point to (the
    //members are 1-hop-stencil read through the arrays). The flat Fmunu and
    //rotation-link buffers are consumed site-locally only (verified against
    //_kernelDFermionWilsonSquareCloverSU3 and the _kernelDFermionKS_PR_Cached_*
    //family, which index them strictly at uiSiteIndex) -- LocalOnly handles,
    //no halo traffic, writes still versioned.
    CHaloBufferSetHandle m_StapleArraySet;
    CHaloBufferSetHandle m_PlaqArraySet;
    CHaloBufferHandle m_FmunuHandle;
    CHaloBufferHandle m_RotationLinksHandle;
};

__CLG_REGISTER_HELPER_HEADER(CStapleCacheSU3)
class CLGAPI CStapleCacheSU3 : public CStapleCacheT<CFieldGaugeSU3>
{
    __CLGDECLARE_CLASS(CStapleCacheSU3)
};

__END_NAMESPACE

#endif //#ifndef _CSTAPLECACHE_H_

//=============================================================================
// END OF FILE
//=============================================================================