//=============================================================================
// FILENAME : CStapleCache.cpp
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [05/17/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CStapleCache.h"

__BEGIN_NAMESPACE

void CStapleCache::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    INT iValue = 1;
    params.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);
    pOwner->m_pStapleCaches[m_byFieldId] = this;

    iValue = 0;
    params.FetchValueINT(_T("CacheStaples"), iValue);
    m_bCacheStaple = (0 != iValue);

    iValue = 0;
    params.FetchValueINT(_T("CachePlaquttes"), iValue);
    m_bCachePlaqutte = (0 != iValue);

    iValue = 1;
    params.FetchValueINT(_T("CacheFmunu"), iValue);
    m_bCacheFmunu = (0 != iValue);

    iValue = 0;
    params.FetchValueINT(_T("CacheRotationLink"), iValue);
    m_bCacheRotationLink = (0 != iValue);

    Real fValue = F(0.0);
    params.FetchValueReal(_T("CacheRotationLinkCharge"), fValue);
    m_fCacheRotationLinkCharge = fValue;

    iValue = 0;
    params.FetchValueINT(_T("CacheRotationLinkPhaseFieldId"), iValue);
    m_byCacheRotationLinkPhaseFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    params.FetchValueINT(_T("FermionUpdate"), iValue);
    m_bFermionUpdate = (0 != iValue);

    iValue = 0;
    params.FetchValueINT(_T("GaugeUpdate"), iValue);
    m_bGaugeUpdate = (0 != iValue);

    iValue = 1;
    params.FetchValueINT(_T("UseEffectiveGauge"), iValue);
    m_bUseEffectiveGauge = (0 != iValue);
}

template<class FieldType>
CStapleCacheT<FieldType>::~CStapleCacheT()
{
    //Improve-1 (I6): deregister every owned extent before it dies.
    m_StapleArraySet.Unbind();
    m_PlaqArraySet.Unbind();
    m_FmunuHandle.Unbind();
    m_RotationLinksHandle.Unbind();
    if (NULL != m_pDeviceStaplePtr)
    {
        checkCudaErrors(__cudaFree(m_pDeviceStaplePtr));
    }
    if (NULL != m_pDevicePlaqPtr)
    {
        checkCudaErrors(__cudaFree(m_pDevicePlaqPtr));
    }
    if (NULL != m_pDeviceFmunuPtr)
    {
        checkCudaErrors(__cudaFree(m_pDeviceFmunuPtr));
    }
    if (NULL != m_pDeviceGaugeRotationLinksPtr)
    {
        checkCudaErrors(__cudaFree(m_pDeviceGaugeRotationLinksPtr));
    }
    for (INT i = 0; i < m_pPooledFields.Num(); ++i)
    {
        m_pPooledFields[i]->Return();
    }
}

template<class FieldType>
void CStapleCacheT<FieldType>::InitialBuffers(BYTE byFieldId)
{
    if (NULL == m_pDeviceStaplePtr && m_bCacheStaple)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceStaplePtr, sizeof(typename FieldType::_Gauge*) * 6));
        typename FieldType::_Gauge* pointers[6];
        TArray<const CHaloBufferHandle*> memberHandles;
        for (INT i = 0; i < 6; ++i)
        {
            FieldType* pPooled = dynamic_cast<FieldType*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
            pointers[i] = (typename FieldType::_Gauge*)pPooled->GetData();
            m_pPooledFields.AddItem(pPooled);
            memberHandles.AddItem(pPooled->GetHaloBufferHandle());
        }
        checkCudaErrors(cudaMemcpy(m_pDeviceStaplePtr, pointers, sizeof(typename FieldType::_Gauge*) * 6, cudaMemcpyHostToDevice));
        //Improve-1 (I6, 3.4.1): the array's entries are the pooled fields above;
        //bind the set so launches receiving the array expand to them. The
        //membership is fixed until this cache is destroyed.
        m_StapleArraySet.Bind(reinterpret_cast<const BYTE*>(m_pDeviceStaplePtr),
            sizeof(typename FieldType::_Gauge*) * 6, memberHandles);
        appParanoiac(_T("CStapleCacheT::Initial staple buffer\n"));
    }
    if (NULL == m_pDevicePlaqPtr && m_bCachePlaqutte)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pDevicePlaqPtr, sizeof(typename FieldType::_Gauge*) * 6));
        typename FieldType::_Gauge* pointers[6];
        TArray<const CHaloBufferHandle*> memberHandles;
        for (INT i = 0; i < 6; ++i)
        {
            FieldType* pPooled = dynamic_cast<FieldType*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
            pointers[i] = (typename FieldType::_Gauge*)pPooled->GetData();
            m_pPooledFields.AddItem(pPooled);
            memberHandles.AddItem(pPooled->GetHaloBufferHandle());
        }
        checkCudaErrors(cudaMemcpy(m_pDevicePlaqPtr, pointers, sizeof(typename FieldType::_Gauge*) * 6, cudaMemcpyHostToDevice));
        m_PlaqArraySet.Bind(reinterpret_cast<const BYTE*>(m_pDevicePlaqPtr),
            sizeof(typename FieldType::_Gauge*) * 6, memberHandles);
        appParanoiac(_T("CStapleCacheT::Initial plaqutte buffer\n"));
    }
    if (NULL == m_pDeviceFmunuPtr && m_bCacheFmunu)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceFmunuPtr, sizeof(typename FieldType::_Gauge) * _HC_Volume * 6));
        //Improve-1 (I6, appendix A.5): consumed strictly site-locally
        //(_kernelDFermionWilsonSquareCloverSU3 reads [component * V + site]);
        //register a LocalOnly extent so writes are versioned, no halo traffic.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceFmunuPtr);
        sInfo.m_uiCapacityBytes = sizeof(typename FieldType::_Gauge) * _HC_Volume * 6;
        sInfo.m_uiBytesPerSite = sizeof(typename FieldType::_Gauge) * 6;
        sInfo.m_uiLocalSiteCount = _HC_Volume;
        sInfo.m_uiHaloSiteCount = 0;
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = byFieldId;
        sInfo.m_bHaloCapable = FALSE;
        m_FmunuHandle.Bind(sInfo);
        appParanoiac(_T("CStapleCacheT::Initial Fmunu buffer\n"));
    }
    if (NULL == m_pDeviceGaugeRotationLinksPtr && m_bCacheRotationLink)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceGaugeRotationLinksPtr, sizeof(typename FieldType::_Gauge) * _HC_Volume * 16));
        //Improve-1 (I6, appendix A.5): consumed strictly site-locally (the
        //_kernelDFermionKS_PR_Cached_* family reads [site * 8 + idx] inside two
        //8*V halves, never a neighbour slot); LocalOnly like Fmunu.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceGaugeRotationLinksPtr);
        sInfo.m_uiCapacityBytes = sizeof(typename FieldType::_Gauge) * _HC_Volume * 16;
        sInfo.m_uiBytesPerSite = sizeof(typename FieldType::_Gauge) * 16;
        sInfo.m_uiLocalSiteCount = _HC_Volume;
        sInfo.m_uiHaloSiteCount = 0;
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = byFieldId;
        sInfo.m_bHaloCapable = FALSE;
        m_RotationLinksHandle.Bind(sInfo);
        appParanoiac(_T("CStapleCacheT::Initial Rotation link buffer\n"));
    }
    _CHECKCUDA;
}

__CLGIMPLEMENT_CLASS(CStapleCacheSU3)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================