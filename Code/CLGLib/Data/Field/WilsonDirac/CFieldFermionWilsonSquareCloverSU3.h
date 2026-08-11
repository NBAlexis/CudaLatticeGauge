//=============================================================================
// FILENAME : CFieldFermionWilsonSquareCloverSU3.h
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [05/17/2025 nbale]
//=============================================================================
#include "CFieldFermionWilsonKernel.h"
#include "CFieldFermionWilsonSquareSU3.h"
#include "CFieldFermionWilsonSquareSU3D.h"
#include "CFieldFermionWilsonSquareSU3EM.h"
#include "CFieldFermionWilsonSquareSU3DR.h"
#include "Update/CStapleCache.h"
#include "GaugeSmearing/CGaugeSmearing.h"

#ifndef _CFIELDFERMIONWILSONSQUARECLOVERSU3_H_
#define _CFIELDFERMIONWILSONSQUARECLOVERSU3_H_

__BEGIN_NAMESPACE

inline class CStapleCache* appGetStapleCache(BYTE byFieldId);

//__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareCloverSU3)

/**
* Improve-1 (I6, multi-GPU-improve1.md 3.1): non-field managed extents still
* need a tag id for the halo exchange MPI messages (SHaloBufferInfo::m_byFieldId
* is a tag, never an identity -- see CHaloBufferHandle.h). Reserve 250+ for
* owner buffers so they cannot collide with real field ids.
*/
static const BYTE CLG_PseudoFieldIdCloverPhidPhi = 250;

class CLGAPI CFieldFermionWilsonSquareCloverBuffer : public CRegisteredBufferCache
{
public:

    CFieldFermionWilsonSquareCloverBuffer()
        : CRegisteredBufferCache()
        , m_pPhidPhiBuffer(NULL)
        , m_uiBufferSize(0)
    {

    }
    ~CFieldFermionWilsonSquareCloverBuffer()
    {
        if (NULL != m_pPhidPhiBuffer)
        {
            //Improve-1 (I6): deregister before the extent dies.
            m_HaloBuffer.Unbind();
            checkCudaErrors(__cudaFree(m_pPhidPhiBuffer));
            m_pPhidPhiBuffer = NULL;
        }
    }

    template<typename T>
    T* GetPhidPhiBuffer(ULONGLONG size)
    {
        if (NULL != m_pPhidPhiBuffer && m_uiBufferSize < size)
        {
            m_HaloBuffer.Unbind();
            checkCudaErrors(__cudaFree(m_pPhidPhiBuffer));
            m_pPhidPhiBuffer = NULL;
            m_uiBufferSize = 0;
        }
        if (NULL == m_pPhidPhiBuffer)
        {
            //Improve-1 (I6, appendix A.5): _kernelCloverForceCacheIndex reads
            //this buffer at 1-hop staple neighbours, so it is a managed
            //HaloCapable extent. The storage is site-major (the sigma index is
            //the fast index inside one site slot) so the generic per-site halo
            //machinery applies, and halo slots are appended past the requested
            //interior bytes. Single-GPU: halo site count is 0 and the
            //allocation is byte-identical to the historical one.
            const UINT uiVolume = _HC_Volume;
            const UINT uiBytesPerSite = (uiVolume > 0) ? static_cast<UINT>(size / uiVolume) : 0;
            const UINT uiHaloSites = static_cast<UINT>(_HC_HaloSiteCount());
            checkCudaErrors(__cudaMalloc((void**)&m_pPhidPhiBuffer, size + uiHaloSites * static_cast<ULONGLONG>(uiBytesPerSite)));
            m_uiBufferSize = size;

            SHaloBufferInfo sInfo;
            sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pPhidPhiBuffer);
            sInfo.m_uiCapacityBytes = static_cast<size_t>(size) + uiHaloSites * static_cast<size_t>(uiBytesPerSite);
            sInfo.m_uiBytesPerSite = uiBytesPerSite;
            sInfo.m_uiLocalSiteCount = uiVolume;
            sInfo.m_uiHaloSiteCount = uiHaloSites;
            sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
            sInfo.m_byFieldId = CLG_PseudoFieldIdCloverPhidPhi;
            sInfo.m_bHaloCapable = TRUE;
            m_HaloBuffer.Bind(sInfo);
        }
        return reinterpret_cast<T*>(m_pPhidPhiBuffer);
    }

    /** Improve-1 (I6): registry/diagnostic access for tests. */
    CHaloBufferHandle* GetHaloBufferHandle() { return &m_HaloBuffer; }

    void* m_pPhidPhiBuffer;
    ULONGLONG m_uiBufferSize;

    //Improve-1 (I6): handle over the whole extent (interior + halo tail).
    CHaloBufferHandle m_HaloBuffer;

    static CFieldFermionWilsonSquareCloverBuffer* GetInstance();
    static CFieldFermionWilsonSquareCloverBuffer* m_pPointer;
};

template<class CFieldWilson>
class __DLL_EXPORT CFieldFermionWilsonSquareClover : public CFieldWilson
{
    //__CLGDECLARE_FIELD(CFieldFermionWilsonSquareCloverSU3)

public:

    CFieldFermionWilsonSquareClover()
        : CFieldWilson()
        , m_fCsw(1.0)
        //, m_pPhidPhi(NULL)
    {
        //6 is for six different sigma types
        //checkCudaErrors(__cudaMalloc((void**)& m_pPhidPhi, this->m_uiSiteCount * 6 * sizeof(deviceWilsonVectorSU3)));
    }

    ~CFieldFermionWilsonSquareClover()
    {
        //if (NULL != m_pPhidPhi)
        //{
        //    checkCudaErrors(__cudaFree(m_pPhidPhi));
        //    m_pPhidPhi = NULL;
        //}
    }

    void CopyParamTo(CField* U) const override
    {
        if (NULL == U || EFT_FermionWilsonSquareSU3 != U->GetFieldType())
        {
            appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
            return;
        }

        CFieldWilson::CopyParamTo(U);

        CFieldFermionWilsonSquareClover<CFieldWilson>* pField = dynamic_cast<CFieldFermionWilsonSquareClover<CFieldWilson>*>(U);
        if (NULL != pField)
        {
            pField->m_fCsw = m_fCsw;

            //m_pPhidPhi is updated in every force calculation, no need to copy it.
            //CCommonKernel<deviceSU3>::CopyBuffer(pField->m_pPhidPhi, m_pPhidPhi, this->m_uiSiteCount * 6);
        }
    }

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldWilson::InitialOtherParameters(params);
        params.FetchValueDOUBLE(_T("Csw"), m_fCsw);
    }

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        _RECORD(CFieldFermionWilsonSquareClover::DOperator);
        CFieldWilson::DOperator(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, bDagger, eOCT, fRealCoeff, cCmpCoeff);

        deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
        const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
        const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
        const deviceSU3* pFmunu = (NULL == appGetStapleCache(byGaugeFieldId)) ? NULL : ((const deviceSU3*)(appGetStapleCache(byGaugeFieldId)->GetFmunu()));
        if (NULL == pFmunu)
        {
            appCrucial(_T("fmunu is not calculated!\n"));
            _FAIL_EXIT;
        }

        CFieldFermionWilsonKernel::DOperatorClover(pTarget, pSource, pGauge, pFmunu, this->m_byFieldId, byGaugeFieldId,
            m_fCsw * this->m_fKai, bDagger, eOCT, fRealCoeff, cCmpCoeff);

        //preparethreadE(6);
        //_LAUNCH_KERNEL(_kernelDFermionWilsonSquareCloverSU3, block, threads,
        //    pSource,
        //    pGauge,
        //    pFmunu,
        //    pTarget,
        //    m_fKai * m_fCsw,
        //    m_byFieldId,
        //    bDagger,
        //    eOCT,
        //    fRealCoeff,
        //    cCmpCoeff);
    }

    void DerivateDOperator(DOUBLE fCoeff, void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        _RECORD(CFieldFermionWilsonSquareClover::DerivateDOperator);

        CFieldWilson::DerivateDOperator(fCoeff, pForce, pDphi, pDDphi, pGaugeBuffer, byGaugeFieldId);

        deviceSU3* pForceSU3 = (deviceSU3*)pForce;
        const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
        const deviceWilsonVectorSU3* phid = (deviceWilsonVectorSU3*)pDphi;
        const deviceWilsonVectorSU3* phi = (deviceWilsonVectorSU3*)pDDphi;

        //if (NULL == m_pPhidPhi)
        //{
        //    checkCudaErrors(__cudaMalloc((void**)&m_pPhidPhi, this->m_uiSiteCount * 6 * sizeof(deviceWilsonVectorSU3)));
        //}
        deviceSU3* pPhidPhi = CFieldFermionWilsonSquareCloverBuffer::GetInstance()->GetPhidPhiBuffer<deviceSU3>(this->m_uiSiteCount * 6 * sizeof(deviceSU3));
        CFieldFermionWilsonKernel::PrepareSigmaMunu(phi, phid, pPhidPhi);
        CFieldFermionWilsonKernel::CloverForce(pPhidPhi, pGauge, pForceSU3, this->m_byFieldId, byGaugeFieldId, fCoeff * m_fCsw);
    }

public:

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldWilson::GetInfos(tab);
        sRet = sRet + tab + _T("Csw : ") + appToString(m_fCsw) + _T("\n");
        return sRet;
    }

    DOUBLE m_fCsw;
    //deviceSU3* m_pPhidPhi;
};


__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareCloverSU3)
class CLGAPI CFieldFermionWilsonSquareCloverSU3 : public CFieldFermionWilsonSquareClover<CFieldFermionWilsonSquareSU3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionWilsonSquareCloverSU3)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareCloverEMSU3)
class CLGAPI CFieldFermionWilsonSquareCloverEMSU3 : public CFieldFermionWilsonSquareClover<CFieldFermionWilsonSquareSU3EM>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionWilsonSquareCloverEMSU3)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareCloverSU3D)
class CLGAPI CFieldFermionWilsonSquareCloverSU3D : public CFieldFermionWilsonSquareClover<CFieldFermionWilsonSquareSU3D>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionWilsonSquareCloverSU3D)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareCloverSU3DR)
class CLGAPI CFieldFermionWilsonSquareCloverSU3DR : public CFieldFermionWilsonSquareClover<CFieldFermionWilsonSquareSU3DR>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionWilsonSquareCloverSU3DR)
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARECLOVERSU3_H_

//=============================================================================
// END OF FILE
//=============================================================================