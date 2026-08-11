//=============================================================================
// FILENAME : CFieldGaugeLink.h
// 
// DESCRIPTION:
// This is the common class for all gauge fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/04/2018 nbale]
//=============================================================================
#pragma once

#include "CFieldGaugeKernel.h"

#ifndef _CFIELDGAUGE_LINK_H_
#define _CFIELDGAUGE_LINK_H_

#define gaugeKernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeKernelFuncionEnd \
    } 

#define __DEFINE_GAUGE_LINK(CLASSNAME, DEVICEDATA, N, FIELDTYPE) \
__CLG_REGISTER_HELPER_HEADER(CLASSNAME) \
class CLGAPI CLASSNAME : public CFieldGaugeLink<DEVICEDATA, N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CLASSNAME) \
public: \
    EFieldType GetFieldType() const override { return FIELDTYPE; } \
};


__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeLink : public CFieldGauge
{

public:

    typedef deviceGauge _Gauge;

    typedef CCommonKernelField<deviceGauge> _FieldKernel;
    typedef CCommonKernelLink<deviceGauge> _LinkKernel;
    typedef CFieldGaugeKernel<deviceGauge, matrixN> _GaugeKernel;

    CFieldGaugeLink() : CFieldGauge()
    {
        //Multi-GPU (Phase 1): append halo link slots after the local links so a
        //split-direction out-of-lattice neighbour can be read from the same flat
        //buffer (Design B, Docs/MultiGPU-Plan.md section 10). m_uiLinkeCount stays
        //the physics volume so CopyDataOut / gather-scatter are unchanged; only the
        //allocation grows. Single-GPU: _HC_HaloLinkCount() == 0 -> identical size.
        m_uiHaloLinkCount = _HC_HaloLinkCount();
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData,
            (m_uiLinkeCount + m_uiHaloLinkCount) * sizeof(deviceGauge)));

        //Improve-1 (multi-GPU-improve1.md 3.1/3.2): bind the halo handle to
        //this exact extent. Link field: Dir matrices per site; site counts are
        //lattice sites (halo slots are per-site too). The field id is assigned
        //after construction and synced through SetFieldId.
        SHaloBufferInfo sInfo;
        sInfo.m_pDeviceData = reinterpret_cast<BYTE*>(m_pDeviceData);
        sInfo.m_uiCapacityBytes = (m_uiLinkeCount + m_uiHaloLinkCount) * sizeof(deviceGauge);
        sInfo.m_uiBytesPerSite = _HC_Dir * sizeof(deviceGauge);
        sInfo.m_uiLocalSiteCount = _HC_Volume;
        sInfo.m_uiHaloSiteCount = _HC_HaloSiteCount();
        sInfo.m_ullLayoutGeneration = appGetLayoutGeneration();
        sInfo.m_byFieldId = 0;
        sInfo.m_bHaloCapable = TRUE;
        m_HaloBuffer.Bind(sInfo);
    }

    ~CFieldGaugeLink()
    {
        checkCudaErrors(__cudaFree(m_pDeviceData));
        m_pDeviceData = NULL;
    }

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override
    {
        if (!CFileSystem::IsFileExist(sFileName))
        {
            appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
            _FAIL_EXIT;
        }

        switch (eFileType)
        {
        case EFFT_CLGBin:
#if _CLG_DOUBLEFLOAT
        case EFFT_CLGBinDouble:
#else
        case EFFT_CLGBinFloat:
#endif
        {
            UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            const UINT uiPerRankSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount);
#if _CLG_MULTI_GPU
            //Multi-GPU: on disk the file is the whole global lattice in global-site
            //order (written by the gather in SaveToFile). Every rank read the full
            //file above; scatter it so each rank keeps only its own sub-lattice in
            //local-site order. Inverse of GatherFieldToRoot (Docs/MultiGPU-Plan.md 8.5).
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                const UINT uiExpectGlobal = uiPerRankSize * appGetComm()->Size();
                if (uiSize != uiExpectGlobal)
                {
                    appCrucial(_T("Loading file size not match (MG): %s, %d, expecting global %d"), sFileName.c_str(), uiSize, static_cast<INT>(uiExpectGlobal));
                }
                const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_Dir);
                BYTE* byLocal = (BYTE*)malloc(uiPerRankSize);
                appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
                free(data);
                data = byLocal;
                uiSize = uiPerRankSize;
            }
            else
#endif
            if (uiSize != uiPerRankSize)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(uiPerRankSize));
            }
            InitialWithByte(data);
            free(data);
            FixBoundary(EFB_Field);
        }
        break;
#if _CLG_DOUBLEFLOAT
        case EFFT_CLGBinFloat:
        {
            //P5-1.1: on disk the file is the whole global lattice in global-site
            //order; read it fully, convert float -> Real, then scatter to the
            //per-rank sub-lattice on multi-GPU (previously each rank read by its
            //LOCAL link count -> misaligned/incomplete data, silent wrong).
            UINT uiSize = static_cast<UINT>(sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            const UINT uiPerRankElems = 2 * MatrixN() * MatrixN() * _HC_LinkCount;
            //Declared outside the #if so the non-MG build (single rank, the file
            //IS the whole lattice) still has uiGlobalElems in scope (P5-1.1 fix).
            UINT uiGlobalElems = uiPerRankElems;
#if _CLG_MULTI_GPU
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                uiGlobalElems = uiPerRankElems * appGetComm()->Size();
                if (uiSize != sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount * appGetComm()->Size())
                {
                    appCrucial(_T("Loading file size not match (MG): %s, %d, expecting global %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount * appGetComm()->Size()));
                }
            }
            else
#endif
            if (uiSize != sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<UINT>(sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
            }
            //P5-1.1: allocate at the GLOBAL element count (the conversion loop
            //writes every element of the full file); scatter shrinks to the
            //per-rank share afterwards. (Allocating local here overflowed the
            //heap on -n2: global elements written into a local-sized buffer.)
            BYTE* data = (BYTE*)malloc(sizeof(Real) * uiGlobalElems);
            Real* rdata = (Real*)data;
            for (UINT i = 0; i < uiGlobalElems; ++i)
            {
                rdata[i] = static_cast<Real>(fdata[i]);
            }
            free(fdata);
#if _CLG_MULTI_GPU
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_Dir);
                BYTE* byLocal = (BYTE*)malloc(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount);
                appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
                free(data);
                data = byLocal;
            }
#endif
            InitialWithByte(data);
            free(data);
            ElementNormalize();
            FixBoundary(EFB_Field);
        }
        break;
#else
        case EFFT_CLGBinDouble:
        {
            //P5-1.1: same as the float-load branch: full-file read, convert
            //double -> Real, scatter on multi-GPU.
            UINT uiSize = static_cast<UINT>(sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
            DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
            const UINT uiPerRankElems = 2 * MatrixN() * MatrixN() * _HC_LinkCount;
            //Declared outside the #if so the non-MG build (single rank, the file
            //IS the whole lattice) still has uiGlobalElems in scope (P5-1.1 fix).
            UINT uiGlobalElems = uiPerRankElems;
#if _CLG_MULTI_GPU
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                uiGlobalElems = uiPerRankElems * appGetComm()->Size();
                if (uiSize != sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount * appGetComm()->Size())
                {
                    appCrucial(_T("Loading file size not match (MG): %s, %d, expecting global %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount * appGetComm()->Size()));
                }
            }
            else
#endif
            if (uiSize != sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
            {
                appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
            }
            //P5-1.1: allocate at the GLOBAL element count (see float branch).
            BYTE* data = (BYTE*)malloc(sizeof(Real) * uiGlobalElems);
            Real* rdata = (Real*)data;
            for (UINT i = 0; i < uiGlobalElems; ++i)
            {
                rdata[i] = static_cast<Real>(ddata[i]);
            }
            free(ddata);
#if _CLG_MULTI_GPU
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                const UINT uiBytesPerSite = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_Dir);
                BYTE* byLocal = (BYTE*)malloc(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount);
                appGetComm()->ScatterFieldFromRoot(data, uiBytesPerSite, byLocal);
                free(data);
                data = byLocal;
            }
#endif
            InitialWithByte(data);
            free(data);
            FixBoundary(EFB_Field);
        }
        break;
#endif
        case EFFT_CLGBinCompressed:
        {
            InitialWithByteCompressed(sFileName);
            FixBoundary(EFB_Field);
        }
        break;
        case EFFT_BridgePPTXT:
        case EFFT_BridgePPBin:
        {
            //Bridge++ (ILDG-style) configuration: every SU(3) link is 18 reals
            //(Re,Im row-major over 3x3), links in site-major then direction
            //order, sites x-fastest -- the same order as _deviceSiteIndexToInt4
            //and _deviceGetLinkIndex(site,dir), so the element stream maps 1:1
            //onto our per-link buffer. Single-GPU for now (external per-test
            //input; a rank-0 gather of arbitrary external formats is out of
            //scope, tests use -n1 for Bridge++ files).
            const UINT uiLinkElems = 2 * MatrixN() * MatrixN();
            UINT uiTotal = uiLinkElems * _HC_LinkCount;
#if _CLG_MULTI_GPU
            //The Bridge++ file is the whole GLOBAL lattice (external input);
            //read it fully on every rank, reorder, then scatter per rank.
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                uiTotal = uiLinkElems * _HC_LinkCount * appGetComm()->Size();
            }
#endif
            Real* rdata = (Real*)malloc(sizeof(Real) * uiTotal);
            if (EFFT_BridgePPTXT == eFileType)
            {
#if _CLG_WIN
                FILE* fp = NULL;
                fopen_s(&fp, sFileName.c_str(), "r");
#else
                FILE* fp = fopen(sFileName.c_str(), "r");
#endif
                if (NULL == fp)
                {
                    free(rdata);
                    appCrucial(_T("Cannot open Bridge++ text %s\n"), sFileName.c_str());
                    _FAIL_EXIT;
                }
                for (UINT i = 0; i < uiTotal; ++i)
                {
                    DOUBLE d = 0.0;
#if _CLG_WIN
                    if (1 != fscanf_s(fp, "%lf", &d))
#else
                    if (1 != fscanf(fp, "%lf", &d))
#endif
                    {
                        fclose(fp);
                        free(rdata);
                        appCrucial(_T("Bridge++ text read failed at element %d of %d\n"), i, uiTotal);
                        _FAIL_EXIT;
                    }
                    rdata[i] = static_cast<Real>(d);
                }
                fclose(fp);
            }
            else
            {
                const UINT uiBytes = sizeof(DOUBLE) * uiTotal;
#if _CLG_WIN
                FILE* fp = NULL;
                fopen_s(&fp, sFileName.c_str(), "rb");
#else
                FILE* fp = fopen(sFileName.c_str(), "rb");
#endif
                if (NULL == fp)
                {
                    free(rdata);
                    appCrucial(_T("Cannot open Bridge++ bin %s\n"), sFileName.c_str());
                    _FAIL_EXIT;
                }
                DOUBLE* ddata = (DOUBLE*)malloc(uiBytes);
                if (uiBytes != fread(ddata, 1, uiBytes, fp))
                {
                    fclose(fp);
                    free(ddata);
                    free(rdata);
                    appCrucial(_T("Bridge++ bin read size mismatch %s\n"), sFileName.c_str());
                    _FAIL_EXIT;
                }
                fclose(fp);
                //Bridge++ binary configs are big-endian (network order); the
                //host is little-endian here, so swap each 8-byte double.
                {
                    BYTE* pb = (BYTE*)ddata;
                    for (UINT i = 0; i < uiTotal; ++i)
                    {
                        BYTE* b = pb + i * 8;
                        BYTE tmp = b[0]; b[0] = b[7]; b[7] = tmp;
                        tmp = b[1]; b[1] = b[6]; b[6] = tmp;
                        tmp = b[2]; b[2] = b[5]; b[5] = tmp;
                        tmp = b[3]; b[3] = b[4]; b[4] = tmp;
                    }
                }
                for (UINT i = 0; i < uiTotal; ++i)
                {
                    rdata[i] = static_cast<Real>(ddata[i]);
                }
                free(ddata);
            }
            //Bridge++ orders sites x-fastest (ILDG); CLGLib orders sites
            //t-fastest (_deviceSiteIndexToInt4Baking). Reorder element blocks:
            //block = one link (18 reals), 4 blocks per site, dirs x,y,z,t both.
            if (4 == _HC_Dir)
            {
                //The reorder runs on the WHOLE global lattice (Bridge++ file is
                //global; on -n2 the buffer was read at the global size), so use
                //GLOBAL lengths here; the per-rank scatter afterwards handles
                //the decomposition. Identity when unsplit.
                const UINT uiLX = static_cast<UINT>(_HC_GlobalLx), uiLY = static_cast<UINT>(_HC_GlobalLy);
                const UINT uiLZ = static_cast<UINT>(_HC_GlobalLz), uiLT = static_cast<UINT>(_HC_GlobalLt);
                Real* rdata2 = (Real*)malloc(sizeof(Real) * uiTotal);
                const UINT uiLinkElems = 2 * MatrixN() * MatrixN();
                for (UINT it = 0; it < uiLT; ++it)
                {
                    for (UINT iz = 0; iz < uiLZ; ++iz)
                    {
                        for (UINT iy = 0; iy < uiLY; ++iy)
                        {
                            for (UINT ix = 0; ix < uiLX; ++ix)
                            {
                                const UINT uiBridge = ((it * uiLZ + iz) * uiLY + iy) * uiLX + ix;
                                const UINT uiMine = ((ix * uiLY + iy) * uiLZ + iz) * uiLT + it;
                                memcpy(rdata2 + uiMine * _HC_Dir * uiLinkElems,
                                    rdata + uiBridge * _HC_Dir * uiLinkElems,
                                    sizeof(Real) * _HC_Dir * uiLinkElems);
                            }
                        }
                    }
                }
                free(rdata);
                rdata = rdata2;
            }
#if _CLG_MULTI_GPU
            //Scatter the reordered (global t-fastest) buffer to the per-rank
            //sub-lattice, mirroring the CLGBin branches (no-op on a lone rank).
            if (NULL != appGetComm() && appGetComm()->Size() > 1)
            {
                const UINT uiPerRank = sizeof(Real) * uiLinkElems * _HC_LinkCount;
                const UINT uiBytesPerSite = sizeof(Real) * uiLinkElems * _HC_Dir;
                BYTE* byLocal = (BYTE*)malloc(uiPerRank);
                appGetComm()->ScatterFieldFromRoot((BYTE*)rdata, uiBytesPerSite, byLocal);
                free(rdata);
                rdata = (Real*)byLocal;
            }
#endif
            InitialWithByte((BYTE*)rdata);
            free(rdata);
            FixBoundary(EFB_Field);
        }
        break;
        default:
            appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eFileType).c_str());
            break;

        }
        NotifyWritten();
    }

    void InitialWithByte(BYTE* byData) override
    {
        CCommonKernelField<deviceGauge>::InitialWithByte(m_pDeviceData, m_uiLinkeCount, byData);
#if _CLG_MULTI_GPU
        //Full interior overwrite (file load / scatter): any previously baked halo is
        //now stale. Without this, an Ensure() that ran earlier (e.g. while the buffer
        //was still zero-init before load) leaves validWidth>=1, so the next stencil's
        //Ensure() early-returns and reads a zero/stale halo -> boundary plaquettes read
        //zero links and the plaquette energy is wrong (~one t-face * N per rank).
        NotifyWritten();
#endif
    }

    void InitialWithByteCompressed(const CCString& fileName) override { appCrucial(_T("CFieldGaugeLink: InitialWithByteCompressed not supoorted!\n")); }

    void InitialField(EFieldInitialType eInitialType) override
    {
        CCommonKernelLink<deviceGauge>::InitialBuffer(m_pDeviceData, m_byFieldId, eInitialType);
        if ((EFIT_RandomGenerator == eInitialType || EFIT_RandomGaussian == eInitialType) && abs(_HC_GaugeMomentumFactor - F(1.0)) > _CLG_FLT_EPSILON)
        {
            CCommonKernelField<deviceGauge>::ScalarMultply(m_pDeviceData, m_uiLinkeCount, sqrt(_HC_GaugeMomentumFactor));
        }
#if _CLG_MULTI_GPU
        //Full interior (re)initialisation: previously baked halo is stale.
        NotifyWritten();
#endif
    }

    void DebugPrintMe() const override
    {
        CCommonKernelLink<deviceGauge>::DebugPrint(m_pDeviceData, m_uiLinkeCount);
    }

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override
    {
        if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        if (NULL != pStaple && GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
            return;
        }

        _RECORD(CFieldGaugeLink::CalculateForceAndStaple);

        CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
        CFieldGaugeLink<deviceGauge, matrixN>* pStapleSU3 = NULL == pStaple ? NULL : dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple(
            m_pDeviceData,
            m_byFieldId,
            pForceSU3->m_pDeviceData,
            NULL == pStapleSU3 ? NULL : pStapleSU3->m_pDeviceData,
            betaOverN);
    }

    void CalculateForceAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const override
    {
        if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
            return;
        }
        _RECORD(CFieldGaugeLink::CalculateForceAnisotropy);
        CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAnisotropy(
            m_pDeviceData,
            m_byFieldId,
            pForceSU3->m_pDeviceData,
            betaOverN,
            xi);
    }

    void CalculateOnlyStaple(CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: staple field is not SU3");
            return;
        }
        CFieldGaugeLink<deviceGauge, matrixN>* pStapleSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple(m_pDeviceData, m_byFieldId, pStapleSU3->m_pDeviceData);
    }

    virtual void CalculateAllStaples(deviceGauge* const * ppDeviceStaple) const
    {
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateAllStaples(m_pDeviceData, m_byFieldId, ppDeviceStaple);
    }

    virtual void CalculateAllPlaquttes(deviceGauge* const* ppDevicePlaq) const
    {
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateAllPlaquttes(m_pDeviceData, m_byFieldId, ppDevicePlaq);
    }

    virtual void CalculateFmunu(deviceGauge* pDeviceFmunu) const
    {
        CFieldGaugeKernel<deviceGauge, matrixN>::CalculateFmunu(m_pDeviceData, m_byFieldId, pDeviceFmunu);
    }

    virtual void CacheRotationKSBuffer(deviceGauge* buffer, UBOOL bHasU1, Real fCharge, const Real* pPhase) const
    {
        CFieldGaugeKernel<deviceGauge, matrixN>::CacheKSRotationGaugeBuffer(m_pDeviceData, m_byFieldId, buffer, pPhase, fCharge, bHasU1);
    }

    void MakeRandomGenerator() override
    {
        _RECORD(CFieldGaugeLink::MakeRandomGenerator);
        InitialField(EFIT_RandomGenerator);
    }

    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergy);
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy(m_pDeviceData, m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergyAnisotropy);
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyAnisotropy(m_pDeviceData, m_byFieldId, betaOverN, xi);
    }

    DOUBLE CalculatePlaqutteEnergyOriginal(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergyOriginal);
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy(m_pDeviceData, m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergyUseClover);
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseClover(m_pDeviceData, m_byFieldId, betaOverN);
    }

    DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(DOUBLE betaOverN, DOUBLE xi) const override
    {
        _RECORD(CFieldGaugeLink::CalculatePlaqutteEnergyUseCloverAnisotropy);
        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseCloverAnisotropy(m_pDeviceData, m_byFieldId, betaOverN, xi);
    }

    DOUBLE CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        if (NULL == pStaple || GetFieldType() != pStaple->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
            return F(0.0);
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pStapleSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pStaple);

        return CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUsingStaple(m_pDeviceData, m_byFieldId, betaOverN, pStapleSU3->m_pDeviceData);
    }

    /**
    * See 10.1007/978-3-642-01850-3
    * It should be 1/2 <P, P>
    * where see Eq.~(A8) <X, Y> = -2 Tr[XY]
    * so it is -1.0 * tr[P^2] and because P is anti-hermitian, it is 1.0 * tr[P^dagger P]
    * LengthSq is tr[P^dagger P]
    */
    DOUBLE CalculateKinematicEnergy() const override
    {
        _RECORD(CFieldGaugeLink::CalculateKinematicEnergy);
        if (abs(_HC_GaugeMomentumFactor - F(1.0)) > _CLG_FLT_EPSILON)
        {
            return CCommonKernelField<deviceGauge>::LengthSq(m_pDeviceData, m_uiLinkeCount) / _HC_GaugeMomentumFactor;
        }
        return CCommonKernelField<deviceGauge>::LengthSq(m_pDeviceData, m_uiLinkeCount);
    }

#pragma endregion

#pragma region BLAS

    void Dagger() override
    {
        CCommonKernelField<deviceGauge>::Dagger(m_pDeviceData, m_uiLinkeCount);
        NotifyWritten();
    }

    void AxpyPlus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::AxpyPlus(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData);
        NotifyWritten();
    }

    void AxpyMinus(const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::AxpyMinus(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(Real a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Axpy(m_pDeviceData, m_uiLinkeCount, a, pSU3x->m_pDeviceData);
        NotifyWritten();
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Axpy(m_pDeviceData, m_uiLinkeCount, a, pSU3x->m_pDeviceData);
        NotifyWritten();
    }

    void Mul(const CField* x, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) override
    {
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::Mul(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void ApplyPhaseC(const CField* other) override;

    void ApplyStaggeredPhase(UINT iType = 0) override
    {
        if (1U == iType)
        {
            CCommonKernelLink<deviceGauge>::ApplyStaggeredPhaseMILC(m_pDeviceData, m_byFieldId);
            NotifyWritten();
            return;
        }

        CCommonKernelLink<deviceGauge>::ApplyStaggeredPhase(m_pDeviceData, m_byFieldId);
        NotifyWritten();
    }

    //Implemented in CFieldGaugeU1Real
    //Can't I just put it into another .h file when using GCC? How wierd...
    //chatGPT says if I put it into another .h file, I need to also explicit instantiation this stupid function alone
    void ApplyPhaseR(const CField* other, Real fCharge) override
    {
        //const class CFieldGaugeU1Real* phaser = dynamic_cast<const class CFieldGaugeU1Real*>(other);
        //if (NULL != phaser)
        if (EFT_GaugeReal == other->GetFieldType())
        {
            CCommonKernelField<deviceGauge>::ApplyPhaseR(m_pDeviceData, m_uiLinkeCount, (const Real*)other->GetData(), fCharge);
            NotifyWritten();
            return;
        }

        appCrucial(_T("ApplyPhase with phase type not supported!\n"));
    }

    void LeftMul(const CField* x, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) override
    {
        _RECORD(CFieldGaugeLink::LeftMul);
        if (NULL == x || GetFieldType() != x->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
            return;
        }
        const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
        CCommonKernelField<deviceGauge>::LeftMul(m_pDeviceData, m_uiLinkeCount, pSU3x->m_pDeviceData, bDaggerLeft, bDaggerRight);
        NotifyWritten();
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        CCommonKernelField<deviceGauge>::ScalarMultply(m_pDeviceData, m_uiLinkeCount, a);
        NotifyWritten();
    }

    void ScalarMultply(Real a) override
    {
        CCommonKernelField<deviceGauge>::ScalarMultply(m_pDeviceData, m_uiLinkeCount, a);
        NotifyWritten();
    }

    void SetOneDirectionUnity(BYTE byDir) override
    {
        if (0 == (byDir & 15))
        {
            return;
        }
        CCommonKernelLink<deviceGauge>::SetOneDirectionUnity(m_pDeviceData, m_byFieldId, byDir);
        NotifyWritten();
    }

    void SetOneDirectionZero(BYTE byDir) override
    {
        if (0 == (byDir & 15))
        {
            return;
        }
        CCommonKernelLink<deviceGauge>::SetOneDirectionZero(m_pDeviceData, m_byFieldId, byDir);
        NotifyWritten();
    }

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override
    {
        if (0 != _HC_ALog)
        {
            CCommonKernelLink<deviceGauge>::StrictLog(m_pDeviceData, m_byFieldId);
        }
        else
        {
            CCommonKernelLink<deviceGauge>::QuickLog(m_pDeviceData, m_byFieldId);
        }
        NotifyWritten();
    }

    void TA() override
    {
        CCommonKernelLink<deviceGauge>::QuickLog(m_pDeviceData, m_byFieldId);
        NotifyWritten();
    }

    /**
     * U=exp(iA)
     */
    void TransformToU() override
    {
        if (0 != _HC_ALog)
        {
            CCommonKernelLink<deviceGauge>::StrictExp(m_pDeviceData, m_byFieldId);
        }
        else
        {
            CCommonKernelLink<deviceGauge>::QuickExp(m_pDeviceData, m_byFieldId);
        }
        NotifyWritten();
    }

    void CalculateE_Using_U(CFieldGauge* pResoult) const override
    {
        if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);
        CCommonKernelLink<deviceGauge>::CalculateE_Using_U(m_pDeviceData, m_byFieldId, pUField->m_pDeviceData);
        pResoult->NotifyWritten();
    }

    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const override
    {
        if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }

        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);
        CCommonKernelLink<deviceGauge>::CalculateNablaE_Using_U(m_pDeviceData, m_byFieldId, pUField->m_pDeviceData, bNaive);
        pResoult->NotifyWritten();
    }

#pragma endregion

    void ExpMult(Real a, CField* U) const override
    {
        if (NULL == U || GetFieldType() != U->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
            return;
        }
        _RECORD(CFieldGaugeLink::ExpMult);
        CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(U);
        CCommonKernelLink<deviceGauge>::ExpMul(pUField->m_pDeviceData, m_byFieldId, m_pDeviceData, a);
#if _CLG_MULTI_GPU
        //Leapfrog gauge update wrote U in place: its halo is now stale.
        U->NotifyWritten();
#endif
    }

    void ElementNormalize() override
    {
        _RECORD(CFieldGaugeLink::ElementNormalize);
        CCommonKernelField<deviceGauge>::Norm(m_pDeviceData, m_uiLinkeCount);
#if _CLG_MULTI_GPU
        //Reunitarisation rewrote every local link: halo is stale.
        NotifyWritten();
#endif
    }

    cuDoubleComplex Dot(const CField* other) const override
    {
        if (NULL == other || GetFieldType() != other->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SUN");
            return make_cuDoubleComplex(0, 0);
        }

        const CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(other);
        return CCommonKernelField<deviceGauge>::Dot(m_pDeviceData, m_uiLinkeCount, pUField->m_pDeviceData);
    }

    DOUBLE GetLength() const override
    {
        return CCommonKernelField<deviceGauge>::LengthSq(m_pDeviceData, m_uiLinkeCount);
    }

    CCString SaveToCompressedFile(const CCString& fileName) const override { appCrucial(_T("CFieldGaugeLink: SaveToCompressedFile not supoorted!\n")); return _T(""); };

    BYTE* CopyDataOut(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOut(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    BYTE* CopyDataOutFloat(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOutFloat(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    BYTE* CopyDataOutDouble(UINT& uiSize) const override
    {
        return CCommonKernelField<deviceGauge>::CopyDataOutDouble(m_pDeviceData, m_uiLinkeCount, uiSize);
    }

    void CopyBufferTo(CField* pTarget) const override
    {
        if (NULL == pTarget || GetFieldType() != pTarget->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: target field is not SUN");
            return;
        }

        //CFieldGauge::CopyBufferTo(pTarget);

        CFieldGaugeLink<deviceGauge, matrixN>* pTargetField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pTarget);
        CCommonKernel<deviceGauge>::CopyBuffer(pTargetField->m_pDeviceData, m_pDeviceData, m_uiLinkeCount);
        pTarget->NotifyWritten();
    }

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer, BYTE byDir = 3) const override
    {
        CCommonKernelLink<deviceGauge>::PolyakovOnSpatialSite(m_pDeviceData, m_byFieldId, buffer, byDir);
    }

    UINT MatrixN() const override { return matrixN; }

    void SetAsConnection(const CField* n, const CField* n_p_m) override
    {
        _RECORD(CFieldGaugeLink::SetAsConnectionField);
        n->Connection(n_p_m, m_pDeviceData);
    }

    void SetAsConnection(const CField* n, Real fCoeff) override
    {
        _RECORD(CFieldGaugeLink::SetAsConnection);
        n->ConnectionSelf(m_pDeviceData, fCoeff);
    }

    void AddConnection(const CField* n, Real fCoeff) override
    {
        _RECORD(CFieldGaugeLink::AddConnection);
        n->AddConnectionSelf(m_pDeviceData, fCoeff);
    }

    void AddConnection(const CField* n, const CField* n_p_m, Real fCoeff) override
    {
        _RECORD(CFieldGaugeLink::AddConnection);
        n->AddConnectionTwoField(m_pDeviceData, n_p_m, fCoeff);
    }

    void AddLinkTo(CField* target, const SCHAR* devicePath, BYTE byPathLen, BYTE mu, Real fCoeff) const override
    {
        if (NULL == target || GetFieldType() != target->GetFieldType())
        {
            appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: AddLinkTo is not SU3");
            return;
        }
        CCommonKernelLink<deviceGauge>::AddLink(m_pDeviceData, (deviceGauge*)target->GetData(), fCoeff, devicePath, byPathLen, mu, m_byFieldId);
    }

    void AddNaikForce(const CFieldGauge* naikf0) override
    {
        _GaugeKernel::NaikForce(naikf0, this);
    }

    deviceGauge* m_pDeviceData;

    //Improve-1: this field's own buffer handle (pool copies each bind theirs).
    CHaloBufferHandle m_HaloBuffer;
    CHaloBufferHandle* GetHaloBufferHandle() override { return &m_HaloBuffer; }
    const CHaloBufferHandle* GetHaloBufferHandle() const override { return &m_HaloBuffer; }

    _GetData

};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1)

class CLGAPI CFieldGaugeU1 : public CFieldGaugeLink<CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeU1)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU2)

class CLGAPI CFieldGaugeSU2 : public CFieldGaugeLink<deviceSU2, 2>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU2)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3)

class CLGAPI CFieldGaugeSU3 : public CFieldGaugeLink<deviceSU3, 3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU3)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }

    void CopyBufferTo(CField* pTarget) const override;
    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToFile(const CCString& fileName, EFieldFileType eType = EFFT_CLGBin) const override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
};

#if _CLG_SU4_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeSU4, deviceSU4, 4, EFT_GaugeSU4)
#endif
#if _CLG_SU5_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeSU5, deviceSU5, 5, EFT_GaugeSU5)
#endif
#if _CLG_SU6_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeSU6, deviceSU6, 6, EFT_GaugeSU6)
#endif
#if _CLG_SU7_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeSU7, deviceSU7, 7, EFT_GaugeSU7)
#endif
#if _CLG_SU8_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeSU8, deviceSU8, 8, EFT_GaugeSU8)
#endif

#if _CLG_Z2_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeZ2, deviceZN<2>, 1, EFT_GaugeZ2)
#endif
#if _CLG_Z3_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeZ3, deviceZN<3>, 1, EFT_GaugeZ3)
#endif
#if _CLG_Z4_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeZ4, deviceZN<4>, 1, EFT_GaugeZ4)
#endif
#if _CLG_Z5_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeZ5, deviceZN<5>, 1, EFT_GaugeZ5)
#endif
#if _CLG_Z6_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeZ6, deviceZN<6>, 1, EFT_GaugeZ6)
#endif

#if _CLG_D3_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeD3, deviceDN<3>, 2, EFT_GaugeD3)
#endif
#if _CLG_D4_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeD4, deviceDN<4>, 2, EFT_GaugeD4)
#endif
#if _CLG_D8_GAUGE
__DEFINE_GAUGE_LINK(CFieldGaugeD8, deviceDN<8>, 2, EFT_GaugeD8)
#endif

template<typename deviceGauge, INT matrixN>
inline void CFieldGaugeLink<deviceGauge, matrixN>::ApplyPhaseC(const CField* other)
{
    const CFieldGaugeU1* phasec = dynamic_cast<const CFieldGaugeU1*>(other);
    if (NULL != phasec)
    {
        CCommonKernelField<deviceGauge>::ApplyPhaseC(m_pDeviceData, m_uiLinkeCount, phasec->m_pDeviceData);
        NotifyWritten();
        return;
    }

    appCrucial(_T("ApplyPhase with phase type not supported!\n"));
}

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_LINK_H_

//=============================================================================
// END OF FILE
//=============================================================================