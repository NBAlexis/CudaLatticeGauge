//=============================================================================
// FILENAME : CActionGaugePlaquetteCylinder.cu
//
// DESCRIPTION:
// This is the gauge action in cylindrical coordinates.
// The lattice directions 0,1,2,3 are interpreted as (r, phi, z, t).
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteCylinder.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteCylinder)

#pragma region kernel

//the weight of a plaquette plane, 1/r for planes containing phi, r for other planes
__device__ __inline__ Real _deviceCylinderW(Real fR, UBOOL bHasPhi)
{
    return bHasPhi ? (F(1.0) / fR) : fR;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_UseCloverCylinder(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE* __restrict__ fBetaOverN,
    DOUBLE* results,
#else
    const Real* __restrict__ fBetaOverN,
    Real* results,
#endif
    Real fRStart,
    Real fDeltaR
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

#if !_CLG_DOUBLEFLOAT
    DOUBLE fRes = 0.0;
#else
    Real fRes = F(0.0);
#endif
    //Multi-GPU: the radial coordinate r AND the per-bin beta index must be
    //GLOBAL (P4-1.5/1.4-R1); identity on single-GPU (offset 0). The device
    //beta array is built over the global radial range (see Initial), so
    //global indexing is always in range. sSite4 itself stays local for the
    //index-table lookups below.
    const INT iRg = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex)).x;
    const Real fR = fRStart + fDeltaR * static_cast<Real>(iRg);
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            const UBOOL bHasPhi = (1 == byDir1) || (1 == byDir2);
            const Real fK = static_cast<Real>(fBetaOverN[iRg]) * _deviceCylinderW(fR, bHasPhi);
            fRes += fK * (F(3.0) - F(0.25) * _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId));
        }
    }
    results[uiSiteIndex] = fRes;
}

//plaquette form energy, every plaquette is weighted by the explicit 4-corner coupling,
//this is different from the clover form at the Dirichlet boundary
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_Cylinder(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
#if !_CLG_DOUBLEFLOAT
    const DOUBLE* __restrict__ fBetaOverN,
    DOUBLE* results,
#else
    const Real* __restrict__ fBetaOverN,
    Real* results,
#endif
    Real fRStart,
    Real fDeltaR
)
{
    intokernalInt4;
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqCountAllSite = plaqCountPerSite * plaqLength;
#endif

#if !_CLG_DOUBLEFLOAT
    DOUBLE fRes = 0.0;
#else
    Real fRes = F(0.0);
#endif
    const UINT indexSkip = plaqCountPerSite * plaqLength * uiSiteIndex;

    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        //count Dirichlet edges, the same as the force kernel does
        BYTE diricCount = 0;
        const SIndex first = pCachedIndex[i * plaqLength + indexSkip];
        const BYTE mu = first.m_byDir;
        deviceSU3 toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.IsDirichlet())
        {
            ++diricCount;
        }
        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        BYTE nu = 0;
        for (BYTE j = 1; j < plaqLength; ++j)
        {
            const SIndex nextlink = pCachedIndex[i * plaqLength + j + indexSkip];
            if (1 == j)
            {
                nu = nextlink.m_byDir;
            }
            deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));
            if (nextlink.IsDirichlet())
            {
                ++diricCount;
            }
            if (nextlink.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        //If more than 3(including 3) of the edges are Dirichlet, the plaqutte dose NOT exist.
        if (diricCount >= plaqLength - 1)
        {
            continue;
        }

        const UBOOL bHasPhi = (1 == mu) || (1 == nu);
        Real fW;
        if (0 == mu || 0 == nu)
        {
            //plane contains r, corners are at r(x) and r(x+1), use the same corrections as the force kernel
            //Multi-GPU: corners and the boundary logic must use GLOBAL
            //coordinates/extent (P4-1.5); identity on single-GPU.
            const INT iR1g = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex)).x;
            INT r2 = iR1g + 1;
            INT r1 = iR1g;
            if (0 == r1)
            {
                if (diricCount > 0)
                {
                    r1 = r2;
                }
            }
            if (r2 >= static_cast<INT>(_DC_GlobalLx))
            {
                if (diricCount > 0)
                {
                    r2 = static_cast<INT>(_DC_GlobalLx) - 1;
                }
                else
                {
                    r2 = 0;
                }
            }
            fW = F(0.5) * (static_cast<Real>(fBetaOverN[r1]) * _deviceCylinderW(fRStart + fDeltaR * static_cast<Real>(r1), bHasPhi)
                + static_cast<Real>(fBetaOverN[r2]) * _deviceCylinderW(fRStart + fDeltaR * static_cast<Real>(r2), bHasPhi));
        }
        else
        {
            //all corners are at r(x)
            const INT iRg = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex)).x;
            fW = static_cast<Real>(fBetaOverN[iRg]) * _deviceCylinderW(fRStart + fDeltaR * static_cast<Real>(iRg), bHasPhi);
        }
        if (diricCount > 0)
        {
            fW = fW * F(0.5);
        }
        fRes += (F(3.0) - toAdd.ReTr()) * fW;
    }

    results[uiSiteIndex] = fRes;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndexCylinder(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceSU3* pStapleData, //can be NULL
    deviceSU3* pForceData,
    const DOUBLE* __restrict__ betaOverN,
    Real fRStart,
    Real fDeltaR
)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, static_cast<BYTE>(idir)))
        {
            continue;
        }

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCountPerLink; ++i)
        {
            //the r slices and whether the plane contains phi are per staple
            INT r1 = 0;
            INT r2 = 0;
            BYTE rNotFound = 1;
            UBOOL bHasPhi = (1 == idir);
            BYTE diricCount = 0;
            const SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAllLink];
            if (first.IsDirichlet())
            {
                ++diricCount;
            }
            deviceSU3 toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

            //Multi-GPU: the radial corner coordinate must be GLOBAL and the
            //link SIndex may be halo-redirected (P4-1.5); identity on
            //single-GPU. Improve-1 (3.7): 32-bit global coordinate.
            r1 = _deviceSIndexToGlobalInt4(first).x;
            if (0 == first.m_byDir)
            {
                r2 = r1 + 1;
                rNotFound = 0;
            }
            if (1 == first.m_byDir)
            {
                bHasPhi = TRUE;
            }

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAllLink];
                if (nextlink.IsDirichlet())
                {
                    ++diricCount;
                }
                deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }

                if (rNotFound && nextlink.m_byDir == 0)
                {
                    rNotFound = 0;
                    r1 = _deviceSIndexToGlobalInt4(nextlink).x;
                    r2 = r1 + 1;
                }
                if (1 == nextlink.m_byDir)
                {
                    bHasPhi = TRUE;
                }
            }

            if (diricCount < plaqLengthm1)
            {
                // If more than 3(including 3) of the edges are Dirichlet,
                // the plaqutte dose NOT exist.
                if (rNotFound)
                {
                    //rNotFound means this is a plaqutte not have r bounds, all corners are at r(n)
                    const Real fR = fRStart + fDeltaR * static_cast<Real>(r1);
                    toAdd.MulReal(static_cast<Real>(betaOverN[r1]) * _deviceCylinderW(fR, bHasPhi) * F(-0.5));
                }
                else
                {
                    if (0 == r1)
                    {
                        if (diricCount > 0)
                        {
                            r1 = r2;
                        }
                    }

                    //Multi-GPU: the boundary logic must use the GLOBAL extent
                    //(P4-1.5); identity on single-GPU (GlobalLx == Lx).
                    if (r2 >= static_cast<INT>(_DC_GlobalLx))
                    {
                        if (diricCount > 0)
                        {
                            r2 = static_cast<INT>(_DC_GlobalLx) - 1;
                        }
                        else
                        {
                            r2 = 0;
                        }
                    }
                    const Real fR1 = fRStart + fDeltaR * static_cast<Real>(r1);
                    const Real fR2 = fRStart + fDeltaR * static_cast<Real>(r2);
                    toAdd.MulReal(F(-0.25) * (static_cast<Real>(betaOverN[r1]) * _deviceCylinderW(fR1, bHasPhi)
                        + static_cast<Real>(betaOverN[r2]) * _deviceCylinderW(fR2, bHasPhi)));
                }
                if (diricCount > 0)
                {
                    toAdd.MulReal(F(0.5));
                }
                res.Add(toAdd);
            }
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, static_cast<BYTE>(idir)))
        {
            //force is additive
            pForceData[linkIndex].Add(res);
        }
    }
}

#pragma endregion

CActionGaugePlaquetteCylinder::CActionGaugePlaquetteCylinder()
    : CAction()
    , m_pDeviceBetaArray(NULL)
    , m_fRStart(F(1.0))
    , m_fREnd(F(1.0))
    , m_fDeltaR(F(1.0))
    , m_bCloverEnergy(FALSE)
    , m_uiPlaqutteCount(0)
{
}

DOUBLE CActionGaugePlaquetteCylinder::CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId]);

    preparethread;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_UseCloverCylinder, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        m_pDeviceBetaArray,
        _D_RealThreadBuffer,
        m_fRStart,
        m_fDeltaR);

    DOUBLE fRet = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    //The ThreadBufferSum reduces only this rank's local sub-lattice; sum
    //across ranks for the global action (P4-1.5). No-op on a lone rank.
    appGlobalSum(fRet);
    return fRet;
}

DOUBLE CActionGaugePlaquetteCylinder::CalculatePlaqutteEnergyUsePlaqutte(const CFieldGaugeSU3* pGauge) const
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId]);

    preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_Cylinder, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        m_pDeviceBetaArray,
        _D_RealThreadBuffer,
        m_fRStart,
        m_fDeltaR);
#else
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_Cylinder, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        m_pDeviceBetaArray,
        _D_RealThreadBuffer,
        m_fRStart,
        m_fDeltaR);
#endif

    DOUBLE fRet = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    appGlobalSum(fRet);
    return fRet;
}

void CActionGaugePlaquetteCylinder::CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const
{
    preparethread;

    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId]);
    deviceSU3* emptyStaple = NULL;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndexCylinder, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        emptyStaple,
        pForce->m_pDeviceData,
        m_pDeviceBetaArray,
        m_fRStart,
        m_fDeltaR);
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndexCylinder, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        emptyStaple,
        pForce->m_pDeviceData,
        m_pDeviceBetaArray,
        m_fRStart,
        m_fDeltaR);
#endif
}

void CActionGaugePlaquetteCylinder::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
        if (NULL == pGaugeSU3)
        {
            appCrucial(_T("CActionGaugePlaquetteCylinder must be CFieldGaugeSU3!\n"));
            return;
        }

        appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId]);
        if (m_bCloverEnergy)
        {
            m_fLastEnergy = CalculatePlaqutteEnergyUseClover(pGaugeSU3);
        }
        else
        {
            m_fLastEnergy = CalculatePlaqutteEnergyUsePlaqutte(pGaugeSU3);
        }
    }
}

void CActionGaugePlaquetteCylinder::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    param.FetchValueReal(_T("RStart"), m_fRStart);
    param.FetchValueReal(_T("DeltaR"), m_fDeltaR);
    if (!param.FetchValueReal(_T("REnd"), m_fREnd))
    {
        m_fREnd = m_fRStart + m_fDeltaR * static_cast<Real>(_HC_Lx);
    }
    if (m_fRStart < _CLG_FLT_EPSILON)
    {
        appCrucial(_T("CActionGaugePlaquetteCylinder: r = 0 is a singularity, RStart (%f) must be positive!\n"), m_fRStart);
    }

    INT iUsing4Plaq = 0;
    if (param.FetchValueINT(_T("CloverEnergy"), iUsing4Plaq))
    {
        if (1 == iUsing4Plaq)
        {
            m_bCloverEnergy = TRUE;
        }
    }

    param.FetchValueArrayDOUBLE(_T("Beta"), m_fBetaArray);
#if _CLG_MULTI_GPU
    //Multi-GPU: kernels index the beta array with the GLOBAL radial coordinate,
    //so every rank builds the array over the GLOBAL radial range (the yaml
    //Beta list is per global bin). Single-GPU: GlobalLattice == local lattice,
    //identical layout (P4-1.5).
    UINT uiGlobalLx = static_cast<UINT>(_HC_Lx);
    if (NULL != appGetComm())
    {
        uiGlobalLx = appGetComm()->GlobalLattice()[0];
    }
    for (UINT i = 0; i < uiGlobalLx; ++i)
#else
    for (INT i = 0; i < _HC_Lxi; ++i)
#endif
    {
        if (static_cast<INT>(i) < m_fBetaArray.Num())
        {
            m_fBetaArray[i] = m_fBetaArray[i] / static_cast<DOUBLE>(GetDefaultMatrixN());
        }
        else
        {
            if (0 == i)
            {
                m_fBetaArray.AddItem(5.0 / static_cast<DOUBLE>(GetDefaultMatrixN()));
            }
            else
            {
                m_fBetaArray.AddItem(m_fBetaArray[m_fBetaArray.Num() - 1]);
            }
        }
    }

    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

#if _CLG_MULTI_GPU
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceBetaArray, sizeof(DOUBLE) * uiGlobalLx));
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(DOUBLE) * uiGlobalLx, cudaMemcpyHostToDevice));
#else
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceBetaArray, sizeof(DOUBLE) * _HC_Lx));
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(DOUBLE) * _HC_Lx, cudaMemcpyHostToDevice));
#endif
}

UBOOL CActionGaugePlaquetteCylinder::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteCylinder must be CFieldGaugeSU3!\n"));
        return FALSE;
    }

    CalculateForceAndStaple(pGaugeSU3, pForceSU3);
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquetteCylinder::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteCylinder must be CFieldGaugeSU3!\n"));
        return m_fNewEnergy;
    }
    m_fNewEnergy = m_bCloverEnergy ? CalculatePlaqutteEnergyUseClover(pGaugeSU3) : CalculatePlaqutteEnergyUsePlaqutte(pGaugeSU3);
    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteCylinder::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("RStart : ") + appToString(m_fRStart) + _T("\n");
    sRet = sRet + tab + _T("REnd : ") + appToString(m_fREnd) + _T("\n");
    sRet = sRet + tab + _T("DeltaR : ") + appToString(m_fDeltaR) + _T("\n");
    sRet = sRet + tab + _T("Clover : ") + appToString(m_bCloverEnergy) + _T("\n");
    sRet = sRet + tab + _T("Beta : ");
    for (INT i = 0; i < _HC_Lxi; ++i)
    {
        sRet = sRet + appToString(m_fBetaArray[i]) + _T(", ");
    }
    sRet = sRet + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
