//=============================================================================
// FILENAME : TestImproveHaloCommon.h
//
// DESCRIPTION:
// Shared deterministic-halo helpers for the Improve-1 MG tests (extracted
// from TestImprove1.cpp in I6 so TestImprove6.cpp can run the same
// global-site-pattern fill / halo-slot verification against owner buffers
// registered outside the field classes (clover PhidPhi, staple-cache sets).
// All functions are static inline: every translation unit keeps its own
// copy and unused helpers stay silent.
//
// REVISION:
//  [08/09/2026 Improve-1 I6 nbale]
//=============================================================================

#pragma once

#ifndef _TESTIMPROVEHALOCOMMON_H_
#define _TESTIMPROVEHALOCOMMON_H_

#pragma region deterministic halo helpers

/** Layout context for the halo slot decode / deterministic fill. */
struct SImprove1HaloCtx
{
    UINT uiLocalL[4];
    UINT uiGrid[4];
    UINT uiOffset[4];
    UINT uiGlobL[4];
    UINT uiVolume;
    UINT uiHaloSites;
    UINT uiHaloWidth;
};

/** Distinguished halo-tail value: never equal to a legal pattern value. */
static inline Real Improve1Sentinel() { return F(-777.0); }

/**
 * Build the layout context. Returns FALSE when nothing is split (halo count 0),
 * in which case every test here trivially passes (single-GPU build or -n 1).
 */
static inline UBOOL Improve1BuildCtx(SImprove1HaloCtx& ctx)
{
    const CIndexData* pIdx = appGetLattice()->m_pIndexCache;
    const CLGComm* pComm = appGetComm();
    if (NULL == pIdx || NULL == pComm)
    {
        return FALSE;
    }
    ctx.uiLocalL[0] = static_cast<UINT>(_HC_Lx);
    ctx.uiLocalL[1] = static_cast<UINT>(_HC_Ly);
    ctx.uiLocalL[2] = static_cast<UINT>(_HC_Lz);
    ctx.uiLocalL[3] = static_cast<UINT>(_HC_Lt);
    const UINT* pGrid = pComm->GpuGrid();
    const UINT* pOffset = pComm->GlobalOffset();
    const UINT* pGlobL = pComm->GlobalLattice();
    for (UINT i = 0; i < 4; ++i)
    {
        ctx.uiGrid[i] = pGrid[i];
        ctx.uiOffset[i] = pOffset[i];
        ctx.uiGlobL[i] = pGlobL[i];
    }
    ctx.uiVolume = _HC_Volume;
    ctx.uiHaloSites = pIdx->m_uiHaloSiteCount;
    ctx.uiHaloWidth = static_cast<UINT>(_HC_HaloWidth);
    return ctx.uiHaloSites > 0;
}

/**
 * Invert the CLGHaloLayout.h slot numbering (host mirror of the device
 * redirect, shared with TestMGHaloSelfExchange): slot -> the list of
 * (dir, side, layer) crossings plus the transverse LOCAL coordinates.
 * Fills tdir/tside/tlayer[0..count) and uiCoord[4] (transverse axes only;
 * crossing axes are left 0). Returns the codimension (1..4), or 0 when the
 * slot is out of layout. I8: codimension-4 slots follow the same rule.
 */
static inline UINT Improve1SlotDecode(UINT slot, const SImprove1HaloCtx& ctx,
    UINT tdir[4], UINT tside[4], UINT tlayer[4], UINT uiCoord[4])
{
    UINT uiScan = 0;
    uiCoord[0] = 0; uiCoord[1] = 0; uiCoord[2] = 0; uiCoord[3] = 0;
    UINT uiCount = 0;
    UBOOL bFound = FALSE;

    //FACE blocks: dir-major, side, layer, faceIdx.
    for (UINT d = 0; d < 4 && !bFound; ++d)
    {
        const UINT uiBlock = _haloFaceBlockSites(ctx.uiLocalL, ctx.uiGrid, ctx.uiVolume, ctx.uiHaloWidth, d);
        for (UINT s = 0; s < 2 && !bFound; ++s)
        {
            if (0 == uiBlock)
            {
                continue;
            }
            if (slot < uiScan + uiBlock)
            {
                const UINT uiInBlock = slot - uiScan;
                const UINT uiFaceVol = uiBlock / ctx.uiHaloWidth;
                const UINT uiLayer = uiInBlock / uiFaceVol;
                UINT uiFaceIdx = uiInBlock % uiFaceVol;
                for (INT a = 3; a >= 0; --a)
                {
                    if (static_cast<UINT>(a) == d)
                    {
                        continue;
                    }
                    uiCoord[a] = uiFaceIdx % ctx.uiLocalL[a];
                    uiFaceIdx /= ctx.uiLocalL[a];
                }
                tdir[0] = d; tside[0] = s; tlayer[0] = uiLayer;
                uiCount = 1;
                bFound = TRUE;
            }
            uiScan += uiBlock;
        }
    }

    //EDGE blocks: dir-pair-major (d1<d2), side-pair 2*s1+s2, layer1-major.
    for (UINT d1 = 0; d1 < 4 && !bFound; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4 && !bFound; ++d2)
        {
            if (ctx.uiGrid[d1] <= 1 || ctx.uiGrid[d2] <= 1)
            {
                continue;
            }
            const UINT uiEdgeVol = ctx.uiVolume / (ctx.uiLocalL[d1] * ctx.uiLocalL[d2]);
            const UINT uiBlock = ctx.uiHaloWidth * ctx.uiHaloWidth * uiEdgeVol;
            for (UINT sp = 0; sp < 4 && !bFound; ++sp)
            {
                if (slot < uiScan + uiBlock)
                {
                    const UINT uiInBlock = slot - uiScan;
                    const UINT uiLayer1 = uiInBlock / (ctx.uiHaloWidth * uiEdgeVol);
                    const UINT uiRest = uiInBlock % (ctx.uiHaloWidth * uiEdgeVol);
                    const UINT uiLayer0 = uiRest / uiEdgeVol;
                    UINT uiEdgeIdx = uiRest % uiEdgeVol;
                    for (INT a = 3; a >= 0; --a)
                    {
                        if (static_cast<UINT>(a) == d1 || static_cast<UINT>(a) == d2)
                        {
                            continue;
                        }
                        uiCoord[a] = uiEdgeIdx % ctx.uiLocalL[a];
                        uiEdgeIdx /= ctx.uiLocalL[a];
                    }
                    tdir[0] = d1; tside[0] = sp / 2; tlayer[0] = uiLayer0;
                    tdir[1] = d2; tside[1] = sp % 2; tlayer[1] = uiLayer1;
                    uiCount = 2;
                    bFound = TRUE;
                }
                uiScan += uiBlock;
            }
        }
    }

    //CORNER blocks: dir-triple-major, side-triple 4*s1+2*s2+s3, layer2-major.
    for (UINT d1 = 0; d1 < 4 && !bFound; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4 && !bFound; ++d2)
        {
            for (UINT d3 = d2 + 1; d3 < 4 && !bFound; ++d3)
            {
                if (ctx.uiGrid[d1] <= 1 || ctx.uiGrid[d2] <= 1 || ctx.uiGrid[d3] <= 1)
                {
                    continue;
                }
                const UINT uiCornerVol = ctx.uiVolume / (ctx.uiLocalL[d1] * ctx.uiLocalL[d2] * ctx.uiLocalL[d3]);
                const UINT uiBlock = ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth * uiCornerVol;
                for (UINT st = 0; st < 8 && !bFound; ++st)
                {
                    if (slot < uiScan + uiBlock)
                    {
                        const UINT uiInBlock = slot - uiScan;
                        const UINT uiLayer2 = uiInBlock / (ctx.uiHaloWidth * ctx.uiHaloWidth * uiCornerVol);
                        const UINT uiRest1 = uiInBlock % (ctx.uiHaloWidth * ctx.uiHaloWidth * uiCornerVol);
                        const UINT uiLayer1 = uiRest1 / (ctx.uiHaloWidth * uiCornerVol);
                        const UINT uiRest0 = uiRest1 % (ctx.uiHaloWidth * uiCornerVol);
                        const UINT uiLayer0 = uiRest0 / uiCornerVol;
                        const UINT uiCornerIdx = uiRest0 % uiCornerVol;
                        for (UINT a = 0; a < 4; ++a)
                        {
                            if (a != d1 && a != d2 && a != d3)
                            {
                                uiCoord[a] = uiCornerIdx;
                                break;
                            }
                        }
                        tdir[0] = d1; tside[0] = st / 4; tlayer[0] = uiLayer0;
                        tdir[1] = d2; tside[1] = (st / 2) % 2; tlayer[1] = uiLayer1;
                        tdir[2] = d3; tside[2] = st % 2; tlayer[2] = uiLayer2;
                        uiCount = 3;
                        bFound = TRUE;
                    }
                    uiScan += uiBlock;
                }
            }
        }
    }

    //HYPER-CORNER blocks: all four dirs out, side-quad 8*s0+4*s1+2*s2+s3,
    //layer3-major then layer2, layer1, layer0; the hyperplane has no free axis.
    if (!bFound && ctx.uiGrid[0] > 1 && ctx.uiGrid[1] > 1 && ctx.uiGrid[2] > 1 && ctx.uiGrid[3] > 1)
    {
        const UINT uiBlock = ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth;
        for (UINT sq = 0; sq < 16 && !bFound; ++sq)
        {
            if (slot < uiScan + uiBlock)
            {
                const UINT uiInBlock = slot - uiScan;
                const UINT uiLayer3 = uiInBlock / (ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth);
                const UINT uiRest2 = uiInBlock % (ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth);
                const UINT uiLayer2 = uiRest2 / (ctx.uiHaloWidth * ctx.uiHaloWidth);
                const UINT uiRest1 = uiRest2 % (ctx.uiHaloWidth * ctx.uiHaloWidth);
                const UINT uiLayer1 = uiRest1 / ctx.uiHaloWidth;
                const UINT uiLayer0 = uiRest1 % ctx.uiHaloWidth;
                tdir[0] = 0; tside[0] = sq / 8; tlayer[0] = uiLayer0;
                tdir[1] = 1; tside[1] = (sq / 4) % 2; tlayer[1] = uiLayer1;
                tdir[2] = 2; tside[2] = (sq / 2) % 2; tlayer[2] = uiLayer2;
                tdir[3] = 3; tside[3] = sq % 2; tlayer[3] = uiLayer3;
                uiCount = 4;
                bFound = TRUE;
            }
            uiScan += uiBlock;
        }
    }

    return uiCount;
}

/** Slot -> the GLOBAL 4-coordinate it represents (each axis in [0, GlobL)). */
static inline UBOOL Improve1SlotGlobalCoord(UINT slot, const SImprove1HaloCtx& ctx, INT gCoord[4])
{
    UINT tdir[4] = { 0, 0, 0, 0 };
    UINT tside[4] = { 0, 0, 0, 0 };
    UINT tlayer[4] = { 0, 0, 0, 0 };
    UINT uiCoord[4] = { 0, 0, 0, 0 };
    const UINT uiCount = Improve1SlotDecode(slot, ctx, tdir, tside, tlayer, uiCoord);
    if (0 == uiCount)
    {
        return FALSE;
    }

    for (UINT d = 0; d < 4; ++d)
    {
        gCoord[d] = static_cast<INT>(ctx.uiOffset[d] + uiCoord[d]);
    }
    for (UINT i = 0; i < uiCount; ++i)
    {
        const UINT d = tdir[i];
        if (0 == tside[i])
        {
            gCoord[d] = static_cast<INT>((ctx.uiOffset[d] + ctx.uiGlobL[d] - 1 - tlayer[i]) % ctx.uiGlobL[d]);
        }
        else
        {
            gCoord[d] = static_cast<INT>((ctx.uiOffset[d] + ctx.uiLocalL[d] + tlayer[i]) % ctx.uiGlobL[d]);
        }
    }
    return TRUE;
}

/**
 * Slot -> global linear site index (x slowest). Returns FALSE when out of
 * layout.
 */
static inline UBOOL Improve1SlotGlobalSite(UINT slot, const SImprove1HaloCtx& ctx, UINT& uiGlobalSiteOut)
{
    INT gCoord[4] = { 0, 0, 0, 0 };
    if (!Improve1SlotGlobalCoord(slot, ctx, gCoord))
    {
        return FALSE;
    }
    uiGlobalSiteOut = ((static_cast<UINT>(gCoord[0]) * ctx.uiGlobL[1] + static_cast<UINT>(gCoord[1]))
        * ctx.uiGlobL[2] + static_cast<UINT>(gCoord[2]))
        * ctx.uiGlobL[3] + static_cast<UINT>(gCoord[3]);
    return TRUE;
}

/**
 * Slot -> the LOCAL interior site the gather table must pack into it
 * (CIndexSquare::_kernalBakeHaloGatherSite): a crossing (dir, side, layer)
 * reads this rank's NEAR-boundary plane (side 0: local coord = layer; side 1:
 * local coord = L-1-layer); transverse axes keep their in-range coordinate.
 */
static inline UBOOL Improve1SlotExpectedSource(UINT slot, const SImprove1HaloCtx& ctx, UINT& uiSourceOut)
{
    UINT tdir[4] = { 0, 0, 0, 0 };
    UINT tside[4] = { 0, 0, 0, 0 };
    UINT tlayer[4] = { 0, 0, 0, 0 };
    UINT uiCoord[4] = { 0, 0, 0, 0 };
    const UINT uiCount = Improve1SlotDecode(slot, ctx, tdir, tside, tlayer, uiCoord);
    if (0 == uiCount)
    {
        return FALSE;
    }

    UINT uiSrc[4] = { uiCoord[0], uiCoord[1], uiCoord[2], uiCoord[3] };
    for (UINT i = 0; i < uiCount; ++i)
    {
        const UINT d = tdir[i];
        uiSrc[d] = (0 == tside[i]) ? tlayer[i] : (ctx.uiLocalL[d] - 1 - tlayer[i]);
    }
    uiSourceOut = ((uiSrc[0] * ctx.uiLocalL[1] + uiSrc[1]) * ctx.uiLocalL[2] + uiSrc[2])
        * ctx.uiLocalL[3] + uiSrc[3];
    return TRUE;
}

/**
 * Fill the local interior of pData with the deterministic global-site pattern
 * value(g) = (Real)g + fBias, and stamp every halo slot with the sentinel.
 * One host upload covers local + halo so the pre-refill halo state is
 * deterministic (no reliance on cudaMalloc scraps).
 */
static inline void Improve1FillLocalDeterministic(void* pData, UINT uiBytesPerSite, UINT uiElemPerSite,
    const SImprove1HaloCtx& ctx, Real fBias)
{
    const UINT uiGlobalVolume = ctx.uiGlobL[0] * ctx.uiGlobL[1] * ctx.uiGlobL[2] * ctx.uiGlobL[3];
    const CLGComm* pComm = appGetComm();

    Real* pGlobal = NULL;
    if (pComm->IsRoot())
    {
        pGlobal = (Real*)malloc(static_cast<size_t>(uiGlobalVolume) * uiBytesPerSite);
        for (UINT g = 0; g < uiGlobalVolume; ++g)
        {
            Real* pSite = pGlobal + static_cast<size_t>(g) * uiElemPerSite;
            for (UINT e = 0; e < uiElemPerSite; ++e)
            {
                pSite[e] = static_cast<Real>(g) + fBias;
            }
        }
    }
    Real* pHost = (Real*)malloc(static_cast<size_t>(ctx.uiVolume + ctx.uiHaloSites) * uiBytesPerSite);
    pComm->ScatterFieldFromRoot((const BYTE*)pGlobal, uiBytesPerSite, (BYTE*)pHost);
    if (NULL != pGlobal)
    {
        free(pGlobal);
    }

    Real* pHalo = pHost + static_cast<size_t>(ctx.uiVolume) * uiElemPerSite;
    for (UINT i = 0; i < ctx.uiHaloSites * uiElemPerSite; ++i)
    {
        pHalo[i] = Improve1Sentinel();
    }
    appSimpleCopyHD(pData, pHost, static_cast<size_t>(ctx.uiVolume + ctx.uiHaloSites) * uiBytesPerSite);
    free(pHost);
}

/** Fill the local interior with a single constant (halo left as-is). */
static inline void Improve1FillLocalConstant(void* pData, UINT uiBytesPerSite, UINT uiElemPerSite,
    const SImprove1HaloCtx& ctx, Real fValue)
{
    Real* pHost = (Real*)malloc(static_cast<size_t>(ctx.uiVolume) * uiBytesPerSite);
    for (UINT i = 0; i < ctx.uiVolume * uiElemPerSite; ++i)
    {
        pHost[i] = fValue;
    }
    appSimpleCopyHD(pData, pHost, static_cast<size_t>(ctx.uiVolume) * uiBytesPerSite);
    free(pHost);
}

/**
 * Verify every halo slot of pData holds value(globalSite) + fBias in ALL
 * elements. Returns the error count.
 */
static inline UINT Improve1CheckHalo(const void* pData, UINT uiBytesPerSite, UINT uiElemPerSite,
    const SImprove1HaloCtx& ctx, Real fBias, const TCHAR* sWhat)
{
    UINT uiError = 0;
    Real* pHost = (Real*)malloc(static_cast<size_t>(ctx.uiVolume + ctx.uiHaloSites) * uiBytesPerSite);
    appSimpleCopyDH(pHost, pData, static_cast<size_t>(ctx.uiVolume + ctx.uiHaloSites) * uiBytesPerSite);

    for (UINT slot = 0; slot < ctx.uiHaloSites; ++slot)
    {
        UINT uiGExp = 0;
        if (!Improve1SlotGlobalSite(slot, ctx, uiGExp))
        {
            ++uiError;
            LastProbem(_T("Halo slot out of layout"));
            break;
        }
        const Real* pSlot = pHost + (static_cast<size_t>(ctx.uiVolume) + slot) * uiElemPerSite;
        for (UINT e = 0; e < uiElemPerSite; ++e)
        {
            if (appAbs(static_cast<DOUBLE>(pSlot[e]) - static_cast<DOUBLE>(uiGExp) - static_cast<DOUBLE>(fBias)) > 0.5)
            {
                ++uiError;
                CCString sProblem;
                sProblem.Format(_T("%s slot %u elem %u: expect %.1f, got %.1f"),
                    sWhat, slot, e, static_cast<DOUBLE>(uiGExp) + static_cast<DOUBLE>(fBias),
                    static_cast<DOUBLE>(pSlot[e]));
                appGeneral(_T("%s\n"), sProblem.c_str());
                LastProbem(sProblem);
                if (uiError >= 10)
                {
                    appGeneral(_T("Too many errors, stopping verification.\n"));
                    free(pHost);
                    return uiError;
                }
            }
        }
    }
    free(pHost);
    return uiError;
}

/**
 * I3/I4 migration point: the refill trigger. Handle-based (multi-GPU-improve1.md
 * 3.6): refill through the field's own bound halo handle -- pointer, stride and
 * MPI tag id all come from the bound descriptor.
 */
static inline void Improve1TriggerRefill(CField* pField)
{
    appGetHaloManager()->RefillHalo(*pField->GetHaloBufferHandle());
}

#pragma endregion

#endif //#ifndef _TESTIMPROVEHALOCOMMON_H_
