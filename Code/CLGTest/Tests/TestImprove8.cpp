//=============================================================================
// FILENAME : TestImprove8.cpp
//
// DESCRIPTION:
// Improve-1 (multi-GPU-improve1.md) I8 verification tests: the baked halo
// tables against analytic references.
//
// - TestMGHaloBakeTables: every halo slot's gather-source SIndex and baked
//   32-bit global coordinate match the analytic mirror (I8 gate 1), and every
//   local table entry is this rank's offset + local coordinate.
// - TestMGHaloRadiusReject: scanning the gauge field's big-cell SIndex table,
//   a split-direction crossing beyond HaloWidth is the invalid sentinel
//   (SIndex::_kInvalidSiteIndex), never a silent local wrap (I8 gate 2). An
//   in-width crossing lands in [Volume, Volume+haloSites) and its slot maps
//   back to the same global site the mirror predicts.
// - TestMGHaloHyperCorner: codimension-4 slots exist exactly when all four
//   directions are split, and their content verifies after a refill
//   (I8 gate 4; run with `runtest.sh TestMGHaloHyperCorner 1` -- metadata grid
//   [2,2,2,2] = 16 oversubscribed ranks, HaloWidth 1 in the yaml block).
// - TestMGGlobalCoordinate32: the whole global-coordinate table on a lattice
//   whose global T extent exceeds 127 (I8 gate 3); the extent is asserted so
//   the test can never silently become vacuous.
//
// The first three run harmlessly single-process / single-GPU (halo count 0
// parts trivially pass; the local-table and big-cell scans still run as
// identity checks). TestMGGlobalCoordinate32 is different: its global T = 160
// only builds when the grid splits T (baked coordinates are 8-bit, local
// extent must stay <= 119), so the yaml sets RequireSplit : 1 and the driver
// skips it on an unsplit grid (single-GPU suite or -n 1) instead of hanging
// the bake kernel.
//
// REVISION:
//  [08/09/2026 Improve-1 I8 nbale]
//=============================================================================

#include "CLGTest.h"

#include "TestImproveHaloCommon.h"

#pragma region shared table verification

/** Fetch this rank's grid/offset/global lengths; all-1/zeros/local on no comm. */
static void Improve8FetchTopology(UINT uiGrid[4], UINT uiOffset[4], UINT uiGlobL[4], UINT uiLocalL[4])
{
    uiLocalL[0] = static_cast<UINT>(_HC_Lx);
    uiLocalL[1] = static_cast<UINT>(_HC_Ly);
    uiLocalL[2] = static_cast<UINT>(_HC_Lz);
    uiLocalL[3] = static_cast<UINT>(_HC_Lt);
    const CLGComm* pComm = appGetComm();
    for (UINT d = 0; d < 4; ++d)
    {
        uiGrid[d] = 1;
        uiOffset[d] = 0;
        uiGlobL[d] = uiLocalL[d];
    }
    if (NULL != pComm)
    {
        const UINT* pGrid = pComm->GpuGrid();
        const UINT* pOffset = pComm->GlobalOffset();
        const UINT* pGlobL = pComm->GlobalLattice();
        for (UINT d = 0; d < 4; ++d)
        {
            uiGrid[d] = pGrid[d];
            uiOffset[d] = pOffset[d];
            uiGlobL[d] = pGlobL[d];
        }
    }
}

/**
 * Verify CIndexData::m_pGlobalCoordinateTable: local entries must be
 * offset + local coordinate; halo entries must equal the analytic mirror.
 * When bIncludeGatherSource is TRUE, also verify every
 * CIndexData::m_pHaloGatherIndex entry against Improve1SlotExpectedSource.
 */
static UINT Improve8VerifyCoordinateTable(UBOOL bIncludeGatherSource, const TCHAR* sWhat)
{
    UINT uiError = 0;
    const CIndexData* pIdx = appGetLattice()->m_pIndexCache;
    if (NULL == pIdx || NULL == pIdx->m_pGlobalCoordinateTable)
    {
        LastProbem(_T("Global coordinate table not baked"));
        return 1;
    }

    UINT uiGrid[4];
    UINT uiOffset[4];
    UINT uiGlobL[4];
    UINT uiLocalL[4];
    Improve8FetchTopology(uiGrid, uiOffset, uiGlobL, uiLocalL);
    const UINT uiVolume = _HC_Volume;
    const UINT uiHalo = pIdx->m_uiHaloSiteCount;

    SInt4* pTable = (SInt4*)malloc(sizeof(SInt4) * (uiVolume + uiHalo));
    appSimpleCopyDH(pTable, pIdx->m_pGlobalCoordinateTable, sizeof(SInt4) * (uiVolume + uiHalo));

    //1. Local entries: offset + local coordinate.
    UINT uiReported = 0;
    for (UINT x = 0; x < uiLocalL[0]; ++x)
    {
        for (UINT y = 0; y < uiLocalL[1]; ++y)
        {
            for (UINT z = 0; z < uiLocalL[2]; ++z)
            {
                for (UINT w = 0; w < uiLocalL[3]; ++w)
                {
                    const UINT uiIdx = ((x * uiLocalL[1] + y) * uiLocalL[2] + z) * uiLocalL[3] + w;
                    const SInt4 sExpected(
                        static_cast<INT>(uiOffset[0] + x), static_cast<INT>(uiOffset[1] + y),
                        static_cast<INT>(uiOffset[2] + z), static_cast<INT>(uiOffset[3] + w));
                    if (!(pTable[uiIdx] == sExpected))
                    {
                        ++uiError;
                        if (uiReported < 10)
                        {
                            ++uiReported;
                            appGeneral(_T("%s local (%u,%u,%u,%u): expect [%d,%d,%d,%d], got [%d,%d,%d,%d]\n"),
                                sWhat, x, y, z, w,
                                sExpected.x, sExpected.y, sExpected.z, sExpected.w,
                                pTable[uiIdx].x, pTable[uiIdx].y, pTable[uiIdx].z, pTable[uiIdx].w);
                        }
                    }
                }
            }
        }
    }

    //2. Halo slot entries (+ gather source when requested).
    SImprove1HaloCtx ctx;
    if (0 == uiError && uiHalo > 0 && Improve1BuildCtx(ctx))
    {
        SIndex* pGather = NULL;
        if (bIncludeGatherSource)
        {
            if (NULL == pIdx->m_pHaloGatherIndex)
            {
                LastProbem(_T("Halo gather index table not baked"));
                free(pTable);
                return 1;
            }
            pGather = (SIndex*)malloc(sizeof(SIndex) * uiHalo);
            appSimpleCopyDH(pGather, pIdx->m_pHaloGatherIndex, sizeof(SIndex) * uiHalo);
        }

        for (UINT slot = 0; slot < uiHalo; ++slot)
        {
            INT gExpected[4] = { 0, 0, 0, 0 };
            if (!Improve1SlotGlobalCoord(slot, ctx, gExpected))
            {
                ++uiError;
                LastProbem(_T("Halo slot out of layout"));
                break;
            }
            const SInt4 sGot = pTable[uiVolume + slot];
            if (!(sGot == SInt4(gExpected[0], gExpected[1], gExpected[2], gExpected[3])))
            {
                ++uiError;
                if (uiReported < 10)
                {
                    ++uiReported;
                    appGeneral(_T("%s slot %u coord: expect [%d,%d,%d,%d], got [%d,%d,%d,%d]\n"),
                        sWhat, slot, gExpected[0], gExpected[1], gExpected[2], gExpected[3],
                        sGot.x, sGot.y, sGot.z, sGot.w);
                }
            }

            if (bIncludeGatherSource)
            {
                UINT uiExpectedSource = 0;
                if (!Improve1SlotExpectedSource(slot, ctx, uiExpectedSource))
                {
                    ++uiError;
                    LastProbem(_T("Halo slot out of layout (source)"));
                    break;
                }
                if (pGather[slot].m_uiSiteIndex != uiExpectedSource)
                {
                    ++uiError;
                    if (uiReported < 10)
                    {
                        ++uiReported;
                        appGeneral(_T("%s slot %u source: expect local %u, got %u\n"),
                            sWhat, slot, uiExpectedSource, pGather[slot].m_uiSiteIndex);
                    }
                }
            }
        }
        if (NULL != pGather)
        {
            free(pGather);
        }
    }
    else if (0 == uiError && 0 == uiHalo)
    {
        appGeneral(_T("Halo count is 0 (nothing split), slot part trivially passes.\n"));
    }

    if (uiError > uiReported)
    {
        LastProbem(_T("Coordinate table verification failed"));
    }
    free(pTable);
    return uiError;
}

#pragma endregion

/**
 * I8 gate 1: all halo slots' source SIndex and target global coordinate
 * against the analytic reference, plus the local part of the table.
 */
UINT TestMGHaloBakeTables(CParameters&)
{
    const UINT uiError = Improve8VerifyCoordinateTable(TRUE, _T("BakeTables"));
    if (0 == uiError)
    {
        appGeneral(_T("Halo bake tables verified against the analytic reference.\n"));
    }
    return uiError;
}
___REGIST_TEST(TestMGHaloBakeTables, MG, TestMGHaloBakeTables, HaloBakeTables, _TEST_MULTIGPU);

/**
 * I8 gate 2: an over-wide split-direction crossing must be the invalid
 * sentinel in the baked big-cell SIndex table, never a silent local wrap.
 * Scans the whole table of the registered gauge field (Torus bake):
 *   - fully in-range cell        -> valid local site index (== linear index);
 *   - crossing within HaloWidth  -> valid halo slot, and the slot maps back
 *                                   to the same global site as the mirror;
 *   - crossing beyond HaloWidth  -> SIndex::_kInvalidSiteIndex.
 * With nothing split every cell resolves locally (classic wrap), so the scan
 * is an identity check that must find zero invalid entries.
 */
UINT TestMGHaloRadiusReject(CParameters&)
{
    UINT uiError = 0;
    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]);
    const CIndexData* pIdx = appGetLattice()->m_pIndexCache;
    if (NULL == pGauge || NULL == pIdx)
    {
        LastProbem(_T("Need a gauge field and an index cache"));
        return 1;
    }
    const SIndex* pDeviceTable = pIdx->m_pIndexPositionToSIndex[pGauge->m_byFieldId];
    if (NULL == pDeviceTable)
    {
        LastProbem(_T("Big-cell SIndex table not baked for the gauge field"));
        return 1;
    }

    UINT uiGrid[4];
    UINT uiOffset[4];
    UINT uiGlobL[4];
    UINT uiLocalL[4];
    Improve8FetchTopology(uiGrid, uiOffset, uiGlobL, uiLocalL);
    const UINT uiVolume = _HC_Volume;
    const UINT uiHalo = pIdx->m_uiHaloSiteCount;
    const UINT uiHaloWidth = static_cast<UINT>(_HC_HaloWidth);

    const INT iEdge = CIndexData::kCacheIndexEdge;
    const UINT uiBigL[4] = {
        uiLocalL[0] + 2 * static_cast<UINT>(iEdge), uiLocalL[1] + 2 * static_cast<UINT>(iEdge),
        uiLocalL[2] + 2 * static_cast<UINT>(iEdge), uiLocalL[3] + 2 * static_cast<UINT>(iEdge) };
    const UINT uiBigCount = uiBigL[0] * uiBigL[1] * uiBigL[2] * uiBigL[3];

    SIndex* pTable = (SIndex*)malloc(sizeof(SIndex) * uiBigCount);
    appSimpleCopyDH(pTable, pDeviceTable, sizeof(SIndex) * uiBigCount);

    SImprove1HaloCtx ctx;
    const UBOOL bSplit = Improve1BuildCtx(ctx);

    UINT uiInvalidCount = 0;
    UINT uiHaloRefCount = 0;
    UINT uiReported = 0;
    for (UINT bx = 0; bx < uiBigL[0]; ++bx)
    {
        for (UINT by = 0; by < uiBigL[1]; ++by)
        {
            for (UINT bz = 0; bz < uiBigL[2]; ++bz)
            {
                for (UINT bw = 0; bw < uiBigL[3]; ++bw)
                {
                    const UINT uiB[4] = { bx, by, bz, bw };
                    INT iRaw[4];
                    UINT uiMaxSplitOut = 0;
                    UBOOL bInRange = TRUE;
                    for (UINT d = 0; d < 4; ++d)
                    {
                        iRaw[d] = static_cast<INT>(uiB[d]) - iEdge;
                        UINT uiOut = 0;
                        if (iRaw[d] < 0)
                        {
                            uiOut = static_cast<UINT>(-iRaw[d]);
                            bInRange = FALSE;
                        }
                        else if (iRaw[d] >= static_cast<INT>(uiLocalL[d]))
                        {
                            uiOut = static_cast<UINT>(iRaw[d] - static_cast<INT>(uiLocalL[d]) + 1);
                            bInRange = FALSE;
                        }
                        if (uiGrid[d] > 1 && uiOut > uiMaxSplitOut)
                        {
                            uiMaxSplitOut = uiOut;
                        }
                    }
                    const UINT uiBigIdx = ((bx * uiBigL[1] + by) * uiBigL[2] + bz) * uiBigL[3] + bw;
                    const SIndex& sEntry = pTable[uiBigIdx];
                    const UBOOL bInvalid = (SIndex::_kInvalidSiteIndex == sEntry.m_uiSiteIndex);

                    if (bInRange)
                    {
                        const UINT uiExpect = ((static_cast<UINT>(iRaw[0]) * uiLocalL[1] + static_cast<UINT>(iRaw[1]))
                            * uiLocalL[2] + static_cast<UINT>(iRaw[2])) * uiLocalL[3] + static_cast<UINT>(iRaw[3]);
                        if (bInvalid || sEntry.m_uiSiteIndex != uiExpect)
                        {
                            ++uiError;
                            if (uiReported < 10)
                            {
                                ++uiReported;
                                appGeneral(_T("RadiusReject in-range cell (%d,%d,%d,%d): expect local %u, got %u\n"),
                                    iRaw[0], iRaw[1], iRaw[2], iRaw[3], uiExpect, sEntry.m_uiSiteIndex);
                            }
                        }
                        continue;
                    }

                    if (uiMaxSplitOut > uiHaloWidth)
                    {
                        //Over-wide split crossing: must be the invalid sentinel.
                        ++uiInvalidCount;
                        if (!bInvalid)
                        {
                            ++uiError;
                            if (uiReported < 10)
                            {
                                ++uiReported;
                                appGeneral(_T("RadiusReject over-wide cell (%d,%d,%d,%d): expect INVALID, got %u (silent wrap!)\n"),
                                    iRaw[0], iRaw[1], iRaw[2], iRaw[3], sEntry.m_uiSiteIndex);
                            }
                        }
                        continue;
                    }

                    //In-width crossing (or non-split wrap): never invalid.
                    if (bInvalid)
                    {
                        ++uiError;
                        if (uiReported < 10)
                        {
                            ++uiReported;
                            appGeneral(_T("RadiusReject cell (%d,%d,%d,%d) unexpectedly INVALID\n"),
                                iRaw[0], iRaw[1], iRaw[2], iRaw[3]);
                        }
                        continue;
                    }

                    if (0 == uiMaxSplitOut)
                    {
                        //Out only in non-split directions: classic local wrap.
                        UINT uiWrapped[4];
                        for (UINT d = 0; d < 4; ++d)
                        {
                            INT iW = iRaw[d] % static_cast<INT>(uiLocalL[d]);
                            if (iW < 0)
                            {
                                iW += static_cast<INT>(uiLocalL[d]);
                            }
                            uiWrapped[d] = static_cast<UINT>(iW);
                        }
                        const UINT uiExpect = ((uiWrapped[0] * uiLocalL[1] + uiWrapped[1])
                            * uiLocalL[2] + uiWrapped[2]) * uiLocalL[3] + uiWrapped[3];
                        if (sEntry.m_uiSiteIndex != uiExpect)
                        {
                            ++uiError;
                            if (uiReported < 10)
                            {
                                ++uiReported;
                                appGeneral(_T("RadiusReject non-split wrap cell (%d,%d,%d,%d): expect %u, got %u\n"),
                                    iRaw[0], iRaw[1], iRaw[2], iRaw[3], uiExpect, sEntry.m_uiSiteIndex);
                            }
                        }
                        continue;
                    }

                    //Split crossing within width: halo slot; the slot maps back
                    //to the same global site as ((offset + raw) mod GlobalL).
                    if (!bSplit || sEntry.m_uiSiteIndex < uiVolume || sEntry.m_uiSiteIndex >= uiVolume + uiHalo)
                    {
                        ++uiError;
                        if (uiReported < 10)
                        {
                            ++uiReported;
                            appGeneral(_T("RadiusReject cell (%d,%d,%d,%d): expect halo slot in [%u,%u), got %u\n"),
                                iRaw[0], iRaw[1], iRaw[2], iRaw[3], uiVolume, uiVolume + uiHalo, sEntry.m_uiSiteIndex);
                        }
                        continue;
                    }
                    ++uiHaloRefCount;
                    const UINT uiSlot = sEntry.m_uiSiteIndex - uiVolume;
                    UINT uiMirrorSite = 0;
                    if (!Improve1SlotGlobalSite(uiSlot, ctx, uiMirrorSite))
                    {
                        ++uiError;
                        LastProbem(_T("RadiusReject: halo slot out of layout"));
                        continue;
                    }
                    UINT uiGCoord[4];
                    for (UINT d = 0; d < 4; ++d)
                    {
                        INT iG = (static_cast<INT>(uiOffset[d]) + iRaw[d]) % static_cast<INT>(uiGlobL[d]);
                        if (iG < 0)
                        {
                            iG += static_cast<INT>(uiGlobL[d]);
                        }
                        uiGCoord[d] = static_cast<UINT>(iG);
                    }
                    const UINT uiExpectSite = ((uiGCoord[0] * uiGlobL[1] + uiGCoord[1])
                        * uiGlobL[2] + uiGCoord[2]) * uiGlobL[3] + uiGCoord[3];
                    if (uiMirrorSite != uiExpectSite)
                    {
                        ++uiError;
                        if (uiReported < 10)
                        {
                            ++uiReported;
                            appGeneral(_T("RadiusReject cell (%d,%d,%d,%d): slot %u maps to global %u, expect %u\n"),
                                iRaw[0], iRaw[1], iRaw[2], iRaw[3], uiSlot, uiMirrorSite, uiExpectSite);
                        }
                    }
                }
            }
        }
    }

    appGeneral(_T("RadiusReject scan: %u cells, %u invalid (over-wide), %u halo references, %u errors.\n"),
        uiBigCount, uiInvalidCount, uiHaloRefCount, uiError);
    if (bSplit && uiHaloWidth < static_cast<UINT>(iEdge) && 0 == uiInvalidCount)
    {
        //The layout SHOULD have produced over-wide cells (HaloWidth < edge);
        //finding none means the scan itself is broken.
        ++uiError;
        LastProbem(_T("RadiusReject: no over-wide cells found though HaloWidth < kCacheIndexEdge"));
    }
    if (uiError > uiReported)
    {
        LastProbem(_T("RadiusReject verification failed"));
    }
    free(pTable);
    return uiError;
}
___REGIST_TEST(TestMGHaloRadiusReject, MG, TestMGHaloRadiusReject, HaloRadiusReject, _TEST_MULTIGPU);

/**
 * I8 gate 4: codimension-4 halo. The hyper-corner slots exist exactly when
 * all four directions are split (16 * HaloWidth^4 sites); their content
 * verifies after a refill like every other block. Run the real case with
 *   mpiexec -n 16 ./CLGTest TestMGHaloHyperCorner --mg-worker --gpu-grid 2,2,2,2 --device-per-node 1
 * (the yaml block pins HaloWidth 1 and a 4^4 global lattice, so each rank
 * holds a 2^4 sub-lattice; 16 is the MPI rank count from the metadata, ranks
 * oversubscribe the local device via --device-per-node 1).
 */
UINT TestMGHaloHyperCorner(CParameters&)
{
    UINT uiError = 0;
    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]);
    if (NULL == pGauge)
    {
        LastProbem(_T("No gauge field"));
        return 1;
    }

    SImprove1HaloCtx ctx;
    if (!Improve1BuildCtx(ctx))
    {
        appGeneral(_T("Halo count is 0 (nothing split), test trivially passes.\n"));
        return 0;
    }

    //Hyper-corner slot count vs layout: exactly 16 * W^4 when all four dirs
    //are split, exactly 0 otherwise.
    const UINT uiFaceEdgeCorner = _haloTotalFaceEdgeCornerSites(ctx.uiLocalL, ctx.uiGrid, ctx.uiVolume, ctx.uiHaloWidth);
    if (uiFaceEdgeCorner > ctx.uiHaloSites)
    {
        LastProbem(_T("Layout face+edge+corner count exceeds total halo count"));
        return 1;
    }
    const UINT uiHyper = ctx.uiHaloSites - uiFaceEdgeCorner;
    const UBOOL bAllSplit = (ctx.uiGrid[0] > 1 && ctx.uiGrid[1] > 1 && ctx.uiGrid[2] > 1 && ctx.uiGrid[3] > 1);
    const UINT uiW4 = ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth * ctx.uiHaloWidth;
    if (bAllSplit && uiHyper != 16 * uiW4)
    {
        ++uiError;
        CCString sProblem;
        sProblem.Format(_T("Hyper-corner slot count %u, expect %u (16 * HaloWidth^4)"), uiHyper, 16 * uiW4);
        appGeneral(_T("%s\n"), sProblem.c_str());
        LastProbem(sProblem);
    }
    if (!bAllSplit && 0 != uiHyper)
    {
        ++uiError;
        CCString sProblem;
        sProblem.Format(_T("Hyper-corner slots present (%u) though not all directions are split"), uiHyper);
        appGeneral(_T("%s\n"), sProblem.c_str());
        LastProbem(sProblem);
    }

    //Content: deterministic fill + refill; the check covers every block kind,
    //hyper-corner included.
    const UINT uiMatrixN = pGauge->MatrixN();
    const UINT uiElemPerSite = _HC_Dir * uiMatrixN * uiMatrixN * 2;
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));
    Improve1FillLocalDeterministic(pGauge->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
    appGetHaloManager()->Ensure(*pGauge->GetHaloBufferHandle(), ctx.uiHaloWidth);
    uiError += Improve1CheckHalo(pGauge->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("HyperCorner"));

    if (0 == uiError)
    {
        appGeneral(_T("HyperCorner verified: %u hyper-corner slots among %u halo slots all hold the expected global site data.\n"),
            uiHyper, ctx.uiHaloSites);
    }
    return uiError;
}
___REGIST_TEST(TestMGHaloHyperCorner, MG, TestMGHaloHyperCorner, HaloHyperCorner, _TEST_MULTIGPU);

/**
 * I8 gate 3: global-coordinate table with a global extent larger than 127
 * (the SCHAR range of the old narrowing path). The extent is asserted, so
 * shrinking the yaml lattice turns this test red instead of vacuous.
 */
UINT TestMGGlobalCoordinate32(CParameters&)
{
    UINT uiGrid[4];
    UINT uiOffset[4];
    UINT uiGlobL[4];
    UINT uiLocalL[4];
    Improve8FetchTopology(uiGrid, uiOffset, uiGlobL, uiLocalL);
    if (uiGlobL[0] <= 127 && uiGlobL[1] <= 127 && uiGlobL[2] <= 127 && uiGlobL[3] <= 127)
    {
        LastProbem(_T("GlobalCoordinate32 requires a global extent > 127 (yaml LatticeLength)"));
        return 1;
    }

    const UINT uiError = Improve8VerifyCoordinateTable(FALSE, _T("GlobalCoord32"));
    if (0 == uiError)
    {
        appGeneral(_T("Global coordinate table verified on a >127 lattice (global %ux%ux%ux%u).\n"),
            uiGlobL[0], uiGlobL[1], uiGlobL[2], uiGlobL[3]);
    }
    return uiError;
}
___REGIST_TEST(TestMGGlobalCoordinate32, MG, TestMGGlobalCoordinate32, GlobalCoordinate32, _TEST_MULTIGPU);

//=============================================================================
// END OF FILE
//=============================================================================
