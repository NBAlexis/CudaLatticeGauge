//=============================================================================
// FILENAME : CBoundaryConditionTorusSquare.cpp
// 
// DESCRIPTION:
// This is the periodic boundary condition
// 
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CBoundaryConditionPeriodicAndDirichletSquare)

#pragma region device functions

//This is XOR...
//static __device__ __inline__ BYTE _deviceToggleBitInverse(BYTE value, BYTE toggle)
//{
//    //return (value & (~toggle)) | ((~value) & toggle);
//    return value ^ toggle;
//}

/**
* This function is not using
*/
//static __device__ __inline__ UBOOL _deviceOnEdge(BYTE regionId, BYTE muLeft)
//{
//    return 0 != ((regionId & muLeft) ^ ((regionId >> 4) & muLeft));
//}

#pragma endregion

#pragma region kernels

/**
* I have to sort this rule because the tree-improve gauge not work with Dirichlet now
* pMapping is a big-index to site4 mapping
* This is for sites
*/
__global__ void _CLG_LAUNCH_BOUND
_kernalBakeEdgePeriodicDirichletBoundary(
    SSmallInt4 bc,
    SIndex* pDeviceData,
    const SSmallInt4* __restrict__ pMapping,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    
    SSmallInt4 realCoord(pMapping[idxAll]);
#if _CLG_MULTI_GPU
    const INT iRaw[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
#endif
    //realCoord.x = static_cast<SCHAR>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SCHAR>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SCHAR>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SCHAR>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;
    
    //SSmallInt4 orig = realCoord;

    SCHAR signchange = 1;
    BYTE byRegionId = 0;
    UBOOL bBoundary = FALSE;
    UBOOL bOutside = FALSE;
    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //Improve-1 (3.8/I8d): the Dirichlet plane/region/outside semantics and
        //the periodic BC sign apply ONLY at the TRUE GLOBAL lattice edge of the
        //crossing direction (see _kernalBakeEdgeTorusBoundary). At an internal
        //split boundary the neighbour is plain halo data (redirect below), so
        //none of the boundary semantics may fire on this rank.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bApplyNeg = (0 == iGlobalOffset);
        const UBOOL bApplyPos = (iGlobalOffset + iLocalLen == iGlobalLen);
#else
        const UBOOL bApplyNeg = TRUE;
        const UBOOL bApplyPos = TRUE;
#endif
        if (realCoord.m_byData4[uiDir] <= 0)
        {
            //printf("-- coord[uiDir]=%d --\n", static_cast<INT>(realCoord.m_byData4[uiDir]));
            UBOOL bPassEdge = FALSE;
            if (realCoord.m_byData4[uiDir] < 0)
            {
                bPassEdge = TRUE;
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];

                if (0 == bc.m_byData4[uiDir])
                {
                    realCoord.m_byData4[uiDir] = 0; // -1;
                    if (bApplyNeg)
                    {
                        bOutside = TRUE;
                    }
                }
                else
                {
                    while (realCoord.m_byData4[uiDir] < 0)
                    {
                        realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
                    }
                }
            }

            if (0 == bc.m_byData4[uiDir])
            {
                if (bApplyNeg)
                {
                    bBoundary = TRUE;
                    //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << uiDir);
                    byRegionId = byRegionId ^ (1 << uiDir);
                }
            }
            else if (bPassEdge && bApplyNeg)
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }

        }
        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1)
        {
            UBOOL bPassEdge = TRUE;
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];

            if (0 == bc.m_byData4[uiDir])
            {
                realCoord.m_byData4[uiDir] = 0; // _constIntegers[ECI_Lx + uiDir];
                if (bApplyPos)
                {
                    bOutside = TRUE;
                }
            }
            else
            {
                while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
                {
                    realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
                }
            }

            if (0 == bc.m_byData4[uiDir])
            {
                if (bApplyPos)
                {
                    bBoundary = TRUE;
                    //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << (uiDir + 4));
                    byRegionId = byRegionId ^ (1 << (uiDir + 4));
                }
            }
            else if (bPassEdge && bApplyPos) //realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }
        }
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
    //if (signchange < 0)
    //{
    //    printf("sign change %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
    //}
    pDeviceData[idxAll].m_byReginId = byRegionId;

    if (bBoundary)
    {
        appAssert(0 != byRegionId);
        //printf("We have dirichlet bc %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
        pDeviceData[idxAll].m_byTag |= _kDirichlet;
    }
    if (bOutside)
    {
        pDeviceData[idxAll].m_byTag |= _kDirichlet;
        pDeviceData[idxAll].m_byTag |= _kOutside;
    }

#if _CLG_MULTI_GPU
    //Multi-GPU: a neighbour crossing a process-grid split (and not on the
    //global edge) lives in halo storage; the halo buffer is filled by
    //CHaloManager with the neighbour rank's boundary data before any read.
    //A GLOBAL-edge crossing keeps the Dirichlet/periodic semantics above.
    {
        const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
        const UINT uiLocalL[4] = { static_cast<UINT>(_constIntegers[ECI_Lx]),
            static_cast<UINT>(_constIntegers[ECI_Ly]), static_cast<UINT>(_constIntegers[ECI_Lz]),
            static_cast<UINT>(_constIntegers[ECI_Lt]) };
        const INT iWrapped[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
        const INT iGlobalOffset[4] = { static_cast<INT>(_constIntegers[ECI_GlobalOffsetX]),
            static_cast<INT>(_constIntegers[ECI_GlobalOffsetY]), static_cast<INT>(_constIntegers[ECI_GlobalOffsetZ]),
            static_cast<INT>(_constIntegers[ECI_GlobalOffsetT]) };
        const INT iGlobalLen[4] = { static_cast<INT>(_constIntegers[ECI_GlobalLx]),
            static_cast<INT>(_constIntegers[ECI_GlobalLy]), static_cast<INT>(_constIntegers[ECI_GlobalLz]),
            static_cast<INT>(_constIntegers[ECI_GlobalLt]) };
        const BYTE byBoundaryMask[4] = { static_cast<BYTE>((0 == bc.x) ? 1 : 0), static_cast<BYTE>((0 == bc.y) ? 1 : 0),
            static_cast<BYTE>((0 == bc.z) ? 1 : 0), static_cast<BYTE>((0 == bc.w) ? 1 : 0) };
        UINT uiHaloSlot = 0;
        const EHaloRedirectResult eRedirect = _deviceHaloRedirectSiteBoundaryAware(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
            iRaw, iWrapped, iGlobalOffset, iGlobalLen, byBoundaryMask, uiHaloSlot);
        if (EHR_Halo == eRedirect)
        {
            pDeviceData[idxAll].m_uiSiteIndex = _DC_Volume + uiHaloSlot;
            pDeviceData[idxAll].m_byTag = _kGlue | (signchange < 0 ? _kDaggerOrOpposite : 0);
            pDeviceData[idxAll].m_byReginId = 0;
        }
        else if (EHR_Invalid == eRedirect)
        {
            //Improve-1 (3.8): a split-direction crossing beyond HaloWidth has
            //NO legal target -- never degrade it to the local handling above.
            pDeviceData[idxAll].m_uiSiteIndex = SIndex::_kInvalidSiteIndex;
            pDeviceData[idxAll].m_byTag = 0;
            pDeviceData[idxAll].m_byReginId = 0;
        }
    }
#endif
}

//__global__ void _CLG_LAUNCH_BOUND
//_kernalBakeBondInfoPeriodicDirichletBoundary(
//    SSmallInt4 bc,
//    BYTE* pDeviceData,
//    const SSmallInt4* __restrict__ pMapping,
//    uint3 mods)
//{
//    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
//
//    SSmallInt4 realCoord(pMapping[idxAll]);
//
//    UBOOL bDirich = FALSE;
//    BYTE byZeroCount = 0;
//    BYTE byZeroDir = 0;
//    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
//    {
//        if (realCoord.m_byData4[uiDir] <= 0 && 0 == bc.m_byData4[uiDir])
//        {
//            if (realCoord.m_byData4[uiDir] < 0)
//            {
//                bDirich = TRUE;
//            }
//            else
//            {
//                byZeroCount++;
//                byZeroDir = uiDir;
//            }
//        }
//        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1 && 0 == bc.m_byData4[uiDir])
//        {
//            bDirich = TRUE;
//        }
//    }
//
//    if (bDirich || byZeroCount > 1)
//    {
//        //all bonds are Dirichlet
//        for (UINT i = 0; i < _DC_Dir; ++i)
//        {
//            pDeviceData[idxAll * _DC_Dir + i] = _kDirichlet;
//        }
//        //printf("idx = %d,%d,%d,%d, dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
//    }
//    else if (1 == byZeroCount)
//    {
//        //except the direction of the zero index, all others are Dirichlet
//        for (UINT i = 0; i < _DC_Dir; ++i)
//        {
//            if (i != byZeroDir)
//            {
//                pDeviceData[idxAll * _DC_Dir + i] = _kDirichlet;
//            }
//            else
//            {
//                pDeviceData[idxAll * _DC_Dir + i] = 0;
//            }
//        }
//        //printf("idx = %d,%d,%d,%d, half dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
//    }
//    else
//    {
//        //not Dirichlet
//        for (UINT i = 0; i < _DC_Dir; ++i)
//        {
//            pDeviceData[idxAll * _DC_Dir + i] = 0;
//        }
//        //printf("idx = %d,%d,%d,%d, not dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
//    }
//}

/**
* I have to sort this rule because the tree-improve gauge not work with Dirichlet now
* pMapping is a big-index to site4 mapping
* This is for links
* 
* For example, if X is Dirichlet, then:
* U_{mu}(n_x<0) are all passedge Dirichlet
* U_{mu!=x}(n_x=0) are Dirichlet
* U_{mu}(n_x>=Lx) are all passedge Dirichlet
*/
__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBoundGlue_DBC(
    SSmallInt4 bc,
    SIndex* pDeviceData,
    const SSmallInt4* __restrict__ pMapping,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;

    SSmallInt4 realCoord(pMapping[idxAll]);
#if _CLG_MULTI_GPU
    const INT iRaw[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
#endif
    //realCoord.x = static_cast<SCHAR>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SCHAR>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SCHAR>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SCHAR>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    //SSmallInt4 orig = realCoord;

    SCHAR signchange = 1;
    BYTE byRegionId = 0;
    UBOOL bBoundary = FALSE;
    UBOOL bPassDirichletEdge = FALSE;
    BYTE tagBoundDir = 0;
    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //Improve-1 (3.8/I8d): Dirichlet plane/region/passthrough semantics and
        //the periodic BC sign apply ONLY at the TRUE GLOBAL edge of the
        //crossing direction; an internal split boundary is halo business (the
        //redirect below), never a physical boundary on this rank.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bApplyNeg = (0 == iGlobalOffset);
        const UBOOL bApplyPos = (iGlobalOffset + iLocalLen == iGlobalLen);
#else
        const UBOOL bApplyNeg = TRUE;
        const UBOOL bApplyPos = TRUE;
#endif
        if (realCoord.m_byData4[uiDir] <= 0)
        {
            UBOOL bPassEdge = FALSE;
            //printf("-- coord[uiDir]=%d --\n", static_cast<INT>(realCoord.m_byData4[uiDir]));
            if (realCoord.m_byData4[uiDir] < 0)
            {
                bPassEdge = TRUE;
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
                if (0 == bc.m_byData4[uiDir])
                {
                    if (bApplyNeg)
                    {
                        bPassDirichletEdge = TRUE;
                    }
                    realCoord.m_byData4[uiDir] = -1;
                }
                else
                {
                    while (realCoord.m_byData4[uiDir] < 0)
                    {
                        realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
                    }
                }
            }

            if (0 == bc.m_byData4[uiDir])
            {
                if (bApplyNeg)
                {
                    bBoundary = TRUE;
                    //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << uiDir);
                    byRegionId = byRegionId ^ (1 << uiDir);
                    tagBoundDir |= (1 << uiDir);
                }
            }
            else if (bPassEdge && bApplyNeg)
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }

        }
        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1)
        {
            //if x > Lx - 1, then x >= Lx
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            if (0 == bc.m_byData4[uiDir])
            {
                if (bApplyPos)
                {
                    bPassDirichletEdge = TRUE;
                    signchange = signchange * bc.m_byData4[uiDir];
                    byRegionId = byRegionId ^ (1 << uiDir);
                    tagBoundDir |= (1 << uiDir);
                }
                realCoord.m_byData4[uiDir] = _constIntegers[ECI_Lx + uiDir];
            }
            else
            {
                while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
                {
                    realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
                }
            }
        }
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);

#if _CLG_MULTI_GPU
    //Multi-GPU: cross-rank neighbour -> halo storage (see edge kernel). All
    //links of this boundary cell then point at the same halo slot and carry no
    //Dirichlet/sign semantics (the halo data is the neighbour's plain field).
    {
        const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
        const UINT uiLocalL[4] = { static_cast<UINT>(_constIntegers[ECI_Lx]),
            static_cast<UINT>(_constIntegers[ECI_Ly]), static_cast<UINT>(_constIntegers[ECI_Lz]),
            static_cast<UINT>(_constIntegers[ECI_Lt]) };
        const INT iWrapped[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
        const INT iGlobalOffset[4] = { static_cast<INT>(_constIntegers[ECI_GlobalOffsetX]),
            static_cast<INT>(_constIntegers[ECI_GlobalOffsetY]), static_cast<INT>(_constIntegers[ECI_GlobalOffsetZ]),
            static_cast<INT>(_constIntegers[ECI_GlobalOffsetT]) };
        const INT iGlobalLen[4] = { static_cast<INT>(_constIntegers[ECI_GlobalLx]),
            static_cast<INT>(_constIntegers[ECI_GlobalLy]), static_cast<INT>(_constIntegers[ECI_GlobalLz]),
            static_cast<INT>(_constIntegers[ECI_GlobalLt]) };
        const BYTE byBoundaryMask[4] = { static_cast<BYTE>((0 == bc.x) ? 1 : 0), static_cast<BYTE>((0 == bc.y) ? 1 : 0),
            static_cast<BYTE>((0 == bc.z) ? 1 : 0), static_cast<BYTE>((0 == bc.w) ? 1 : 0) };
        UINT uiHaloSlot = 0;
        const EHaloRedirectResult eRedirect = _deviceHaloRedirectSiteBoundaryAware(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
            iRaw, iWrapped, iGlobalOffset, iGlobalLen, byBoundaryMask, uiHaloSlot);
        if (EHR_Halo == eRedirect)
        {
            for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
            {
                pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(_DC_Volume + uiHaloSlot);
                pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;
                pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = _kGlue | (signchange < 0 ? _kDaggerOrOpposite : 0);
                pDeviceData[idxAll * _DC_Dir + byDir].m_byReginId = 0;
            }
            return;
        }
        if (EHR_Invalid == eRedirect)
        {
            //Improve-1 (3.8): a split-direction crossing beyond HaloWidth has
            //NO legal target -- never degrade it to the local handling below.
            for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
            {
                pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(SIndex::_kInvalidSiteIndex);
                pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;
                pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = 0;
                pDeviceData[idxAll * _DC_Dir + byDir].m_byReginId = 0;
            }
            return;
        }
    }
#endif

    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex);
        pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;

        //Bound should never have anti-periodic boundary condition?
        pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
        //if (signchange < 0)
        //{
        //    printf("sign change\n");
        //}
        pDeviceData[idxAll * _DC_Dir + byDir].m_byReginId = byRegionId;

        if (bPassDirichletEdge)
        {
            appAssert(0 != byRegionId);
            //printf("We have dirichlet bc %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag |= (_kDirichlet | _kOutside);
        }
        else if (bBoundary)
        {
            appAssert(0 != byRegionId);

            //if mu are Dirichlet, and if me is mu on nu-boundary with nu != mu, I am Dirichlet
            UBOOL bReallyOnBoundary = FALSE;
            for (BYTE byDir2 = 0; byDir2 < _DC_Dir; ++byDir2)
            {
                if (byDir2 != byDir)
                {
                    bReallyOnBoundary = bReallyOnBoundary || (tagBoundDir & (1 << byDir2));
                }
            }

            if (bReallyOnBoundary)
            {
                pDeviceData[idxAll * _DC_Dir + byDir].m_byTag |= _kDirichlet;
            }
        }
    }
}

#pragma endregion

CBoundaryConditionPeriodicAndDirichletSquare::CBoundaryConditionPeriodicAndDirichletSquare() : CBoundaryCondition()
{
    for (UINT i = 0; i < kMaxFieldCount; ++i)
    {
        m_FieldBC[i].x = 0;
        m_FieldBC[i].y = 0;
        m_FieldBC[i].z = 1;
        m_FieldBC[i].w = -1;
    }
    m_FieldBC[0].w = 1;
    m_FieldBC[1].w = 1;
}

//void CBoundaryConditionPeriodicAndDirichletSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
//{
//    appAssert(byFieldId < kMaxFieldCount);
//    m_FieldBC[byFieldId] = bc.m_sPeriodic;
//}

void CBoundaryConditionPeriodicAndDirichletSquare::BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    uint3 biggerLatticeMod;

    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    const UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);

    //appGeneral(_T("block=%d, %d, %d, thread= %d, %d, %d\n"), blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _LAUNCH_KERNEL(_kernalBakeEdgePeriodicDirichletBoundary, blocks, threads, m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
}

void CBoundaryConditionPeriodicAndDirichletSquare::BakeRegionTable(UINT* deviceTable) const
{
    UINT regionTable[256];

    for (UINT i = 0; i < 256; ++i)
    {
        regionTable[i] = 0;
        if (0 != (i & byXLeft))
        {
            regionTable[i] = 0;
        }
        else if (0 != (i & byYLeft))
        {
            regionTable[i] = 1;
        }
        else if (0 != (i & byZLeft))
        {
            regionTable[i] = 2;
        }
        else if (0 != (i & byTLeft))
        {
            regionTable[i] = 3;
        }
        else if (0 != (i & byXRight))
        {
            regionTable[i] = 4;
        }
        else if (0 != (i & byYRight))
        {
            regionTable[i] = 5;
        }
        else if (0 != (i & byZRight))
        {
            regionTable[i] = 6;
        }
        else if (0 != (i & byTRight))
        {
            regionTable[i] = 7;
        }
    }

    checkCudaErrors(cudaMemcpy(deviceTable, regionTable, sizeof(UINT) * 256, cudaMemcpyHostToDevice));
}

//void CBoundaryConditionPeriodicAndDirichletSquare::BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceBuffer, BYTE byFieldId) const
//{
//    uint4 biggerLattice;
//    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
//    uint3 biggerLatticeMod;
//
//    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
//    const UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
//    dim3 threads(threadPerSite, 1, 1);
//    dim3 blocks(uiVolumn / threadPerSite, 1, 1);
//
//    //appGeneral(_T("block=%d, %d, %d, thread= %d, %d, %d\n"), blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);
//
//    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
//    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
//    biggerLatticeMod.z = biggerLattice.w;
//
//    _LAUNCH_KERNEL(_kernalBakeBondInfoPeriodicDirichletBoundary, blocks, threads, m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
//}

void CBoundaryConditionPeriodicAndDirichletSquare::BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    uint3 biggerLatticeMod;

    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    const UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);

    //appGeneral(_T("block=%d, %d, %d, thread= %d, %d, %d\n"), blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _LAUNCH_KERNEL(_kernalBakeBoundGlue_DBC, blocks, threads, m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
