//=============================================================================
// FILENAME : CBoundaryConditionTorusSquare.cpp
// 
// DESCRIPTION:
// This is the periodic boundary condition
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CBoundaryConditionTorusSquare)

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeEdgeTorusBoundary(
    SSmallInt4 bc, 
    const SSmallInt4* __restrict__ pMapping,
    SIndex* pDeviceData, 
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 realCoord(pMapping[idxAll]);

#if _CLG_MULTI_GPU
    //Raw (pre-wrap) neighbour coordinate, needed to tell a split-direction
    //out-of-lattice neighbour (lives on another rank -> halo) from an ordinary
    //periodic wrap within this rank.
    const INT iRaw[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
#endif

    SCHAR signchange = 1;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //The BC sign (e.g. antiperiodic-t, bc.w = -1) must toggle ONLY when the
        //neighbour crosses the TRUE GLOBAL lattice boundary, never an internal
        //split boundary. _constIntegers[ECI_Lx+dir] here is this rank's LOCAL
        //length, so an unconditional local wrap wrongly flips the sign on the
        //split-adjacent plane (bug signature: 2*kappa error at each internal
        //t-face). Gate each wrap on whether this rank sits at the global edge.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bNegEdge = (0 == iGlobalOffset);              // at global -edge
        const UBOOL bPosEdge = (iGlobalOffset + iLocalLen == iGlobalLen); // +edge
        while (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + iLocalLen;
            if (bNegEdge) { signchange = signchange * bc.m_byData4[uiDir]; }
        }
        while (realCoord.m_byData4[uiDir] >= iLocalLen)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - iLocalLen;
            if (bPosEdge) { signchange = signchange * bc.m_byData4[uiDir]; }
        }
#else
        while (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }

        while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }
#endif
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;

#if _CLG_MULTI_GPU
    //Redirect split-direction face neighbours to halo storage (Design B). The
    //wrapped local site index above stays as the DEFAULT (and is exactly the
    //correct value when this direction's process-grid neighbour is this rank
    //itself, i.e. the single-card self-exchange oracle). When the neighbour is a
    //different rank the halo buffer is filled by CHaloManager before any read.
    const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
    const UINT uiLocalL[4] = { static_cast<UINT>(_constIntegers[ECI_Lx]),
        static_cast<UINT>(_constIntegers[ECI_Ly]), static_cast<UINT>(_constIntegers[ECI_Lz]),
        static_cast<UINT>(_constIntegers[ECI_Lt]) };
    const INT iWrapped[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
    UINT uiHaloSlot = 0;
    const EHaloRedirectResult eRedirect = _deviceHaloRedirectSite(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
        iRaw, iWrapped, uiHaloSlot);
    if (EHR_Halo == eRedirect)
    {
        pDeviceData[idxAll].m_uiSiteIndex = _DC_Volume + uiHaloSlot;
        pDeviceData[idxAll].m_byTag |= _kGlue;
    }
    else if (EHR_Invalid == eRedirect)
    {
        //Improve-1 (3.8): a split-direction crossing beyond HaloWidth has NO
        //legal target -- never degrade it to the silent local wrap above.
        pDeviceData[idxAll].m_uiSiteIndex = SIndex::_kInvalidSiteIndex;
        pDeviceData[idxAll].m_byTag = 0;
    }
#endif
}

/**
* Nothing to write, just initial as 0
*/
//__global__ void _CLG_LAUNCH_BOUND
//_kernalBakeBondInfo_Torus(BYTE* pDeviceData)
//{
//    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
//    for (UINT i = 0; i < _DC_Dir; ++i)
//    {
//        pDeviceData[idxAll * _DC_Dir + i] = 0;
//    }
//}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBoundGlueTorusBoundary(
    SSmallInt4 bc,
    const SSmallInt4* __restrict__ pMapping,
    SIndex* pDeviceData,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 realCoord(pMapping[idxAll]);

#if _CLG_MULTI_GPU
    const INT iRaw[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
#endif

    SCHAR signchange = 1;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //See _kernalBakeEdgeTorusBoundary: only toggle the BC sign on a true
        //GLOBAL boundary crossing, not an internal split boundary.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bNegEdge = (0 == iGlobalOffset);
        const UBOOL bPosEdge = (iGlobalOffset + iLocalLen == iGlobalLen);
        while (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + iLocalLen;
            if (bNegEdge) { signchange = signchange * bc.m_byData4[uiDir]; }
        }
        while (realCoord.m_byData4[uiDir] >= iLocalLen)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - iLocalLen;
            if (bPosEdge) { signchange = signchange * bc.m_byData4[uiDir]; }
        }
#else
        while (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }

        while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }
#endif
    }

    UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);

#if _CLG_MULTI_GPU
    //Same split-direction face redirect as the site edge kernel: point the link
    //at the halo site's link block. A halo site owns _DC_Dir links laid out just
    //like a local site, so the link index is haloSite * Dir + byDir.
    const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
    const UINT uiLocalL[4] = { static_cast<UINT>(_constIntegers[ECI_Lx]),
        static_cast<UINT>(_constIntegers[ECI_Ly]), static_cast<UINT>(_constIntegers[ECI_Lz]),
        static_cast<UINT>(_constIntegers[ECI_Lt]) };
    const INT iWrapped[4] = { realCoord.x, realCoord.y, realCoord.z, realCoord.w };
    UINT uiHaloSlot = 0;
    const EHaloRedirectResult eGlue = _deviceHaloRedirectSite(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
        iRaw, iWrapped, uiHaloSlot);
    if (EHR_Halo == eGlue)
    {
        uiSiteIndex = _DC_Volume + uiHaloSlot;
    }
    else if (EHR_Invalid == eGlue)
    {
        //Improve-1 (3.8): a split-direction crossing beyond HaloWidth has NO
        //legal target -- never degrade it to the silent local wrap.
        uiSiteIndex = SIndex::_kInvalidSiteIndex;
    }
#endif

    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex);
        pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;

        //Bound should never have anti-periodic boundary condition?
        pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;

#if _CLG_MULTI_GPU
        if (EHR_Halo == eGlue)
        {
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag |= _kGlue;
        }
        else if (EHR_Invalid == eGlue)
        {
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = 0;
        }
#endif
    }
}

#pragma endregion

CBoundaryConditionTorusSquare::CBoundaryConditionTorusSquare() : CBoundaryCondition()
{
    for (UINT i = 0; i < kMaxFieldCount; ++i)
    {
        m_FieldBC[i].x = 1;
        m_FieldBC[i].y = 1;
        m_FieldBC[i].z = 1;
        m_FieldBC[i].w = -1;
    }
    m_FieldBC[0].w = 1;
    m_FieldBC[1].w = 1;
}

//void CBoundaryConditionTorusSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
//{
//    appAssert(byFieldId < kMaxFieldCount);
//    m_FieldBC[byFieldId] = bc.m_sPeriodic;
//}

void CBoundaryConditionTorusSquare::BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
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
    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _LAUNCH_KERNEL(_kernalBakeEdgeTorusBoundary, blocks, threads, m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

//void CBoundaryConditionTorusSquare::BakeBondInfo(const SSmallInt4*, BYTE* deviceTable, BYTE byFieldId) const
//{
//    uint4 biggerLattice;
//    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
//    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
//
//    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
//    const UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
//    dim3 threads(threadPerSite, 1, 1);
//    dim3 blocks(uiVolumn / threadPerSite, 1, 1);
//
//    _LAUNCH_KERNEL(_kernalBakeBondInfo_Torus, blocks, threads, deviceTable);
//}

void CBoundaryConditionTorusSquare::BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
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
    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _LAUNCH_KERNEL(_kernalBakeBoundGlueTorusBoundary, blocks, threads, m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
