//=============================================================================
// FILENAME : CBoundaryConditionProjectivePlaneSquare.cpp
// 
// DESCRIPTION:
// This is the periodic boundary condition
//
// REVISION:
//  [09/10/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CBoundaryConditionProjectivePlaneSquare)

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeEdgeProjectivePlaneBoundary(
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
    //realCoord.x = static_cast<SCHAR>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SCHAR>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SCHAR>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SCHAR>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    SCHAR signchange = 1;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //Improve-1 (3.8/I8d): the projective reflection and the BC sign apply
        //ONLY at the TRUE GLOBAL edge of the crossing direction; an internal
        //split boundary is plain neighbour data (halo redirect below), so the
        //transverse axis must NOT be mirrored on this rank.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bApplyNeg = (0 == iGlobalOffset);
        const UBOOL bApplyPos = (iGlobalOffset + iLocalLen == iGlobalLen);
#else
        const UBOOL bApplyNeg = TRUE;
        const UBOOL bApplyPos = TRUE;
#endif
        if (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            if (bApplyNeg)
            {
                signchange = signchange * bc.m_byData4[uiDir];
            }

            if (0 == uiDir)
            {
                if (bApplyNeg)
                {
                    realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                }
            }
            if (1 == uiDir)
            {
                if (bApplyNeg)
                {
                    realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
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
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            if (bApplyPos)
            {
                signchange = signchange * bc.m_byData4[uiDir];
            }

            if (0 == uiDir)
            {
                if (bApplyPos)
                {
                    realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                }
            }
            if (1 == uiDir)
            {
                if (bApplyPos)
                {
                    realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
                }
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
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;

#if _CLG_MULTI_GPU
    //Multi-GPU: cross-rank neighbour -> halo storage (see the Dirichlet bake).
    //The projective reflection stays active on the GLOBAL edge only.
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
        //Projective plane: x/y are reflection directions (boundary at the
        //global edge), z/t are periodic/anti-periodic (torus wrap between ranks).
        const BYTE byBoundaryMask[4] = { 1, 1, 0, 0 };
        UINT uiHaloSlot = 0;
        const EHaloRedirectResult eRedirect = _deviceHaloRedirectSiteBoundaryAware(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
            iRaw, iWrapped, iGlobalOffset, iGlobalLen, byBoundaryMask, uiHaloSlot);
        if (EHR_Halo == eRedirect)
        {
            pDeviceData[idxAll].m_uiSiteIndex = _DC_Volume + uiHaloSlot;
            pDeviceData[idxAll].m_byTag = _kGlue | (signchange < 0 ? _kDaggerOrOpposite : 0);
        }
        else if (EHR_Invalid == eRedirect)
        {
            //Improve-1 (3.8): a split-direction crossing beyond HaloWidth has
            //NO legal target -- never degrade it to the local handling above.
            pDeviceData[idxAll].m_uiSiteIndex = SIndex::_kInvalidSiteIndex;
            pDeviceData[idxAll].m_byTag = 0;
        }
    }
#endif
}

/**
* Nothing to write, just initial as 0
*/
//__global__ void _CLG_LAUNCH_BOUND
//_kernalBakeBondInfo_ProjectivePlane(BYTE* pDeviceData)
//{
//    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
//    for (UINT i = 0; i < _DC_Dir; ++i)
//    {
//        pDeviceData[idxAll * _DC_Dir + i] = 0;
//    }
//}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBoundGlueProjectivePlaneBoundary(
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
    //realCoord.x = static_cast<SCHAR>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SCHAR>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SCHAR>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SCHAR>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    //SCHAR signchange = 1;
    UBOOL bDaggerX = FALSE;
    UBOOL bDaggerY = FALSE;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
#if _CLG_MULTI_GPU
        //Improve-1 (3.8/I8d): the projective reflection / dagger applies ONLY
        //at the TRUE GLOBAL edge of the crossing direction; an internal split
        //boundary is plain neighbour data (halo redirect below), so the
        //transverse axis must NOT be mirrored on this rank.
        const INT iLocalLen = _constIntegers[ECI_Lx + uiDir];
        const INT iGlobalOffset = _constIntegers[ECI_GlobalOffsetX + uiDir];
        const INT iGlobalLen = _constIntegers[ECI_GlobalLx + uiDir];
        const UBOOL bApplyNeg = (0 == iGlobalOffset);
        const UBOOL bApplyPos = (iGlobalOffset + iLocalLen == iGlobalLen);
#else
        const UBOOL bApplyNeg = TRUE;
        const UBOOL bApplyPos = TRUE;
#endif
        if (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            //signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                if (bApplyNeg)
                {
                    realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                    bDaggerY = TRUE;
                }
            }
            if (1 == uiDir)
            {
                if (bApplyNeg)
                {
                    realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
                    bDaggerX = TRUE;
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
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            //signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                if (bApplyPos)
                {
                    realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                    bDaggerY = TRUE;
                }
            }
            if (1 == uiDir)
            {
                if (bApplyPos)
                {
                    realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
                    bDaggerX = TRUE;
                }
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

    //Note that, dagger X, means U_{-x} which is U(n-x)_x^+, same is Y
    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);

#if _CLG_MULTI_GPU
    //Multi-GPU: cross-rank neighbour -> halo storage (all links, no dagger).
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
        //Projective plane: x/y are reflection directions (boundary at the
        //global edge), z/t are periodic/anti-periodic (torus wrap between ranks).
        const BYTE byBoundaryMask[4] = { 1, 1, 0, 0 };
        UINT uiHaloSlot = 0;
        const EHaloRedirectResult eRedirect = _deviceHaloRedirectSiteBoundaryAware(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
            iRaw, iWrapped, iGlobalOffset, iGlobalLen, byBoundaryMask, uiHaloSlot);
        if (EHR_Halo == eRedirect)
        {
            for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
            {
                pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(_DC_Volume + uiHaloSlot);
                pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;
                //PP glue carries no sign/dagger semantics (x/y reflection only
                //applies on the global edge; cross-rank halo data is plain).
                pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = _kGlue;
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
            }
            return;
        }
    }
#endif

    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex);
        pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;

        if (0 == byDir && bDaggerX)
        {
            //1 uiSiteIndex to site4
            SSmallInt4 oldSite4 = __deviceSiteIndexToInt4Baking(uiSiteIndex);

            //2 shift site
            oldSite4.x = oldSite4.x - 1;

            //3 use the same logic to find new site index
            if (oldSite4.x < 0)
            {
                oldSite4.x = oldSite4.x + _DC_Lx;
                oldSite4.y = _DC_Ly - oldSite4.y - 1;
            }
            else if (oldSite4.x >= _DC_Lx)
            {
                oldSite4.x = oldSite4.x - _DC_Lx;
                oldSite4.y = _DC_Ly - oldSite4.y - 1;
            }

            const UINT uiSiteIndex2 = _deviceGetSiteIndex(oldSite4);
            pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex2);
            pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = _kDaggerOrOpposite;
        }

        if (1 == byDir && bDaggerY)
        {
            //1 uiSiteIndex to site4
            SSmallInt4 oldSite4 = __deviceSiteIndexToInt4Baking(uiSiteIndex);

            //2 shift site
            oldSite4.y = oldSite4.y - 1;

            //3 use the same logic to find new site index
            if (oldSite4.y < 0)
            {
                oldSite4.y = oldSite4.y + _DC_Ly;
                oldSite4.x = _DC_Lx - oldSite4.x - 1;
            }
            else if (oldSite4.y >= _DC_Ly)
            {
                oldSite4.y = oldSite4.y - _DC_Ly;
                oldSite4.x = _DC_Lx - oldSite4.x - 1;
            }

            const UINT uiSiteIndex2 = _deviceGetSiteIndex(oldSite4);
            pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex2);
            pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = _kDaggerOrOpposite;
        }

        //Bound should never have anti-periodic boundary condition?
        //pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
    }
}

#pragma endregion

CBoundaryConditionProjectivePlaneSquare::CBoundaryConditionProjectivePlaneSquare() : CBoundaryCondition()
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

void CBoundaryConditionProjectivePlaneSquare::BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
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

    _LAUNCH_KERNEL(_kernalBakeEdgeProjectivePlaneBoundary, blocks, threads, m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

//void CBoundaryConditionProjectivePlaneSquare::BakeBondInfo(const SSmallInt4*, BYTE* deviceTable, BYTE byFieldId) const
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
//    _LAUNCH_KERNEL(_kernalBakeBondInfo_ProjectivePlane, blocks, threads, deviceTable);
//}

void CBoundaryConditionProjectivePlaneSquare::BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const
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

    _LAUNCH_KERNEL(_kernalBakeBoundGlueProjectivePlaneBoundary, blocks, threads, m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
