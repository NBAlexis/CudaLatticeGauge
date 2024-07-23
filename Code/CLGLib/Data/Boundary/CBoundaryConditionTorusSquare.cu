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
    //realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    //UBOOL bDebug = FALSE;
    //if (realCoord.m_byData4[3] < 0 || realCoord.m_byData4[3] > 1)
    //{
    //    bDebug = TRUE;
    //}
    //SSmallInt4 old(realCoord);

    SBYTE signchange = 1;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
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
    }

    //if (bDebug)
    //{
    //    printf("%d %d %d %d to %d %d %d %d\n", old.x, old.y, old.z, old.w, realCoord.x, realCoord.y, realCoord.z, realCoord.w);
    //}

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
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
    //realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    SBYTE signchange = 1;
    for (UINT uiDir = 0; uiDir < 4; ++uiDir)
    {
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
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex);
        pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;

        //Bound should never have anti-periodic boundary condition?
        pDeviceData[idxAll * _DC_Dir + byDir].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
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
//    assert(byFieldId < kMaxFieldCount);
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

    _kernalBakeEdgeTorusBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
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
//    _kernalBakeBondInfo_Torus << <blocks, threads >> > (deviceTable);
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

    _kernalBakeBoundGlueTorusBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
