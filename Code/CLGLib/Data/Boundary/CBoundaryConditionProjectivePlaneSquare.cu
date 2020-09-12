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
    //realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    SBYTE signchange = 1;
    for (UINT uiDir = 4 - _DC_Dir; uiDir < _DC_Dir; ++uiDir)
    {
        if (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
            }
            if (1 == uiDir)
            {
                realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
            }
        }
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
            }
            if (1 == uiDir)
            {
                realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
            }
        }
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
}

/**
* Nothing to write, just initial as 0
*/
__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBondInfo_ProjectivePlane(BYTE* pDeviceData)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    for (UINT i = 0; i < _DC_Dir; ++i)
    {
        pDeviceData[idxAll * _DC_Dir + i] = 0;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBoundGlueProjectivePlaneBoundary(
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

    //SBYTE signchange = 1;
    UBOOL bDaggerX = FALSE;
    UBOOL bDaggerY = FALSE;
    for (UINT uiDir = 4 - _DC_Dir; uiDir < _DC_Dir; ++uiDir)
    {
        if (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            //signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                bDaggerY = TRUE;
            }
            if (1 == uiDir)
            {
                realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
                bDaggerX = TRUE;
            }
        }
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            //signchange = signchange * bc.m_byData4[uiDir];

            if (0 == uiDir)
            {
                realCoord.m_byData4[1] = _DC_Ly - realCoord.m_byData4[1] - 1;
                bDaggerY = TRUE;
            }
            if (1 == uiDir)
            {
                realCoord.m_byData4[0] = _DC_Lx - realCoord.m_byData4[0] - 1;
                bDaggerX = TRUE;
            }
        }
    }

    //Note that, dagger X, means U_{-x} which is U(n-x)_x^+, same is Y
    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        pDeviceData[idxAll * _DC_Dir + byDir] = SIndex(uiSiteIndex);
        pDeviceData[idxAll * _DC_Dir + byDir].m_byDir = byDir;

        if (0 == byDir && bDaggerX)
        {
            //1 uiSiteIndex to site4
            SSmallInt4 oldSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

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
            SSmallInt4 oldSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

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

void CBoundaryConditionProjectivePlaneSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
{
    assert(byFieldId < kMaxFieldCount);
    m_FieldBC[byFieldId] = bc.m_sPeriodic;
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

    _kernalBakeEdgeProjectivePlaneBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

void CBoundaryConditionProjectivePlaneSquare::BakeBondInfo(const SSmallInt4*, BYTE* deviceTable) const
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;

    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    const UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);

    _kernalBakeBondInfo_ProjectivePlane << <blocks, threads >> > (deviceTable);
}

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

    _kernalBakeBoundGlueProjectivePlaneBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceMappingTable, deviceBuffer, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
