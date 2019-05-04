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
    SIndex* pDeviceData, 
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 coord;
    coord.x = static_cast<SBYTE>(idxAll / mods.x);
    coord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y);
    coord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z);
    coord.w = static_cast<SBYTE>(idxAll % mods.z);

    SSmallInt4 realCoord = coord;
    realCoord.x -= CIndexData::kCacheIndexEdge;
    realCoord.y -= CIndexData::kCacheIndexEdge;
    realCoord.z -= CIndexData::kCacheIndexEdge;
    realCoord.w -= CIndexData::kCacheIndexEdge;

    SBYTE signchange = 1;
    for (UINT uiDir = 4 - _DC_Dir; uiDir < _DC_Dir; ++uiDir)
    {
        if (realCoord.m_byData4[uiDir] < 0)
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
        {
            realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            signchange = signchange * bc.m_byData4[uiDir];
        }
    }

    UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
}

#pragma endregion

CBoundaryConditionTorusSquare::CBoundaryConditionTorusSquare() : CBoundaryCondition()
{
    for (UINT i = 0; i < _kMaxFieldCount; ++i)
    {
        m_FieldBC[i].x = 1;
        m_FieldBC[i].y = 1;
        m_FieldBC[i].z = 1;
        m_FieldBC[i].w = -1;
    }
    m_FieldBC[0].w = 1;
    m_FieldBC[1].w = 1;
}

void CBoundaryConditionTorusSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
{
    assert(byFieldId < _kMaxFieldCount);
    m_FieldBC[byFieldId] = bc.m_sPeriodic;
}

void CBoundaryConditionTorusSquare::BakeEdgePoints(BYTE byFieldId, SIndex* deviceBuffer) const
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    uint3 biggerLatticeMod;

    UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    UINT threadPerSite = CIndexSquare::GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);
    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _kernalBakeEdgeTorusBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceBuffer, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
