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

__CLGIMPLEMENT_CLASS(CBoundaryConditionPeriodicAndDirichletSquare)

#pragma region device functions

static __device__ __inline__ BYTE _deviceToggleBitInverse(BYTE value, BYTE toggle)
{
    return (value & (~toggle)) | ((~value) & toggle);
}

/**
* This function is not using
*/
//static __device__ __inline__ UBOOL _deviceOnEdge(BYTE regionId, BYTE muLeft)
//{
//    return 0 != ((regionId & muLeft) ^ ((regionId >> 4) & muLeft));
//}

#pragma endregion

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeEdgePeriodicDirichletBoundary(
    SSmallInt4 bc,
    SIndex* pDeviceData,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    
    SSmallInt4 realCoord;
    realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;
    
    SSmallInt4 orig = realCoord;

    SBYTE signchange = 1;
    BYTE byRegionId = 0;
    for (BYTE uiDir = static_cast<BYTE>(4 - _DC_Dir); uiDir < _DC_Dir; ++uiDir)
    {
        if (realCoord.m_byData4[uiDir] <= 0)
        {
            if (realCoord.m_byData4[uiDir] < 0)
            {
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] + _constIntegers[ECI_Lx + uiDir];
            }
            
            if (0 == bc.m_byData4[uiDir])
            {
                byRegionId = _deviceToggleBitInverse(byRegionId, 1 << uiDir);
            }
            else if (realCoord.m_byData4[uiDir] < 0)
            {
                signchange = signchange * bc.m_byData4[uiDir];
            }
            
        }
        else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir] - 1)
        {
            if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
            }
            
            if (0 == bc.m_byData4[uiDir])
            {
                byRegionId = _deviceToggleBitInverse(byRegionId, 1 << (uiDir + 4));
            }
            else if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                signchange = signchange * bc.m_byData4[uiDir];
            }
        }
    }

    UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);
    pDeviceData[idxAll] = SIndex(uiSiteIndex);
    pDeviceData[idxAll].m_byTag = signchange < 0 ? _kDaggerOrOpposite : 0;
    pDeviceData[idxAll].m_byReginId = byRegionId;
    if (0 != byRegionId)
    {
        //printf("We have dirichlet bc %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
        pDeviceData[idxAll].m_byTag |= _kDirichlet;
    }
}

#pragma endregion

CBoundaryConditionPeriodicAndDirichletSquare::CBoundaryConditionPeriodicAndDirichletSquare() : CBoundaryCondition()
{
    for (UINT i = 0; i < _kMaxFieldCount; ++i)
    {
        m_FieldBC[i].x = 0;
        m_FieldBC[i].y = 0;
        m_FieldBC[i].z = 1;
        m_FieldBC[i].w = -1;
    }
    m_FieldBC[0].w = 1;
    m_FieldBC[1].w = 1;
}

void CBoundaryConditionPeriodicAndDirichletSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
{
    assert(byFieldId < _kMaxFieldCount);
    m_FieldBC[byFieldId] = bc.m_sPeriodic;
}

void CBoundaryConditionPeriodicAndDirichletSquare::BakeEdgePoints(BYTE byFieldId, SIndex* deviceBuffer) const
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

    //appGeneral(_T("block=%d, %d, %d, thread= %d, %d, %d\n"), blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    _kernalBakeEdgePeriodicDirichletBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceBuffer, biggerLatticeMod);
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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================