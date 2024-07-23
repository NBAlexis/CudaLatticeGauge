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

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeEdgePeriodicDirichletBoundary(
    SSmallInt4 bc,
    SIndex* pDeviceData,
    const SSmallInt4* __restrict__ pMapping,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    
    SSmallInt4 realCoord(pMapping[idxAll]);
    //realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;
    
    //SSmallInt4 orig = realCoord;

    SBYTE signchange = 1;
    BYTE byRegionId = 0;
    UBOOL bBoundary = FALSE;
    
    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
    {
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
                    realCoord.m_byData4[uiDir] = 0;
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
                bBoundary = TRUE;
                //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << uiDir);
                byRegionId = byRegionId ^ (1 << uiDir);
            }
            else if (bPassEdge)
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }

        }
        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1)
        {
            UBOOL bPassEdge = FALSE;
            if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                bPassEdge = TRUE;
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];

                if (0 == bc.m_byData4[uiDir])
                {
                    realCoord.m_byData4[uiDir] = _constIntegers[ECI_Lx + uiDir] - 1;
                }
                else
                {
                    while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
                    {
                        realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
                    }
                }
            }

            if (0 == bc.m_byData4[uiDir])
            {
                bBoundary = TRUE;
                //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << (uiDir + 4));
                byRegionId = byRegionId ^ (1 << (uiDir + 4));
            }
            else if (bPassEdge) //realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
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
    //    printf("sign change\n");
    //}
    pDeviceData[idxAll].m_byReginId = byRegionId;

    if (bBoundary)
    {
        assert(0 != byRegionId);
        //printf("We have dirichlet bc %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
        pDeviceData[idxAll].m_byTag |= _kDirichlet;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBondInfoPeriodicDirichletBoundary(
    SSmallInt4 bc,
    BYTE* pDeviceData,
    const SSmallInt4* __restrict__ pMapping,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;

    SSmallInt4 realCoord(pMapping[idxAll]);

    UBOOL bDirich = FALSE;
    BYTE byZeroCount = 0;
    BYTE byZeroDir = 0;
    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
    {
        if (realCoord.m_byData4[uiDir] <= 0 && 0 == bc.m_byData4[uiDir])
        {
            if (realCoord.m_byData4[uiDir] < 0)
            {
                bDirich = TRUE;
            }
            else
            {
                byZeroCount++;
                byZeroDir = uiDir;
            }
        }
        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1 && 0 == bc.m_byData4[uiDir])
        {
            bDirich = TRUE;
        }
    }

    if (bDirich || byZeroCount > 1)
    {
        //all bonds are Dirichlet
        for (UINT i = 0; i < _DC_Dir; ++i)
        {
            pDeviceData[idxAll * _DC_Dir + i] = _kDirichlet;
        }
        //printf("idx = %d,%d,%d,%d, dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
    }
    else if (1 == byZeroCount)
    {
        //except the direction of the zero index, all others are Dirichlet
        for (UINT i = 0; i < _DC_Dir; ++i)
        {
            if (i != byZeroDir)
            {
                pDeviceData[idxAll * _DC_Dir + i] = _kDirichlet;
            }
            else
            {
                pDeviceData[idxAll * _DC_Dir + i] = 0;
            }
        }
        //printf("idx = %d,%d,%d,%d, half dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
    }
    else
    {
        //not Dirichlet
        for (UINT i = 0; i < _DC_Dir; ++i)
        {
            pDeviceData[idxAll * _DC_Dir + i] = 0;
        }
        //printf("idx = %d,%d,%d,%d, not dirichlet\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernalBakeBoundGlue_DBC(
    SSmallInt4 bc,
    SIndex* pDeviceData,
    const SSmallInt4* __restrict__ pMapping,
    uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;

    SSmallInt4 realCoord(pMapping[idxAll]);
    //realCoord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    //realCoord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    //realCoord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    //realCoord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    //SSmallInt4 orig = realCoord;

    SBYTE signchange = 1;
    BYTE byRegionId = 0;
    UBOOL bBoundary = FALSE;
    UBOOL bPassDirichletEdge = FALSE;
    BYTE tagBoundDir = 0;
    for (BYTE uiDir = 0; uiDir < 4; ++uiDir)
    {
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
                    bPassDirichletEdge = TRUE;
                    realCoord.m_byData4[uiDir] = 0;
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
                bBoundary = TRUE;
                //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << uiDir);
                byRegionId = byRegionId ^ (1 << uiDir);
                tagBoundDir |= (1 << uiDir);
            }
            else if (bPassEdge)
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }

        }
        else if (realCoord.m_byData4[uiDir] > _constIntegers[ECI_Lx + uiDir] - 1)
        {
            UBOOL bPassEdge = FALSE;
            if (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                bPassEdge = TRUE;
                realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
                if (0 == bc.m_byData4[uiDir])
                {
                    bPassDirichletEdge = TRUE;
                    realCoord.m_byData4[uiDir] = _constIntegers[ECI_Lx + uiDir] - 1;
                }
                else
                {
                    while (realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
                    {
                        realCoord.m_byData4[uiDir] = realCoord.m_byData4[uiDir] - _constIntegers[ECI_Lx + uiDir];
                    }
                }
            }

            if (0 == bc.m_byData4[uiDir])
            {
                bBoundary = TRUE;
                //byRegionId = _deviceToggleBitInverse(byRegionId, 1 << (uiDir + 4));
                byRegionId = byRegionId ^ (1 << (uiDir + 4));
                tagBoundDir |= (1 << uiDir);
            }
            else if (bPassEdge) //realCoord.m_byData4[uiDir] >= _constIntegers[ECI_Lx + uiDir])
            {
                //printf("bc=%d\n", static_cast<INT>(bc.m_byData4[uiDir]));
                signchange = signchange * bc.m_byData4[uiDir];
            }
        }
    }

    const UINT uiSiteIndex = _deviceGetSiteIndex(realCoord);

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
            assert(0 != byRegionId);
            //printf("We have dirichlet bc %d %d %d %d\n", orig.x, orig.y, orig.z, orig.w);
            pDeviceData[idxAll * _DC_Dir + byDir].m_byTag |= _kDirichlet;
        }
        else if (bBoundary)
        {
            assert(0 != byRegionId);

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
//    assert(byFieldId < kMaxFieldCount);
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

    _kernalBakeEdgePeriodicDirichletBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
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

void CBoundaryConditionPeriodicAndDirichletSquare::BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceBuffer, BYTE byFieldId) const
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

    _kernalBakeBondInfoPeriodicDirichletBoundary << <blocks, threads >> > (m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
}

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

    _kernalBakeBoundGlue_DBC << <blocks, threads >> > (m_FieldBC[byFieldId], deviceBuffer, deviceMappingTable, biggerLatticeMod);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
