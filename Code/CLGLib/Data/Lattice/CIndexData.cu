//=============================================================================
// FILENAME : CIndexData.cu
// 
// DESCRIPTION:
// This is for testing the index
//
// REVISION:
//  [05/04/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE


#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelPrintWalkingTable()
{
    UINT bigIndex = threadIdx.x + blockDim.x * blockIdx.x;

    SSmallInt4 coord = __idx->_deviceBigIndexToInt4(bigIndex);

    UINT neightbouridx1 = __idx->m_pWalkingTable[bigIndex * _DC_Dir * 2 + 0];
    SSmallInt4 mappingIdx1 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][neightbouridx1].m_uiSiteIndex);
    UINT neightbouridx2 = __idx->m_pWalkingTable[bigIndex * _DC_Dir * 2 + 2];
    SSmallInt4 mappingIdx2 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][neightbouridx2].m_uiSiteIndex);
    UINT neightbouridx3 = __idx->m_pWalkingTable[bigIndex * _DC_Dir * 2 + 5];
    SSmallInt4 mappingIdx3 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][neightbouridx3].m_uiSiteIndex);
    UINT neightbouridx4 = __idx->m_pWalkingTable[bigIndex * _DC_Dir * 2 + 7];
    SSmallInt4 mappingIdx4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][neightbouridx4].m_uiSiteIndex);

    printf("coord: (%d,%d,%d,%d) - \n%d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d)\n",
        coord.x, coord.y, coord.z, coord.w, 
        0, mappingIdx1.x, mappingIdx1.y, mappingIdx1.z, mappingIdx1.w,
        2, mappingIdx2.x, mappingIdx2.y, mappingIdx2.z, mappingIdx2.w,
        5, mappingIdx3.x, mappingIdx3.y, mappingIdx3.z, mappingIdx3.w,
        7, mappingIdx4.x, mappingIdx4.y, mappingIdx4.z, mappingIdx4.w
    );
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDebugPlaqutteTable()
{
    intokernal;

    if (123 == uiSiteIndex)
    {
        SSmallInt4 myself = __deviceSiteIndexToInt4(uiSiteIndex);
        printf("me = %d, %d, %d, %d\n", myself.x, myself.y, myself.z, myself.w);

        UINT uiListIdx = uiSiteIndex
            * __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx]
            * __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];

        for (UINT i = 0; 
            i < __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] 
            * __idx->m_pSmallData[CIndexData::kPlaqLengthIdx]; ++i)
        {
            SSmallInt4 coord = __deviceSiteIndexToInt4(__idx->m_pPlaqutteCache[uiListIdx + i].m_uiSiteIndex);
            printf("%s%s(%d,%d,%d,%d)_(%d)\n",
                (__idx->m_pPlaqutteCache[uiListIdx + i].m_byTag & _kDaggerOrOpposite ? "-" : "+"),
                (__idx->m_pPlaqutteCache[uiListIdx + i].IsDirichlet() ? "D" : ""),
                coord.x, coord.y, coord.z, coord.w,
                __idx->m_pPlaqutteCache[uiListIdx + i].m_byDir);
        }
    }
}

#pragma endregion


UINT inline _GetDecompose(UINT volumn)
{
    TArray<UINT> factors = _getFactors(volumn);
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    UINT maxThreadPerBlock = deviceConstraints[0];

    UINT uiMax = 1;
    for (INT i = 0; i < factors.Num(); ++i)
    {
        if (factors[i] <= maxThreadPerBlock && factors[i] > uiMax)
        {
            uiMax = factors[i];
        }
    }

    return uiMax;
}

void CIndexData::DebugPrintWalkingTable()
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    //uint3 biggerLatticeMod;

    UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    UINT threadPerSite = _GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);

    appGeneral(_T("decomp %d x %d = %d\n"), threads.x, blocks.x, uiVolumn);

    _kernelPrintWalkingTable << <blocks, threads >> > ();
}

void CIndexData::DebugPlaqutteTable()
{
    preparethread;
    _kernelDebugPlaqutteTable << <block, threads >> > ();
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================