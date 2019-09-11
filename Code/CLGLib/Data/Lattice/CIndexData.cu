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


__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateLinkCount(
    INT* res,
    const BYTE* __restrict__ pBoundInfo,
    const UINT* __restrict__ pSmallData
)
{
    intokernalOnlyInt4;
    const UINT uiBigIdx = _deviceGetBigIndex(sSite4, pSmallData);

    INT uiCount = 0;
    for (BYTE byDir = 0; byDir < _DC_Dir; ++byDir)
    {
        if (0 == (pBoundInfo[uiBigIdx * _DC_Dir + byDir] & _kDirichlet))
        {
            ++uiCount;
        }
    }
    atomicAdd(res, uiCount);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateSiteCount(
    INT* res,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallData
)
{
    intokernalOnlyInt4;
    const UINT uiBigIdx = _deviceGetBigIndex(sSite4, pSmallData);

    if (!pMappingTable[uiBigIdx].IsDirichlet())
    {
        atomicAdd(&res[0], 1);
        if (0 == sSite4.w)
        {
            atomicAdd(&res[1], 1);
        }
    }
}

#pragma endregion


UINT inline _GetDecompose(UINT volumn)
{
    TArray<UINT> factors = _getFactors(volumn);
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    const UINT maxThreadPerBlock = deviceConstraints[0];

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

    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    const UINT threadPerSite = _GetDecompose(uiVolumn);
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

void CIndex::CalculateSiteCount(class CIndexData* pData) const
{
    INT hostres[2] = { 0, 0 };
    INT* deviceRes = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceRes, sizeof(INT) * 2));

    preparethread;

    assert(NULL != pData->m_pSmallData);

    for (BYTE i = 1; i < kMaxFieldCount; ++i)
    {
        if (NULL != pData->m_pIndexPositionToSIndex[i])
        {
            hostres[0] = 0;
            hostres[1] = 0;
            checkCudaErrors(cudaMemcpy(deviceRes, hostres, sizeof(INT) * 2, cudaMemcpyHostToDevice));

            _kernelCalculateSiteCount << <block, threads >> > (deviceRes, pData->m_pIndexPositionToSIndex[i], pData->m_pSmallData);
            checkCudaErrors(cudaMemcpy(hostres, deviceRes, sizeof(INT) * 2, cudaMemcpyDeviceToHost));
            pData->m_uiSiteNumber[i] = static_cast<UINT>(hostres[0]);

            if (1 == i)
            {
                pData->m_uiSiteXYZT = static_cast<UINT>(hostres[0]);
                pData->m_uiSiteXYZ = static_cast<UINT>(hostres[1]);
                appGeneral(_T("============== Real Volume = %d ============\n"), pData->m_uiSiteXYZT);
                appGeneral(_T("============== Real Spatial Volume = %d ============\n"), pData->m_uiSiteXYZ);
            }
        }
    }

    hostres[0] = 0;
    hostres[1] = 0;
    
    checkCudaErrors(cudaMemcpy(deviceRes, hostres, sizeof(INT) * 2, cudaMemcpyHostToDevice));
    _kernelCalculateLinkCount << <block, threads >> > (deviceRes, pData->m_pBondInfoTable, pData->m_pSmallData);
    checkCudaErrors(cudaMemcpy(hostres, deviceRes, sizeof(INT) * 2, cudaMemcpyDeviceToHost));
    pData->m_uiLinkNumber = static_cast<UINT>(hostres[0]);
    checkCudaErrors(cudaFree(deviceRes));
    appGeneral(_T("============== Real Link Count = %d ============\n"), pData->m_uiLinkNumber);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================