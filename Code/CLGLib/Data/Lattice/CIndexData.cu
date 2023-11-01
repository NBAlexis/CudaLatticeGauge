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

//void CIndexData::DebugPrintWalkingTable()
//{
//    const UINT uiSize = (_HC_Lx + 2 * CIndexData::kCacheIndexEdge)
//                      * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
//                      * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
//                      * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);
//    UINT* tb = (UINT*)malloc(sizeof(UINT) * uiSize * 2 * _HC_Dir);
//    SIndex* stb = (SIndex*)malloc(sizeof(SIndex) * uiSize);
//    checkCudaErrors(cudaMemcpy(tb, appGetLattice()->m_pIndexCache->m_pWalkingTable, sizeof(UINT) * uiSize * 2 * _HC_Dir, cudaMemcpyDeviceToHost));
//    checkCudaErrors(cudaMemcpy(stb, appGetLattice()->m_pIndexCache->m_pIndexPositionToSIndex[1], sizeof(SIndex) * uiSize, cudaMemcpyDeviceToHost));
//    appPushLogDate(FALSE);
//    for (UINT uiBigIdx = 0; uiBigIdx < uiSize; ++uiBigIdx)
//    {
//        const SSmallInt4 coord = _hostBigIndexToInt4(uiBigIdx);
//
//        const UINT neightbouridx1 = tb[uiBigIdx * _HC_Dir * 2 + 0];
//        const SSmallInt4 mappingIdx1 = __hostSiteIndexToInt4(stb[neightbouridx1].m_uiSiteIndex);
//        const UINT neightbouridx2 = tb[uiBigIdx * _HC_Dir * 2 + 2];
//        const SSmallInt4 mappingIdx2 = __hostSiteIndexToInt4(stb[neightbouridx2].m_uiSiteIndex);
//        const UINT neightbouridx3 = tb[uiBigIdx * _HC_Dir * 2 + 5];
//        const SSmallInt4 mappingIdx3 = __hostSiteIndexToInt4(stb[neightbouridx3].m_uiSiteIndex);
//        const UINT neightbouridx4 = tb[uiBigIdx * _HC_Dir * 2 + 7];
//        const SSmallInt4 mappingIdx4 = __hostSiteIndexToInt4(stb[neightbouridx4].m_uiSiteIndex);
//
//        appGeneral(_T("coord: (%d,%d,%d,%d)  -  %d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d) %d=(%d,%d,%d,%d)\n"),
//            coord.x, coord.y, coord.z, coord.w,
//            0, mappingIdx1.x, mappingIdx1.y, mappingIdx1.z, mappingIdx1.w,
//            2, mappingIdx2.x, mappingIdx2.y, mappingIdx2.z, mappingIdx2.w,
//            5, mappingIdx3.x, mappingIdx3.y, mappingIdx3.z, mappingIdx3.w,
//            7, mappingIdx4.x, mappingIdx4.y, mappingIdx4.z, mappingIdx4.w
//        );
//    }
//    appPopLogDate();
//
//    appSafeFree(tb);
//    appSafeFree(stb);
//}

void CIndexData::DebugPlaqutteTable()
{
    appPushLogDate(FALSE);
    const UINT uiPlaqLength = 4;
    const UINT uiPlaqPerSite = _HC_Dim * (_HC_Dim - 1) / 2;
    const UINT uiSiteCount = _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt;
    SIndex* cache = (SIndex*)malloc(sizeof(SIndex) * uiSiteCount * uiPlaqLength * uiPlaqPerSite);
    checkCudaErrors(cudaMemcpy(cache, appGetLattice()->m_pIndexCache->m_pPlaqutteCache, sizeof(SIndex) * uiSiteCount * uiPlaqLength * uiPlaqPerSite, cudaMemcpyDeviceToHost));

    for (UINT uiSite = 0; uiSite < uiSiteCount; ++uiSite)
    {
        const SSmallInt4 myself = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T("me = %d, %d, %d, %d\n"), myself.x, myself.y, myself.z, myself.w);

        const UINT uiListIdx = uiSite * uiPlaqPerSite * uiPlaqLength;

        for (UINT i = 0; i < uiPlaqPerSite; ++i)
        {
            for (UINT j = 0; j < uiPlaqLength; ++j)
            {
                const UINT uiIdx = i * uiPlaqLength + j;
                const SIndex& idx = cache[uiListIdx + uiIdx];
                const SSmallInt4 coord = __hostSiteIndexToInt4(idx.m_uiSiteIndex);

                appGeneral(_T("%s%s(%d,%d,%d,%d)_(%d)  "),
                    (idx.m_byTag & _kDaggerOrOpposite ? "-" : "+"),
                    (idx.m_byTag & _kDirichlet ? "D" : ""),
                    coord.x, coord.y, coord.z, coord.w,
                    idx.m_byDir);
            }
            appGeneral(_T("\n"));
        }
    }
    appSafeFree(cache);
    appPopLogDate();
}

void CIndexData::DebugPlaqutteTable(const SSmallInt4& sSite)
{
    const UINT uiPlaqLength = 4;
    const UINT uiPlaqPerSite = _HC_Dim * (_HC_Dim - 1) / 2;
    const UINT uiSiteCount = _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt;

    SIndex* cache = (SIndex*)malloc(sizeof(SIndex) * uiSiteCount * uiPlaqLength * uiPlaqPerSite);
    checkCudaErrors(cudaMemcpy(cache, appGetLattice()->m_pIndexCache->m_pStappleCache, sizeof(SIndex) * uiSiteCount * uiPlaqLength * uiPlaqPerSite, cudaMemcpyDeviceToHost));
    UINT uiSite = _hostGetSiteIndex(sSite);
    appGeneral(_T("me = %d, %d, %d, %d\n"), sSite.x, sSite.y, sSite.z, sSite.w);

    const UINT uiListIdx = uiSite * uiPlaqPerSite * uiPlaqLength;

    for (UINT i = 0; i < uiPlaqPerSite; ++i)
    {
        for (UINT j = 0; j < uiPlaqLength; ++j)
        {
            const UINT uiIdx = uiListIdx + i * uiPlaqLength + j;
            const SIndex& idx = cache[uiListIdx + uiIdx];
            const SSmallInt4 coord = __hostSiteIndexToInt4(idx.m_uiSiteIndex);

            appGeneral(_T("%s%s(%d,%d,%d,%d)_(%d)  "),
                (idx.m_byTag & _kDaggerOrOpposite ? "-" : "+"),
                (idx.m_byTag & _kDirichlet ? "D" : ""),
                coord.x, coord.y, coord.z, coord.w,
                idx.m_byDir);
        }
        appGeneral(_T("\n"));
    }
    appSafeFree(cache);
}

void CIndexData::DebugStapleTable()
{
    appPushLogDate(FALSE);
    const UINT uiPlaqLength = 4;
    const UINT uiPlaqPerLink = 2 * (_HC_Dim - 1);
    const UINT uiSiteCount = _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt;
    SIndex* cache = (SIndex*)malloc(sizeof(SIndex) * uiSiteCount * (uiPlaqLength - 1) * uiPlaqPerLink);
    checkCudaErrors(cudaMemcpy(cache, appGetLattice()->m_pIndexCache->m_pStappleCache, 
        sizeof(SIndex) * uiSiteCount * (uiPlaqLength - 1) * uiPlaqPerLink, cudaMemcpyDeviceToHost));

    for (UINT uiSite = 0; uiSite < uiSiteCount; ++uiSite)
    {
        const SSmallInt4 myself = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T("me = %d, %d, %d, %d\n"), myself.x, myself.y, myself.z, myself.w);

        const UINT uiListIdx = uiSite * (uiPlaqLength - 1) * uiPlaqPerLink;

        for (UINT i = 0; i < uiPlaqPerLink; ++i)
        {
            for (UINT j = 0; j < uiPlaqLength - 1; ++j)
            {
                const UINT uiIdx = i * (uiPlaqLength - 1) + j;
                const SIndex& idx = cache[uiListIdx + uiIdx];
                const SSmallInt4 coord = __hostSiteIndexToInt4(idx.m_uiSiteIndex);

                appGeneral(_T("%s%s(%d,%d,%d,%d)_(%d)  "),
                    (idx.m_byTag & _kDaggerOrOpposite ? "-" : "+"),
                    (idx.m_byTag & _kDirichlet ? "D" : ""),
                    coord.x, coord.y, coord.z, coord.w,
                    idx.m_byDir);
            }
            appGeneral(_T("\n"));
        }
    }
    appSafeFree(cache);
    appPopLogDate();
}

void CIndexData::DebugEdgeMapping(BYTE byFieldId, const SSmallInt4& xyzt)
{
    const UINT uiSize = (_HC_Lx + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);

    SIndex* tb = (SIndex*)malloc(sizeof(SIndex) * uiSize);

    checkCudaErrors(cudaMemcpy(tb, 
        appGetLattice()->m_pIndexCache->m_pIndexPositionToSIndex[byFieldId], 
        sizeof(SIndex) * uiSize, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiBigIdx = 0; uiBigIdx < uiSize; ++uiBigIdx)
    {
        const SSmallInt4 coord = _hostBigIndexToInt4(uiBigIdx);
        INT x = static_cast<INT>(coord.x);
        INT y = static_cast<INT>(coord.y);
        INT z = static_cast<INT>(coord.z);
        INT t = static_cast<INT>(coord.w);
        if (x < 0 || x >= static_cast<INT>(_HC_Lx)
         || y < 0 || y >= static_cast<INT>(_HC_Ly)
         || z < 0 || z >= static_cast<INT>(_HC_Lz)
         || t < 0 || t >= static_cast<INT>(_HC_Lt))
        {
            if ( (-1 == xyzt.x || coord.x == xyzt.x)
              && (-1 == xyzt.y || coord.y == xyzt.y)
              && (-1 == xyzt.z || coord.z == xyzt.z)
              && (-1 == xyzt.w || coord.w == xyzt.w)
                )
            {
                const SSmallInt4 mappingIdx = __hostSiteIndexToInt4(tb[uiBigIdx].m_uiSiteIndex);
                appGeneral(_T("coord: (%d,%d,%d,%d)=(%d,%d,%d,%d)\n"),
                    coord.x, coord.y, coord.z, coord.w,
                    mappingIdx.x, mappingIdx.y, mappingIdx.z, mappingIdx.w);
            }
        }
    }
    appPopLogDate();

    appSafeFree(tb);
}

void CIndexData::DebugEdgeGlue(BYTE byFieldId, const SSmallInt4& xyzt)
{
    const UINT uiSize = (_HC_Lx + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);

    SIndex* tb = (SIndex*)malloc(sizeof(SIndex) * uiSize * _HC_Dir);
    checkCudaErrors(cudaMemcpy(tb,
        appGetLattice()->m_pIndexCache->m_pIndexLinkToSIndex[byFieldId],
        sizeof(SIndex) * uiSize * _HC_Dir, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiBigIdx = 0; uiBigIdx < uiSize; ++uiBigIdx)
    {
        const SSmallInt4 coord = _hostBigIndexToInt4(uiBigIdx);
        INT x = static_cast<INT>(coord.x);
        INT y = static_cast<INT>(coord.y);
        INT z = static_cast<INT>(coord.z);
        INT t = static_cast<INT>(coord.w);
        if (x < 0 || x >= static_cast<INT>(_HC_Lx)
            || y < 0 || y >= static_cast<INT>(_HC_Ly)
            || z < 0 || z >= static_cast<INT>(_HC_Lz)
            || t < 0 || t >= static_cast<INT>(_HC_Lt))
        {
            if ((-1 == xyzt.x || coord.x == xyzt.x)
                && (-1 == xyzt.y || coord.y == xyzt.y)
                && (-1 == xyzt.z || coord.z == xyzt.z)
                && (-1 == xyzt.w || coord.w == xyzt.w)
                )
            {
                for (BYTE byDir = 0; byDir < _HC_Dir; ++byDir)
                {
                    const SIndex& idx = tb[uiBigIdx * _HC_Dir + byDir];
                    const SSmallInt4 mappingIdx = __hostSiteIndexToInt4(idx.m_uiSiteIndex);
                    appGeneral(_T("(%d,%d,%d,%d)_%d=(%d,%d,%d,%d)_%d%s "),
                        coord.x, coord.y, coord.z, coord.w, byDir,
                        mappingIdx.x, mappingIdx.y, mappingIdx.z, mappingIdx.w,
                        idx.m_byDir, (idx.m_byTag & _kDaggerOrOpposite) ? _T("+") : _T(""));
                }
                appGeneral(_T("\n"));
            }
        }
    }
    appPopLogDate();

    appSafeFree(tb);
}

void CIndexData::DebugLinkDirichletOrDagger(BYTE byFieldId)
{
    const UINT uiSize = (_HC_Lx + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);

    SIndex* tb = (SIndex*)malloc(sizeof(SIndex) * uiSize * _HC_Dir);
    checkCudaErrors(cudaMemcpy(tb,
        appGetLattice()->m_pIndexCache->m_pIndexLinkToSIndex[byFieldId],
        sizeof(SIndex) * uiSize * _HC_Dir, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiBigIdx = 0; uiBigIdx < uiSize; ++uiBigIdx)
    {
        const SSmallInt4 coord = _hostBigIndexToInt4(uiBigIdx);
        INT x = static_cast<INT>(coord.x);
        INT y = static_cast<INT>(coord.y);
        INT z = static_cast<INT>(coord.z);
        INT t = static_cast<INT>(coord.w);
        if (x >= 0 && x < _HC_Lxi
            && y >= 0 && y < _HC_Lyi
            && z >= 0 && z < _HC_Lzi
            && t >= 0 && t < _HC_Lti)
        {
            for (BYTE byDir = 0; byDir < _HC_Dir; ++byDir)
            {
                const SIndex& idx = tb[uiBigIdx * _HC_Dir + byDir];
                if (0 != (_kDirichlet & idx.m_byTag))
                {
                    appGeneral(_T("(%d, %d, %d, %d)_%d is dirichlet\n"), x, y, z, t, byDir);
                }
                if (0 != (_kDaggerOrOpposite & idx.m_byTag))
                {
                    appGeneral(_T("(%d, %d, %d, %d)_%d is reversed\n"), x, y, z, t, byDir);
                }
            }
        }
    }
    appPopLogDate();

    appSafeFree(tb);
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