//=============================================================================
// FILENAME : CMeasure.cpp
// 
// DESCRIPTION:
// Some common functions
//
// REVISION:
//  [09/26/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial as zero
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero_XYPlane(Real* pBuffer)
{
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    pBuffer[_ixy] = F(0.0);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero_XYPlaneC(CLGComplex* pBuffer)
{
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    pBuffer[_ixy] = _zeroc;
}

/**
* Average over z and t
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverZT_XYPlane(Real* pBuffer)
{
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    pBuffer[_ixy] = pBuffer[_ixy] / (_DC_Lz * _DC_Lt);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverZT_XYPlaneC(CLGComplex* pBuffer)
{
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    pBuffer[_ixy].x = pBuffer[_ixy].x / (_DC_Lz * _DC_Lt);
    pBuffer[_ixy].y = pBuffer[_ixy].y / (_DC_Lz * _DC_Lt);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_R_C(
    const CLGComplex* __restrict__ jgsXY,
    SSmallInt4 sCenter, 
    UINT uiMax, 
    BYTE byFieldId,
    CLGComplex* result,
    UINT* pCount)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    SBYTE uiX = static_cast<SBYTE>(uiXY / _DC_Ly);
    SBYTE uiY = static_cast<SBYTE>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
             + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = uiX;
    sSite4.y = uiY;
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (NULL != pCount)
        {
            atomicAdd(&pCount[uiC], 1);
        }
        atomicAdd(&result[uiC].x, jgsXY[uiXY].x);
        atomicAdd(&result[uiC].y, jgsXY[uiXY].y);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_RAverage_C(const UINT* __restrict__ pCount, CLGComplex* pValue)
{
    const UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pValue[uiIdx].x = pValue[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
        pValue[uiIdx].y = pValue[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_R_R(
    const Real* __restrict__ jgsXY,
    SSmallInt4 sCenter,
    UINT uiMax,
    BYTE byFieldId,
    Real* result,
    UINT* pCount)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    INT uiX = static_cast<INT>(uiXY / _DC_Ly);
    INT uiY = static_cast<INT>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
             + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = static_cast<SBYTE>(uiX);
    sSite4.y = static_cast<SBYTE>(uiY);
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (NULL != pCount)
        {
            atomicAdd(&pCount[uiC], 1);
        }
        atomicAdd(&result[uiC], jgsXY[uiXY]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_RAverage_R(const UINT* __restrict__ pCount, Real* pValue)
{
    const UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pValue[uiIdx] = pValue[uiIdx] / static_cast<Real>(pCount[uiIdx]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialDist(UINT* pCount, Real* pValue, CLGComplex* pValueC)
{
    if (NULL != pCount)
    {
        pCount[threadIdx.x] = 0;
    }
    if (NULL != pValue)
    {
        pValue[threadIdx.x] = F(0.0);
    }
    if (NULL != pValueC)
    {
        pValueC[threadIdx.x] = _zeroc;
    }
}

#pragma endregion

void CMeasure::FillDataWithR_R(
    TArray<Real>& arrData,
    TArray<Real>& arrInner,
    TArray<Real>& arrFull,
    TArray<UINT>& arrR,
    Real* hostData,
    UINT* hostR,
    UINT uiConfig,
    UINT uiMaxR,
    UINT uiEdgeR,
    UBOOL bFillR)
{
    Real fAverageJGInner = F(0.0);
    Real fAverageJGAll = F(0.0);
    UINT uiInnerPointsAll = 0;
    UINT uiInnerPointsInner = 0;
    if (0 == uiConfig)
    {
        assert(!bFillR || 0 == arrR.Num());
        assert(0 == arrData.Num());

        for (UINT uiL = 0; uiL <= uiMaxR; ++uiL)
        {
            if (hostR[uiL] > 0)
            {
                if (bFillR)
                {
                    arrR.AddItem(uiL);
                }
                
                arrData.AddItem(hostData[uiL]);

                uiInnerPointsAll += hostR[uiL];
                fAverageJGAll += hostR[uiL] * hostData[uiL];
                if (uiL < uiEdgeR)
                {
                    uiInnerPointsInner += hostR[uiL];
                    fAverageJGInner += hostR[uiL] * hostData[uiL];
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < arrR.Num(); ++i)
        {
            assert(hostR[arrR[i]] > 0);
            arrData.AddItem(hostData[arrR[i]]);

            uiInnerPointsAll += hostR[arrR[i]];
            fAverageJGAll += hostR[arrR[i]] * hostData[arrR[i]];
            if (arrR[i] < uiEdgeR)
            {
                uiInnerPointsInner += hostR[arrR[i]];
                fAverageJGInner += hostR[arrR[i]] * hostData[arrR[i]];
            }
        }
    }

    if (uiInnerPointsAll > 0)
    {
        fAverageJGAll = fAverageJGAll / uiInnerPointsAll;
    }
    if (uiInnerPointsInner > 0)
    {
        fAverageJGInner = fAverageJGInner / uiInnerPointsInner;
    }
    arrFull.AddItem(fAverageJGAll);
    arrInner.AddItem(fAverageJGInner);
}

void CMeasure::FillDataWithR_C(
    TArray<CLGComplex>& arrData,
    TArray<CLGComplex>& arrInner,
    TArray<CLGComplex>& arrFull,
    TArray<UINT>& arrR,
    CLGComplex* hostData,
    UINT* hostR,
    UINT uiConfig,
    UINT uiMaxR,
    UINT uiEdgeR,
    UBOOL bFillR)
{
    CLGComplex cAverageJGInner = _zeroc;
    CLGComplex cAverageJGAll = _zeroc;
    UINT uiInnerPointsAll = 0;
    UINT uiInnerPointsInner = 0;
    if (0 == uiConfig)
    {
        assert(!bFillR || 0 == arrR.Num());
        assert(0 == arrData.Num());

        for (UINT uiL = 0; uiL <= uiMaxR; ++uiL)
        {
            if (hostR[uiL] > 0)
            {
                if (bFillR)
                {
                    arrR.AddItem(uiL);
                }

                arrData.AddItem(hostData[uiL]);

                uiInnerPointsAll += hostR[uiL];
                cAverageJGAll.x += hostR[uiL] * hostData[uiL].x;
                cAverageJGAll.y += hostR[uiL] * hostData[uiL].y;
                if (uiL < uiEdgeR)
                {
                    uiInnerPointsInner += hostR[uiL];
                    cAverageJGInner.x += hostR[uiL] * hostData[uiL].x;
                    cAverageJGInner.y += hostR[uiL] * hostData[uiL].y;
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < arrR.Num(); ++i)
        {
            assert(hostR[arrR[i]] > 0);
            arrData.AddItem(hostData[arrR[i]]);

            uiInnerPointsAll += hostR[arrR[i]];
            cAverageJGAll.x += hostR[arrR[i]] * hostData[arrR[i]].x;
            cAverageJGAll.y += hostR[arrR[i]] * hostData[arrR[i]].y;
            if (arrR[i] < uiEdgeR)
            {
                uiInnerPointsInner += hostR[arrR[i]];
                cAverageJGInner.x += hostR[arrR[i]] * hostData[arrR[i]].x;
                cAverageJGInner.y += hostR[arrR[i]] * hostData[arrR[i]].y;
            }
        }
    }

    if (uiInnerPointsAll > 0)
    {
        cAverageJGAll.x = cAverageJGAll.x / uiInnerPointsAll;
        cAverageJGAll.y = cAverageJGAll.y / uiInnerPointsAll;
    }
    if (uiInnerPointsInner > 0)
    {
        cAverageJGInner.x = cAverageJGInner.x / uiInnerPointsInner;
        cAverageJGInner.y = cAverageJGInner.y / uiInnerPointsInner;
    }
    arrFull.AddItem(cAverageJGAll);
    arrInner.AddItem(cAverageJGInner);
}


void CMeasure::ReportDistributionXY_R(UINT uiConfig, const TArray<Real>& arrayRes)
{
    assert(uiConfig * (_HC_Lx - 1) * (_HC_Ly - 1)
        == static_cast<UINT>(arrayRes.Num()));

    TArray<Real> tmpjgs;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < uiConfig; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                const UINT idx = k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j;

                if (0 == k)
                {
                    tmpjgs.AddItem(arrayRes[idx]);
                }
                else
                {
                    tmpjgs[i * (_HC_Lx - 1) + j] += arrayRes[idx];
                }

                if (0 == j)
                {
                    appGeneral(_T("%2.12f"), arrayRes[idx]);
                }
                else
                {
                    appGeneral(_T(", %2.12f"), arrayRes[idx]);
                }
            }
            appGeneral(_T("}, "));
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n -------------------- Average -------------------------\n\n"));

    for (UINT i = 0; i < _HC_Ly - 1; ++i)
    {
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("(x=%d,y=%d)%2.8f,   "),
                j + 1, i + 1,
                tmpjgs[i * (_HC_Lx - 1) + j] / uiConfig);
        }
        appGeneral(_T("\n"));
    }
}

void CMeasure::ReportDistributeWithR_R(UINT uiConf, UINT uiR, const TArray<Real>& arrayData)
{
    assert(uiConf * uiR == static_cast<UINT>(arrayData.GetCount()));
    appGeneral(_T("{\n"));
    for (UINT conf = 0; conf < uiConf; ++conf)
    {
        for (UINT r = 0; r < uiR; ++r)
        {
            if (0 == r)
            {
                appGeneral(_T("{ %2.12f"), arrayData[uiR * conf + r]);
            }
            else
            {
                appGeneral(_T(", %2.12f"), arrayData[uiR * conf + r]);
            }
        }

        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));
}

void CMeasure::ReportDistributionXY_C(UINT uiConfig, const TArray<CLGComplex>& arrayRes)
{
    assert(uiConfig * (_HC_Lx - 1) * (_HC_Ly - 1)
        == static_cast<UINT>(arrayRes.Num()));

    TArray<CLGComplex> tmpjgs;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < uiConfig; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                const UINT idx = k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j;

                if (0 == k)
                {
                    tmpjgs.AddItem(arrayRes[idx]);
                }
                else
                {
                    tmpjgs[i * (_HC_Lx - 1) + j] = _cuCaddf(tmpjgs[i * (_HC_Lx - 1) + j], arrayRes[idx]);
                }

                if (0 == j)
                {
                    appGeneral(_T("%2.12f %s %2.12f I"), 
                        arrayRes[idx].x, 
                        arrayRes[idx].y > 0 ? _T("+") : _T("-"), 
                        appAbs(arrayRes[idx].y));
                }
                else
                {
                    appGeneral(_T(", %2.12f %s %2.12f I"),
                        arrayRes[idx].x,
                        arrayRes[idx].y > 0 ? _T("+") : _T("-"),
                        appAbs(arrayRes[idx].y));
                }
            }
            appGeneral(_T("},\n  "));
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n -------------------- Average -------------------------\n\n"));

    for (UINT i = 0; i < _HC_Ly - 1; ++i)
    {
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("(x=%d,y=%d)%2.8f + %2.8f I,   "), 
                j + 1, i + 1, 
                tmpjgs[i * (_HC_Lx - 1) + j].x / uiConfig, 
                tmpjgs[i * (_HC_Lx - 1) + j].y / uiConfig);
        }
        appGeneral(_T("\n"));
    }
}

/**
* array[x, y] = array[x, y] / (lz * lt)
*/
void CMeasure::_AverageXYPlane(Real* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1);
    const dim3 threads(_HC_DecompLx, 1, 1);
    _kernelAverageOverZT_XYPlane << <block, threads >> > (pDeviceRes);
}

/**
* array[x, y] = 0
*/
void CMeasure::_ZeroXYPlane(Real* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1);
    const dim3 threads(_HC_DecompLx, 1, 1);
    _kernelInitialZero_XYPlane << <block, threads >> > (pDeviceRes);
}

/**
* array[x, y] = array[x, y] / (lz * lt)
*/
void CMeasure::_AverageXYPlaneC(CLGComplex* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1);
    const dim3 threads(_HC_DecompLx, 1, 1);
    _kernelAverageOverZT_XYPlaneC << <block, threads >> > (pDeviceRes);
}

void CMeasure::_ZeroXYPlaneC(CLGComplex* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1); 
    const dim3 threads(_HC_DecompLx, 1, 1);
    _kernelInitialZero_XYPlaneC << <block, threads >> > (pDeviceRes);
}

void CMeasure::XYDataToRdistri_R(
    const Real* __restrict__ source,
    UINT* count,
    Real* result,
    UINT uiMaxR,
    UBOOL bCalculateCounter,
    BYTE byFieldId)
{
    const dim3 block2(_HC_DecompX, 1, 1);
    const dim3 threads2(_HC_DecompLx, 1, 1);
    const dim3 block3(1, 1, 1);
    const dim3 threads3(uiMaxR + 1, 1, 1);

    _kernelInitialDist << <block3, threads3 >> > (
        bCalculateCounter ? count : NULL,
        result,
        NULL
        );

    _kernelXY_To_R_R << <block2, threads2 >> > (
        source,
        CCommonData::m_sCenter,
        uiMaxR,
        byFieldId,
        result,
        bCalculateCounter ? count : NULL
        );

    _kernelXY_To_RAverage_R << <block3, threads3 >> > (count, result);
}

void CMeasure::XYDataToRdistri_C(
    const CLGComplex* __restrict__ source,
    UINT* count,
    CLGComplex* result,
    UINT uiMaxR,
    UBOOL bCalculateCounter,
    BYTE byFieldId)
{
    const dim3 block2(_HC_DecompX, 1, 1);
    const dim3 threads2(_HC_DecompLx, 1, 1);
    const dim3 block3(1, 1, 1);
    const dim3 threads3(uiMaxR + 1, 1, 1);

    _kernelInitialDist << <block3, threads3 >> > (
        bCalculateCounter ? count : NULL,
        NULL,
        result
        );

    _kernelXY_To_R_C << <block2, threads2 >> > (
        source,
        CCommonData::m_sCenter,
        uiMaxR,
        byFieldId,
        result,
        bCalculateCounter ? count : NULL
        );

    _kernelXY_To_RAverage_C << <block3, threads3 >> > (count, result);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================