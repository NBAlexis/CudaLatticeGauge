//=============================================================================
// FILENAME : CMeasure.cpp
// 
// DESCRIPTION:
// Some common functions
//
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
    UINT uiMax, 
    BYTE byFieldId,
    CLGComplex* result,
    UINT* pCount,
    UBOOL bShiftCenter)
{
    const UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    const INT iX = static_cast<INT>(uiXY / _DC_Ly);
    const INT iY = static_cast<INT>(uiXY % _DC_Ly);
    INT iC;
    const INT iCenterX = _DC_Centerx;
    const INT iCenterY = _DC_Centery;
    if (bShiftCenter)
    {
        iC = (((iCenterX - iX) * 2) - 1) * (((iCenterX - iX) * 2) - 1)
           + (((iCenterY - iY) * 2) - 1) * (((iCenterY - iY) * 2) - 1);
    }
    else
    {
        iC = (iCenterX - iX) * (iCenterX - iX)
           + (iCenterY - iY) * (iCenterY - iY);
    }

    //In the following, sSite is only used for Dirichlet check
    SSmallInt4 sSite4;
    sSite4.z = _DC_Centerz;
    sSite4.w = _DC_Centert;
    sSite4.x = static_cast<SBYTE>(iX);
    sSite4.y = static_cast<SBYTE>(iY);
    if (iC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (NULL != pCount)
        {
            atomicAdd(&pCount[iC], 1);
        }
        atomicAdd(&result[iC].x, jgsXY[uiXY].x);
        atomicAdd(&result[iC].y, jgsXY[uiXY].y);
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
    UINT uiMax,
    BYTE byFieldId,
    Real* result,
    UINT* pCount,
    UBOOL bShiftCenter)
{
    const UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    const INT iX = static_cast<INT>(uiXY / _DC_Ly);
    const INT iY = static_cast<INT>(uiXY % _DC_Ly);
    INT iC;
    const INT iCenterX = _DC_Centerx;
    const INT iCenterY = _DC_Centery;
    if (bShiftCenter)
    {
        iC = (((iCenterX - iX) * 2) - 1) * (((iCenterX - iX) * 2) - 1)
           + (((iCenterY - iY) * 2) - 1) * (((iCenterY - iY) * 2) - 1);
    }
    else
    {
        iC = (iCenterX - iX) * (iCenterX - iX)
           + (iCenterY - iY) * (iCenterY - iY);
    }

    SSmallInt4 sSite4;
    sSite4.z = _DC_Centerz;
    sSite4.w = _DC_Centert;
    sSite4.x = static_cast<SBYTE>(iX);
    sSite4.y = static_cast<SBYTE>(iY);
    if (iC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (NULL != pCount)
        {
            atomicAdd(&pCount[iC], 1);
        }
        atomicAdd(&result[iC], jgsXY[uiXY]);
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

void CMeasure::Average()
{
    appPushLogDate(FALSE);
    if (m_lstRealResults.Num() > 0)
    {
        assert(m_uiConfigurationCount == static_cast<UINT>(m_lstRealResults.Num()));
        m_fAverageRealRes = F(0.0);
        for (INT i = 0; i < m_lstRealResults.Num(); ++i)
        {
            m_fAverageRealRes += m_lstRealResults[i];
        }
        m_fAverageRealRes = m_fAverageRealRes / m_uiConfigurationCount;
        appParanoiac(_T(" === Averaged (%d measures) === %f\n"), m_uiConfigurationCount, m_fAverageRealRes);
    }

    if (m_lstComplexResults.Num() > 0)
    {
        assert(m_uiConfigurationCount == static_cast<UINT>(m_lstComplexResults.Num()));
        m_cAverageCmpRes = _zeroc;
        for (INT i = 0; i < m_lstRealResults.Num(); ++i)
        {
            m_cAverageCmpRes.x += m_lstComplexResults[i].x;
            m_cAverageCmpRes.y += m_lstComplexResults[i].y;
        }
        m_cAverageCmpRes.x = m_cAverageCmpRes.x / m_uiConfigurationCount;
        m_cAverageCmpRes.y = m_cAverageCmpRes.y / m_uiConfigurationCount;
        appParanoiac(_T(" === Averaged (%d measures) === %f + %f \n"), m_uiConfigurationCount, m_cAverageCmpRes.x, m_cAverageCmpRes.y);
    }
    appPopLogDate();
}

void CMeasure::WriteRealListToFile(const CCString& sFileName) const
{
    WriteRealArray(sFileName, m_lstRealResults);
}

void CMeasure::WriteCmpListToFile(const CCString& sFileName) const
{
    WriteComplexArray(sFileName, m_lstComplexResults);
}

void CMeasure::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_pLatticeData = pLatticeData;
    m_byId = byId;

    INT iNeedGaugeSmearing = 0;
    param.FetchValueINT(_T("GaugeSmearing"), iNeedGaugeSmearing);
    m_bNeedSmearing = 0 != iNeedGaugeSmearing;

    INT iValue = 0;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFermionFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    param.FetchValueArrayBYTE(_T("GaugeFields"), m_lstGaugeFieldIds);
    param.FetchValueArrayBYTE(_T("BosonFields"), m_lstBosonFieldIds);

    if (0 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        m_lstGaugeFieldIds.AddItem(1);
    }
}

void CMeasure::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        OnConfigurationAcceptedSingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx]);
    }
}

void CMeasure::SourceSanning(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        SourceSanningSingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx], sources, site);
    }
}

void CMeasure::OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        OnConfigurationAcceptedZ4SingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx], pZ4, pInverseZ4, bStart, bEnd);
    }
}

void CMeasure::FillDataWithR_R(
    TArray<Real>& arrData,
    TArray<Real>* arrInner,
    TArray<Real>& arrFull,
    TArray<UINT>& arrR,
    Real* hostData,
    UINT* hostR,
    UINT uiConfig,
    UINT uiMaxR,
    UINT uiEdgeR,
    Real fDivider,
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
                
                arrData.AddItem(hostData[uiL] * fDivider);

                uiInnerPointsAll += hostR[uiL];
                fAverageJGAll += hostR[uiL] * hostData[uiL] * fDivider;
                if (NULL != arrInner && uiL < uiEdgeR)
                {
                    uiInnerPointsInner += hostR[uiL];
                    fAverageJGInner += hostR[uiL] * hostData[uiL] * fDivider;
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < arrR.Num(); ++i)
        {
            assert(hostR[arrR[i]] > 0);
            arrData.AddItem(hostData[arrR[i]] * fDivider);

            uiInnerPointsAll += hostR[arrR[i]];
            fAverageJGAll += hostR[arrR[i]] * hostData[arrR[i]] * fDivider;
            if (NULL != arrInner && arrR[i] < uiEdgeR)
            {
                uiInnerPointsInner += hostR[arrR[i]];
                fAverageJGInner += hostR[arrR[i]] * hostData[arrR[i]] * fDivider;
            }
        }
    }

    if (uiInnerPointsAll > 0)
    {
        fAverageJGAll = fAverageJGAll / uiInnerPointsAll;
    }
    if (NULL != arrInner && uiInnerPointsInner > 0)
    {
        fAverageJGInner = fAverageJGInner / uiInnerPointsInner;
    }
    arrFull.AddItem(fAverageJGAll);
    if (NULL != arrInner)
    {
        arrInner->AddItem(fAverageJGInner);
    }
}

void CMeasure::FillDataWithR_C(
    TArray<CLGComplex>& arrData,
    TArray<CLGComplex>* arrInner,
    TArray<CLGComplex>& arrFull,
    TArray<UINT>& arrR,
    CLGComplex* hostData,
    UINT* hostR,
    UINT uiConfig,
    UINT uiMaxR,
    UINT uiEdgeR,
    Real fDivider,
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

                arrData.AddItem(cuCmulf_cr(hostData[uiL], fDivider));

                uiInnerPointsAll += hostR[uiL];
                cAverageJGAll.x += hostR[uiL] * hostData[uiL].x * fDivider;
                cAverageJGAll.y += hostR[uiL] * hostData[uiL].y * fDivider;
                if (NULL != arrInner && uiL < uiEdgeR)
                {
                    uiInnerPointsInner += hostR[uiL];
                    cAverageJGInner.x += hostR[uiL] * hostData[uiL].x * fDivider;
                    cAverageJGInner.y += hostR[uiL] * hostData[uiL].y * fDivider;
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < arrR.Num(); ++i)
        {
            assert(hostR[arrR[i]] > 0);
            arrData.AddItem(cuCmulf_cr(hostData[arrR[i]], fDivider));

            uiInnerPointsAll += hostR[arrR[i]];
            cAverageJGAll.x += hostR[arrR[i]] * hostData[arrR[i]].x * fDivider;
            cAverageJGAll.y += hostR[arrR[i]] * hostData[arrR[i]].y * fDivider;
            if (NULL != arrInner && arrR[i] < uiEdgeR)
            {
                uiInnerPointsInner += hostR[arrR[i]];
                cAverageJGInner.x += hostR[arrR[i]] * hostData[arrR[i]].x * fDivider;
                cAverageJGInner.y += hostR[arrR[i]] * hostData[arrR[i]].y * fDivider;
            }
        }
    }

    if (uiInnerPointsAll > 0)
    {
        cAverageJGAll.x = cAverageJGAll.x / uiInnerPointsAll;
        cAverageJGAll.y = cAverageJGAll.y / uiInnerPointsAll;
    }
    if (NULL != arrInner && uiInnerPointsInner > 0)
    {
        cAverageJGInner.x = cAverageJGInner.x / uiInnerPointsInner;
        cAverageJGInner.y = cAverageJGInner.y / uiInnerPointsInner;
    }
    arrFull.AddItem(cAverageJGAll);
    if (NULL != arrInner)
    {
        arrInner->AddItem(cAverageJGInner);
    }
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
    UBOOL bShiftCenter,
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
        uiMaxR,
        byFieldId,
        result,
        bCalculateCounter ? count : NULL,
        bShiftCenter
        );

    _kernelXY_To_RAverage_R << <block3, threads3 >> > (count, result);
}

void CMeasure::XYDataToRdistri_C(
    UBOOL bShiftCenter,
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
        uiMaxR,
        byFieldId,
        result,
        bCalculateCounter ? count : NULL,
        bShiftCenter
        );

    _kernelXY_To_RAverage_C << <block3, threads3 >> > (count, result);
}

/**
* Sometimes, we need to set bMinus = TRUE, because
* <qbar M q> = - tr[MD^{-1}]
*/
void CMeasure::TransformFromXYDataToRData_C(
    UBOOL bMinus,
    UBOOL bShiftCenter,
    UINT uiMaxR,
    UINT uiEdgeR,
    BYTE byFieldId,
    UINT uiFieldCount,
    UINT uiMeasureCount,
    UINT uiConfig,
    const CLGComplex* const* pXYBuffers,
    UINT* pCountBuffer,
    CLGComplex* pValueBuffer,
    UINT* pHostCountBuffer,
    CLGComplex* pHostValueBuffer,
    TArray<UINT>& lstR,
    TArray<CLGComplex>* lstValues,
    TArray<CLGComplex>* lstAll,
    TArray<CLGComplex>* lstInner)
{
    for (UINT i = 0; i < uiMeasureCount; ++i)
    {
        XYDataToRdistri_C(bShiftCenter, pXYBuffers[i], pCountBuffer, pValueBuffer, uiMaxR, 0 == i, byFieldId);
        if (0 == i)
        {
            checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
        }

        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(CLGComplex) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        FillDataWithR_C(
            lstValues[i],
            NULL == lstInner ? NULL : &(lstInner[i]),
            lstAll[i],
            lstR,
            pHostValueBuffer,
            pHostCountBuffer,
            uiConfig,
            uiMaxR,
            uiEdgeR,
            (bMinus ? F(-1.0) : F(1.0)) / static_cast<Real>(uiFieldCount * _HC_Lz * _HC_Lt),
            0 == i
        );
    }
}

void CMeasure::TransformFromXYDataToRData_R(
    UBOOL bShiftCenter,
    UINT uiMaxR,
    UINT uiEdgeR,
    BYTE byFieldId,
    UINT uiFieldCount,
    UINT uiMeasureCount,
    UINT uiConfig,
    const Real* const* pXYBuffers,
    UINT* pCountBuffer,
    Real* pValueBuffer,
    UINT* pHostCountBuffer,
    Real* pHostValueBuffer,
    TArray<UINT>& lstR,
    TArray<Real>* lstValues,
    TArray<Real>* lstAll,
    TArray<Real>* lstInner)
{
    for (UINT i = 0; i < uiMeasureCount; ++i)
    {
        XYDataToRdistri_R(bShiftCenter, pXYBuffers[i], pCountBuffer, pValueBuffer, uiMaxR, 0 == i, byFieldId);
        if (0 == i)
        {
            checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
        }

        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(Real) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        FillDataWithR_R(
            lstValues[i],
            NULL == lstInner ? NULL : &(lstInner[i]),
            lstAll[i],
            lstR,
            pHostValueBuffer,
            pHostCountBuffer,
            uiConfig,
            uiMaxR,
            uiEdgeR,
            F(1.0) / static_cast<Real>(uiFieldCount * _HC_Lz * _HC_Lt),
            0 == i
        );
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================