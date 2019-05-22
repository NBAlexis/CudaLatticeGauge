//=============================================================================
// FILENAME : CMeasureAMomemtumJG.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAMomentumJG)

#pragma region kernles

/**
* Initial as zero
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero_XYPlane(Real* pBuffer)
{
    intokernalOnlyInt4;

    if (0 == sSite4.z && 0 == sSite4.w)
    {
        pBuffer[sSite4.x * _DC_Ly + sSite4.y] = F(0.0);
    }
}

/**
* Average over z and t
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverZT_XYPlane(Real* pBuffer)
{
    intokernalOnlyInt4;

    if (0 == sSite4.z && 0 == sSite4.w)
    {
        UINT uiIdx = sSite4.x * _DC_Ly + sSite4.y;
        pBuffer[uiIdx] = pBuffer[uiIdx] / (_DC_Lz * _DC_Lt);
    }
}

/**
* calculate momentum, and sum over y and z
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumJG(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer, 
    SSmallInt4 sCenter, 
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {

        //======================================================
        //4-chair terms except for the last one
        betaOverN = F(0.125) * betaOverN;

        Real fX = (sSite4.x - sCenter.x);

        //===============
        //+x Omega V412
        Real fV412 = fX * _deviceChairTerm(pDeviceData, 3, 0, 1, uiN);

        //===============
        //+x Omega V432
        Real fV432 = fX * _deviceChairTerm(pDeviceData, 3, 2, 1, uiN);

        Real fY = -(sSite4.y - sCenter.y);

        //===============
        //-y Omega V421
        Real fV421 = fY * _deviceChairTerm(pDeviceData, 3, 1, 0, uiN);

        //===============
        //-y Omega V431
        Real fV431 = fY * _deviceChairTerm(pDeviceData, 3, 2, 0, uiN);

        fRes = (fV412 + fV432 + fV421 + fV431) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes);
}

//sum over all config
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverConfigurations(
    Real* pBufferAllConfig,
    const Real* __restrict__ pBufferOneConfig,
    UINT uiConfigCount)
{
    intokernalOnlyInt4;

    if (0 == sSite4.z && 0 == sSite4.w)
    {
        Real fConfigCount = static_cast<Real>(uiConfigCount);
        Real fFactor1 = (fConfigCount - F(1.0)) / fConfigCount;
        Real fFactor2 = F(1.0) / fConfigCount;
        UINT uiIdx = sSite4.x * _DC_Ly + sSite4.y;
        pBufferAllConfig[uiIdx]
            = pBufferAllConfig[uiIdx] * fFactor1
            + pBufferOneConfig[uiIdx] * fFactor2;
    }
}

#pragma endregion

/**
* array[x, y] = array[x, y] / (lz * lt)
*/
static void _AverageXYPlane(Real* pDeviceRes)
{
    preparethread;
    _kernelAverageOverZT_XYPlane << <block, threads >> > (pDeviceRes);
}

/**
* array[x, y] = 0
*/
static void _ZeroXYPlane(Real* pDeviceRes)
{
    preparethread;
    _kernelInitialZero_XYPlane << <block, threads >> > (pDeviceRes);
}

/**
* array[x, y] = (array[x, y] * (N-1) + oneconfiguration[x, y]) / N
*/
static void _AverageXYPlaneOverConf(Real* pDeviceRes, const Real* __restrict__ pDeviceResOneConfig, UINT uiConfigCount)
{
    preparethread;
    _kernelAverageOverConfigurations << <block, threads >> > (pDeviceRes, pDeviceResOneConfig, uiConfigCount);
}

CMeasureAMomentumJG::~CMeasureAMomentumJG()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBuffer));
    }
    if (NULL != m_pDeviceDataBufferOneConfig)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBufferOneConfig));
    }
}

void CMeasureAMomentumJG::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly));
    Reset();

    //get center
    TArray<INT> centerArray;
    param.FetchValueArrayINT(_T("Center"), centerArray);
    if (centerArray.Num() > 3)
    {
        m_sCenter.x = static_cast<SBYTE>(centerArray[0]);
        m_sCenter.y = static_cast<SBYTE>(centerArray[1]);
        m_sCenter.z = static_cast<SBYTE>(centerArray[2]);
        m_sCenter.w = static_cast<SBYTE>(centerArray[3]);
    }

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasureAMomentumJG::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    Real fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(_HC_SUN);

    _ZeroXYPlane(m_pDeviceDataBufferOneConfig);

    preparethread;

    _kernelCalculateAngularMomentumJG << <block, threads >> > (
        pGaugeSU3->m_pDeviceData, 
        m_pDeviceDataBufferOneConfig,
        m_sCenter,
        fBetaOverN,
        m_byFieldId);

    _AverageXYPlane(m_pDeviceDataBufferOneConfig);

    ++m_uiConfigurationCount;
    if (1 == m_uiConfigurationCount)
    {
        checkCudaErrors(cudaMemcpy(m_pDeviceDataBuffer, m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToDevice));
    }
    else
    {
        _AverageXYPlaneOverConf(m_pDeviceDataBuffer, m_pDeviceDataBufferOneConfig, m_uiConfigurationCount);
    }

    if (m_bShowResult)
    {
        checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
        appGeneral(_T(" === Angular Momentum JG of site y=%d ======\n"), m_sCenter.y);
        for (UINT i = 1; i < _HC_Lx; ++i)
        {
            appGeneral(_T("%d=%1.6f  "), i, m_pHostDataBuffer[i * _HC_Ly + m_sCenter.y]);
        }
        appGeneral(_T("\n"));
    }

}

void CMeasureAMomentumJG::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumJG::Report()
{
    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    appGeneral(_T("\n=========== Angular Momentum JG of sites ==========\n"), m_sCenter.x);
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            appGeneral(_T("(%d,%d)=%1.6f  "), j, i, m_pHostDataBuffer[j * _HC_Ly + i]);
        }
        appGeneral(_T("\n"));
    }
    appGeneral(_T("===================================================\n"));
}

void CMeasureAMomentumJG::Reset()
{
    m_uiConfigurationCount = 0;
    _ZeroXYPlane(m_pDeviceDataBuffer);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================