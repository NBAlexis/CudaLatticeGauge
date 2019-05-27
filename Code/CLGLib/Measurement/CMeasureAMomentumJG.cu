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

CMeasureAMomentumJG::~CMeasureAMomentumJG()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
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
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly));
    Reset();

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
        CCommonData::m_sCenter,
        fBetaOverN,
        m_byFieldId);

    _AverageXYPlane(m_pDeviceDataBufferOneConfig);

    ++m_uiConfigurationCount;
    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

    if (m_bShowResult)
    {
        appDetailed(_T(" === Angular Momentum JG of site y=%d ======\n"), CCommonData::m_sCenter.y);
    }
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstRes.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
        }
    }

    for (UINT i = 1; i < _HC_Lx; ++i)
    {
        if (m_bShowResult)
        {
            appDetailed(_T("%d=%1.6f  "), i, m_pHostDataBuffer[i * _HC_Ly + CCommonData::m_sCenter.y]);
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n"));
    }

}

void CMeasureAMomentumJG::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumJG::Report()
{
    appSetLogDate(FALSE);

    assert(m_uiConfigurationCount * (_HC_Lx - 1) * (_HC_Ly - 1) 
        == static_cast<UINT>(m_lstRes.Num()));
    appGeneral(_T("\n===================================================\n"));
    appGeneral(_T("=========== Angular Momentum JG of sites ==========\n"), CCommonData::m_sCenter.x);
    appGeneral(_T("===================================================\n"));

    TArray<Real> tmp;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                appGeneral(_T("%1.6f, "), m_lstRes[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                if (0 == k)
                {
                    tmp.AddItem(m_lstRes[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                }
                else
                {
                    tmp[i * (_HC_Lx - 1) + j] += m_lstRes[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j];
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
            appGeneral(_T("(x=%d,y=%d)%2.8f,   "), j, i, tmp[i * (_HC_Lx - 1) + j] / m_uiConfigurationCount);
        }
        appGeneral(_T("\n"));
    }

    appGeneral(_T("===================================================\n"));
    appGeneral(_T("===================================================\n"));

    appSetLogDate(TRUE);
}

void CMeasureAMomentumJG::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstRes.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================