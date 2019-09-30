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

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {

        //======================================================
        //4-chair terms except for the last one
        betaOverN = F(0.125) * betaOverN;

        const Real fX = (sSite4.x - sCenter.x);

        //===============
        //+x Omega V412
        const Real fV412 = fX * _deviceChairTerm(pDeviceData, 3, 0, 1, uiN);

        //===============
        //+x Omega V432
        const Real fV432 = fX * _deviceChairTerm(pDeviceData, 3, 2, 1, uiN);

        const Real fY = -(sSite4.y - sCenter.y);

        //===============
        //-y Omega V421
        const Real fV421 = fY * _deviceChairTerm(pDeviceData, 3, 1, 0, uiN);

        //===============
        //-y Omega V431
        const Real fV431 = fY * _deviceChairTerm(pDeviceData, 3, 2, 0, uiN);

        fRes = (fV412 + fV432 + fV421 + fV431) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGaugeSpin(
    const deviceSU3* __restrict__ pE, 
    const deviceSU3* __restrict__ piAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    
    if (!site.IsDirichlet())
    {
        //Note i A phys is i x A phys
        UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
        UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
        deviceSU3 beforeTrace = pE[uiLinkX].MulC(piAphys[uiLinkY]);
        beforeTrace.Sub(pE[uiLinkY].MulC(piAphys[uiLinkX]));
        CLGComplex cRes = cuCmulf_cr(beforeTrace.Tr(), -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], cRes.x);
    }
}

/**
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChen(
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ piUpure,
    const deviceSU3* __restrict__ piAphys,
    SSmallInt4 sCenter,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;
    
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    CLGComplex res = _zeroc;
    const Real fX = (sSite4.x - sCenter.x);
    const Real fY = (sSite4.y - sCenter.y);
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 DyAphys = _deviceDPureMuUpure(piAphys, piUpure, uiBigIdx, 1, dir);
            deviceSU3 DxAphys = _deviceDPureMuUpure(piAphys, piUpure, uiBigIdx, 0, dir);
            DyAphys.MulReal(fX);
            DxAphys.MulReal(fY);
            DyAphys.Sub(DxAphys);

            DyAphys = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(DyAphys);
            _cuCaddf(res, DyAphys.Tr());
        }

        res = cuCmulf_cr(res, -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

#pragma endregion

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

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionJG)
    {
        checkCudaErrors(cudaFree(m_pDistributionJG));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionJG)
    {
        free(m_pHostDistributionJG);
    }

    appSafeFree(m_pHostSpinBuffer);
    cudaSafeFree(m_pDeviceSpinBuffer);
    appSafeFree(m_pHostDistributionJGS);
    cudaSafeFree(m_pDistributionJGS);
    appSafeDelete(m_pE);
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

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureSpin"), iValue);
    m_bMeasureSpin = iValue != 0;

    if (m_bMeasureSpin)
    {
        m_pHostSpinBuffer = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
        checkCudaErrors(cudaMalloc((void**)& m_pDeviceSpinBuffer, sizeof(Real) * _HC_Lx * _HC_Ly));
        m_pE = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(m_byFieldId)->GetCopy());
    }

    if (m_bMeasureDistribution)
    {
        //assuming the center is really at center
        m_uiMaxR = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
                 + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2);

        m_uiEdgeR = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1);

        checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionJG = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));

        if (m_bMeasureSpin)
        {
            m_pHostDistributionJGS = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
            checkCudaErrors(cudaMalloc((void**)& m_pDistributionJGS, sizeof(Real) * (m_uiMaxR + 1)));
        }
    }
}

void CMeasureAMomentumJG::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    const Real fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(_HC_SUN);

    _ZeroXYPlane(m_pDeviceDataBufferOneConfig);

    preparethread;

    _kernelCalculateAngularMomentumJG << <block, threads >> > (
        pGaugeSU3->m_pDeviceData, 
        m_pDeviceDataBufferOneConfig,
        CCommonData::m_sCenter,
        fBetaOverN,
        m_byFieldId);

    _AverageXYPlane(m_pDeviceDataBufferOneConfig);

    checkCudaErrors(cudaGetLastError());

    if (m_bMeasureDistribution)
    {
        XYDataToRdistri_R(m_pDeviceDataBufferOneConfig, m_pDistributionR, m_pDistributionJG,
            m_uiMaxR, TRUE, m_byFieldId);

        checkCudaErrors(cudaGetLastError());

        //extract res
        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        FillDataWithR_R(
            m_lstJG, m_lstJGInner, m_lstJGAll, m_lstR, 
            m_pHostDistributionJG, m_pHostDistributionR, 
            m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, TRUE
            );
    }

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
            appDetailed(_T("(%d,%d)=%1.6f  "), 
                i, 
                CCommonData::m_sCenter.y, 
                m_pHostDataBuffer[i * _HC_Ly + CCommonData::m_sCenter.y]);
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n"));
    }


    if (m_bMeasureSpin)
    {
        _ZeroXYPlane(m_pDeviceSpinBuffer);
        pGaugeSU3->CalculateE_Using_U(m_pE);
        CFieldGaugeSU3* pESU3 = dynamic_cast<CFieldGaugeSU3*>(m_pE);
        CFieldGaugeSU3* pAphysSU3 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pAphys);
        if (NULL == pAphysSU3)
        {
            appCrucial(_T("CMeasureAMomentumJG: A phys not calculated\n"));
        }
        else
        {
            _kernelCalculateGaugeSpin << <block, threads >> > (pESU3->m_pDeviceData, pAphysSU3->m_pDeviceData, m_pDeviceSpinBuffer, fBetaOverN);

            _AverageXYPlane(m_pDeviceSpinBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostSpinBuffer, m_pDeviceSpinBuffer, sizeof(Real)* _HC_Lx* _HC_Ly, cudaMemcpyDeviceToHost));
            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGS.AddItem(m_pHostSpinBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_pDeviceSpinBuffer, m_pDistributionR, m_pDistributionJGS, 
                    m_uiMaxR, FALSE, m_byFieldId);
                
                checkCudaErrors(cudaMemcpy(m_pHostDistributionJGS, m_pDistributionJGS, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGS, m_lstJGSInner, m_lstJGSAll, m_lstR,
                    m_pHostDistributionJGS, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                );
            }

            CFieldGaugeSU3* pUpure = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pUpure);
            if (NULL == pUpure)
            {
                appCrucial(_T("CMeasureAMomentumJG: A pure not calculated\n"));
            }
            else
            {
                _ZeroXYPlane(m_pDeviceSpinBuffer);
                _kernelMomemtumJGChen << <block, threads >> > (
                    pESU3->m_pDeviceData,
                    pUpure->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    CCommonData::m_sCenter,
                    m_pDeviceSpinBuffer,
                    fBetaOverN);

                _AverageXYPlane(m_pDeviceSpinBuffer);
                checkCudaErrors(cudaMemcpy(m_pHostSpinBuffer, m_pDeviceSpinBuffer, sizeof(Real)* _HC_Lx* _HC_Ly, cudaMemcpyDeviceToHost));
                for (UINT i = 1; i < _HC_Ly; ++i)
                {
                    for (UINT j = 1; j < _HC_Lx; ++j)
                    {
                        m_lstResJGChen.AddItem(m_pHostSpinBuffer[j * _HC_Ly + i]);
                    }
                }

                if (m_bMeasureDistribution)
                {
                    XYDataToRdistri_R(
                        m_pDeviceSpinBuffer, m_pDistributionR, m_pDistributionJGS, 
                        m_uiMaxR, FALSE, m_byFieldId);

                    //extract res
                    checkCudaErrors(cudaMemcpy(m_pHostDistributionJGS, m_pDistributionJGS, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                    FillDataWithR_R(
                        m_lstJGChen, m_lstJGChenInner, m_lstJGChenAll, m_lstR,
                        m_pHostDistributionJGS, m_pHostDistributionR,
                        m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                    );
                }
            }
        }
    }

    ++m_uiConfigurationCount;
}

void CMeasureAMomentumJG::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumJG::Report()
{
    appSetLogDate(FALSE);

    appGeneral(_T("\n===================================================\n"));
    appGeneral(_T("=========== Angular Momentum JG of sites ==========\n"), CCommonData::m_sCenter.x);
    appGeneral(_T("===================================================\n"));

    ReportDistributionXY_R(m_uiConfigurationCount, m_lstRes);

    appGeneral(_T("===================================================\n"));

    if (m_bMeasureSpin)
    {
        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JGS of sites ==========\n"), CCommonData::m_sCenter.x);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGS);

        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JG Chen of sites ==========\n"), CCommonData::m_sCenter.x);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChen);
    }

    appGeneral(_T("===================================================\n"));
    appGeneral(_T("===================================================\n"));

    appSetLogDate(TRUE);
}

void CMeasureAMomentumJG::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstRes.RemoveAll();
    m_lstResJGS.RemoveAll();
    m_lstResJGChen.RemoveAll();

    m_lstR.RemoveAll();
    m_lstJG.RemoveAll();
    m_lstJGS.RemoveAll();
    m_lstJGChen.RemoveAll();

    m_lstJGAll.RemoveAll();
    m_lstJGInner.RemoveAll();
    m_lstJGSAll.RemoveAll();
    m_lstJGSInner.RemoveAll();
    m_lstJGChenAll.RemoveAll();
    m_lstJGChenInner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================