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
 * Use Apure directly
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChenApprox(
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pApure,
    const deviceSU3* __restrict__ pAphys,
    SSmallInt4 sCenter,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const Real fmY = -static_cast<Real>(sSite4.y - sCenter.y);
    const Real fmX = -static_cast<Real>(sSite4.x - sCenter.x);

    CLGComplex res = _zeroc;
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 DxAphys = _deviceDPureMu(pAphys, pApure, uiBigIdx, 0, dir);
            DxAphys.MulReal(fmY);
            deviceSU3 DyAphys = _deviceDPureMu(pAphys, pApure, uiBigIdx, 1, dir);
            DyAphys.MulReal(fmX);
            DyAphys.Sub(DxAphys);
            res = _cuCaddf(res, _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(DyAphys).Tr());
        }

        res = cuCmulf_cr(res, -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChenApprox2(
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pApure,
    const deviceSU3* __restrict__ pAphys,
    SSmallInt4 sCenter,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const Real fmY = -static_cast<Real>(sSite4.y - sCenter.y);
    const Real fmX = -static_cast<Real>(sSite4.x - sCenter.x);

    CLGComplex res = _zeroc;
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 DxAphys = _deviceDPureMu2(pAphys, pApure, uiBigIdx, 0, dir);
            DxAphys.MulReal(fmY);
            deviceSU3 DyAphys = _deviceDPureMu2(pAphys, pApure, uiBigIdx, 1, dir);
            DyAphys.MulReal(fmX);
            DyAphys.Sub(DxAphys);
            res = _cuCaddf(res, _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(DyAphys).Tr());
        }

        res = cuCmulf_cr(res, -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

/**
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChen(
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pXcrossDpureA,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalInt4;
    
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    CLGComplex res = _zeroc;
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 beforeTrace = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(pXcrossDpureA[_deviceGetLinkIndex(uiSiteIndex, dir)]);
            res = _cuCaddf(res, beforeTrace.Tr());
        }

        res = cuCmulf_cr(res, fBetaOverN * F(0.5));
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomentumJGChenDpureA(
    deviceSU3* pXcrossDpureA,
    const deviceSU3* __restrict__ pGauge,
    const deviceSU3* __restrict__ pAphys,
    SSmallInt4 sCenter)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const Real fmY = -static_cast<Real>(sSite4.y - sCenter.y);
    const Real fmX = -static_cast<Real>(sSite4.x - sCenter.x);

    const UINT x_m_x_Gauge = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 0];
    const UINT x_m_y_Gauge = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 1];
    const UINT x_p_x_Gauge = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 0 + uiDir];
    const UINT x_p_y_Gauge = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 1 + uiDir];

    //idir = mu
    for (UINT idir = 0; idir < uiDir - 1; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pXcrossDpureA[uiLinkIndex] = deviceSU3::makeSU3Zero();

        //U_x(n) A_dir(n+x)U_x^+(n)
        deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 0));
        deviceSU3 a(_deviceGetGaugeBCSU3DirZero(pAphys, x_p_x_Gauge, idir));
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmY);
        pXcrossDpureA[uiLinkIndex].Add(u);

        //U_x^+(n-x) A_dir(n-x) U_x(n-x)
        u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_x_Gauge, 0);
        a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_x_Gauge, idir);
        a.Mul(u);
        u.DaggerMul(a);
        u.MulReal(fmY);
        pXcrossDpureA[uiLinkIndex].Sub(u);

        //U_y(n) A_dir(n+y)U_y^+(n)
        u = _deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 1);
        a = _deviceGetGaugeBCSU3DirZero(pAphys, x_p_y_Gauge, idir);
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmX);
        pXcrossDpureA[uiLinkIndex].Sub(u);

        //U_y^+(n-y) A_dir(n-y) U_y(n-y)
        u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_y_Gauge, 1);
        a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_y_Gauge, idir);
        a.Mul(u);
        u.DaggerMul(a);
        u.MulReal(fmX);
        pXcrossDpureA[uiLinkIndex].Add(u);
    }
}

#pragma endregion

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

    appSafeDelete(m_pE);
    appSafeDelete(m_pDpureA);
}

void CMeasureAMomentumJG::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly));
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

    iValue = 0;
    param.FetchValueINT(_T("MeasureApprox"), iValue);
    m_bMeasureApprox = iValue != 0;

    if (m_bMeasureSpin)
    {
        m_pE = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(m_byFieldId)->GetCopy());
        m_pDpureA = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(m_byFieldId)->GetCopy());
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

    _ZeroXYPlane(m_pDeviceDataBuffer);

    preparethread;

    _kernelCalculateAngularMomentumJG << <block, threads >> > (
        pGaugeSU3->m_pDeviceData, 
        m_pDeviceDataBuffer,
        CCommonData::m_sCenter,
        fBetaOverN,
        m_byFieldId);

    _AverageXYPlane(m_pDeviceDataBuffer);

    checkCudaErrors(cudaGetLastError());

    if (m_bMeasureDistribution)
    {
        XYDataToRdistri_R(m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
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

        checkCudaErrors(cudaGetLastError());
    }

    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

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
        _ZeroXYPlane(m_pDeviceDataBuffer);
        pGaugeSU3->CalculateE_Using_U(m_pE);
        CFieldGaugeSU3* pESU3 = dynamic_cast<CFieldGaugeSU3*>(m_pE);
        const CFieldGaugeSU3* pAphysSU3 = dynamic_cast<const CFieldGaugeSU3*>(appGetLattice()->m_pAphys);
        CFieldGaugeSU3* pDpureA = dynamic_cast<CFieldGaugeSU3*>(m_pDpureA);
        if (NULL == pAphysSU3)
        {
            appCrucial(_T("CMeasureAMomentumJG: A phys not calculated\n"));
        }
        else
        {
            _kernelCalculateGaugeSpin << <block, threads >> > (
                pESU3->m_pDeviceData, 
                pAphysSU3->m_pDeviceData, 
                m_pDeviceDataBuffer,
                fBetaOverN);

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real)* _HC_Lx* _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGS.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, m_byFieldId);
                
                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGS, m_lstJGSInner, m_lstJGSAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                );
            }

            _ZeroXYPlane(m_pDeviceDataBuffer);

            _kernelMomentumJGChenDpureA << <block, threads >> > (
                pDpureA->m_pDeviceData,
                pGaugeSU3->m_pDeviceData,
                pAphysSU3->m_pDeviceData,
                CCommonData::m_sCenter
                );
            _kernelMomemtumJGChen << <block, threads >> > (
                pESU3->m_pDeviceData,
                pDpureA->m_pDeviceData,
                m_pDeviceDataBuffer,
                fBetaOverN);

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGChen.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, m_byFieldId);

                //extract res
                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGChen, m_lstJGChenInner, m_lstJGChenAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                );

                checkCudaErrors(cudaGetLastError());
            }

            if (m_bMeasureApprox)
            {
                pGaugeSU3->CopyTo(pDpureA);
                pDpureA->TransformToIA();
                //[A, Aphys] = [Apure, Aphys], so we do not need to minus Aphys
                //pDpureA->AxpyMinus(pAphysSU3);

                _ZeroXYPlane(m_pDeviceDataBuffer);

                _kernelMomemtumJGChenApprox << <block, threads >> > (
                    pESU3->m_pDeviceData,
                    pDpureA->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    CCommonData::m_sCenter,
                    m_pDeviceDataBuffer,
                    fBetaOverN
                    );

                _AverageXYPlane(m_pDeviceDataBuffer);
                checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real)* _HC_Lx* _HC_Ly, cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaGetLastError());

                for (UINT i = 1; i < _HC_Ly; ++i)
                {
                    for (UINT j = 1; j < _HC_Lx; ++j)
                    {
                        m_lstResJGChenApprox.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                    }
                }

                if (m_bMeasureDistribution)
                {
                    XYDataToRdistri_R(
                        m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                        m_uiMaxR, FALSE, m_byFieldId);

                    //extract res
                    checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                    FillDataWithR_R(
                        m_lstJGChenApprox, m_lstJGChenApproxInner, m_lstJGChenApproxAll, m_lstR,
                        m_pHostDistributionJG, m_pHostDistributionR,
                        m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                    );
                    checkCudaErrors(cudaGetLastError());
                }

                _ZeroXYPlane(m_pDeviceDataBuffer);

                _kernelMomemtumJGChenApprox2 << <block, threads >> > (
                    pESU3->m_pDeviceData,
                    pDpureA->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    CCommonData::m_sCenter,
                    m_pDeviceDataBuffer,
                    fBetaOverN
                    );

                _AverageXYPlane(m_pDeviceDataBuffer);
                checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaGetLastError());

                for (UINT i = 1; i < _HC_Ly; ++i)
                {
                    for (UINT j = 1; j < _HC_Lx; ++j)
                    {
                        m_lstResJGChenApprox2.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                    }
                }

                if (m_bMeasureDistribution)
                {
                    XYDataToRdistri_R(
                        m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                        m_uiMaxR, FALSE, m_byFieldId);

                    //extract res
                    checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                    FillDataWithR_R(
                        m_lstJGChenApprox2, m_lstJGChenApprox2Inner, m_lstJGChenApprox2All, m_lstR,
                        m_pHostDistributionJG, m_pHostDistributionR,
                        m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, FALSE
                    );
                    checkCudaErrors(cudaGetLastError());
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

        if (m_bMeasureApprox)
        {
            appGeneral(_T("\n========================================================\n"));
            appGeneral(_T("=========== Angular Momentum JG Chen Approx of sites ==========\n"), CCommonData::m_sCenter.x);
            appGeneral(_T("========================================================\n"));

            ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChenApprox);

            appGeneral(_T("\n========================================================\n"));
            appGeneral(_T("=========== Angular Momentum JG Chen Approx 2 of sites ==========\n"), CCommonData::m_sCenter.x);
            appGeneral(_T("========================================================\n"));

            ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChenApprox2);
        }
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
    m_lstResJGChenApprox.RemoveAll();
    m_lstResJGChenApprox2.RemoveAll();

    m_lstR.RemoveAll();
    m_lstJG.RemoveAll();
    m_lstJGS.RemoveAll();
    m_lstJGChen.RemoveAll();
    m_lstJGChenApprox.RemoveAll();
    m_lstJGChenApprox2.RemoveAll();

    m_lstJGAll.RemoveAll();
    m_lstJGInner.RemoveAll();
    m_lstJGSAll.RemoveAll();
    m_lstJGSInner.RemoveAll();
    m_lstJGChenAll.RemoveAll();
    m_lstJGChenInner.RemoveAll();
    m_lstJGChenApproxAll.RemoveAll();
    m_lstJGChenApproxInner.RemoveAll();
    m_lstJGChenApprox2All.RemoveAll();
    m_lstJGChenApprox2Inner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================