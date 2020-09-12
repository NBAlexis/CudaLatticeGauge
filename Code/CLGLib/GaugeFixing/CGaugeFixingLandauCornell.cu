//=============================================================================
// FILENAME : CGaugeFixingLandauCornell.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/18/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Cornell Steepest Descend

/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateA(
    const deviceSU3* __restrict__ pU,
    Real* pA11,
    CLGComplex* pA12, 
    CLGComplex* pA13,
    Real* pA22,
    CLGComplex* pA23)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex] = su3A.m_me[0].y;
            pA12[uiLinkIndex] = su3A.m_me[1];
            pA13[uiLinkIndex] = su3A.m_me[2];
            pA22[uiLinkIndex] = su3A.m_me[4].y;
            pA23[uiLinkIndex] = su3A.m_me[5];
        }
        else
        {
            pA11[uiLinkIndex] = F(0.0);
            pA12[uiLinkIndex] = _zeroc;
            pA13[uiLinkIndex] = _zeroc;
            pA22[uiLinkIndex] = F(0.0);
            pA23[uiLinkIndex] = _zeroc;
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelCalculateALog(
    const deviceSU3* __restrict__ pU,
    Real* pA11,
    CLGComplex* pA12,
    CLGComplex* pA13,
    Real* pA22,
    CLGComplex* pA23)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A = su3A.Log();
            pA11[uiLinkIndex] = su3A.m_me[0].y;
            pA12[uiLinkIndex] = su3A.m_me[1];
            pA13[uiLinkIndex] = su3A.m_me[2];
            pA22[uiLinkIndex] = su3A.m_me[4].y;
            pA23[uiLinkIndex] = su3A.m_me[5];
        }
        else
        {
            pA11[uiLinkIndex] = F(0.0);
            pA12[uiLinkIndex] = _zeroc;
            pA13[uiLinkIndex] = _zeroc;
            pA22[uiLinkIndex] = F(0.0);
            pA23[uiLinkIndex] = _zeroc;
        }
    }
}

/**
 * Gamma(n) = Delta _{-mu} A(n) = \sum _mu (A_mu(n - mu) - A_mu(n))
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAGradient(
    BYTE byFieldId,
    Real* pGamma11,
    CLGComplex* pGamma12,
    CLGComplex* pGamma13,
    Real* pGamma22,
    CLGComplex* pGamma23,
    const Real* __restrict__ pA11,
    const CLGComplex* __restrict__ pA12,
    const CLGComplex* __restrict__ pA13,
    const Real* __restrict__ pA22,
    const CLGComplex* __restrict__ pA23)
{
    intokernalInt4;

    pGamma11[uiSiteIndex] = F(0.0);
    pGamma12[uiSiteIndex] = _zeroc;
    pGamma13[uiSiteIndex] = _zeroc;
    pGamma22[uiSiteIndex] = F(0.0);
    pGamma23[uiSiteIndex] = _zeroc;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        return;
    }
    
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            pGamma11[uiSiteIndex] = pGamma11[uiSiteIndex] - pA11[uiLinkIndex];
            pGamma12[uiSiteIndex] = _cuCsubf(pGamma12[uiSiteIndex], pA12[uiLinkIndex]);
            pGamma13[uiSiteIndex] = _cuCsubf(pGamma13[uiSiteIndex], pA13[uiLinkIndex]);
            pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] - pA22[uiLinkIndex];
            pGamma23[uiSiteIndex] = _cuCsubf(pGamma23[uiSiteIndex], pA23[uiLinkIndex]);
        }

        //const SIndex site_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_mu];
        //const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_mu.m_uiSiteIndex, dir);
        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& p_m_mu_dir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        //if (!__idx->_deviceIsBondOnSurface(p_m_mu, dir))
        if (!p_m_mu_dir.IsDirichlet())
        {
            const UINT uiLinkIndex2 = _deviceGetLinkIndex(p_m_mu_dir.m_uiSiteIndex, dir);
            if (p_m_mu_dir.NeedToDagger())
            {
                pGamma11[uiSiteIndex] = pGamma11[uiSiteIndex] - pA11[uiLinkIndex2];
                pGamma12[uiSiteIndex] = _cuCsubf(pGamma12[uiSiteIndex], pA12[uiLinkIndex2]);
                pGamma13[uiSiteIndex] = _cuCsubf(pGamma13[uiSiteIndex], pA13[uiLinkIndex2]);
                pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] - pA22[uiLinkIndex2];
                pGamma23[uiSiteIndex] = _cuCsubf(pGamma23[uiSiteIndex], pA23[uiLinkIndex2]);
            }
            else
            {
                pGamma11[uiSiteIndex] = pGamma11[uiSiteIndex] + pA11[uiLinkIndex2];
                pGamma12[uiSiteIndex] = _cuCaddf(pGamma12[uiSiteIndex], pA12[uiLinkIndex2]);
                pGamma13[uiSiteIndex] = _cuCaddf(pGamma13[uiSiteIndex], pA13[uiLinkIndex2]);
                pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] + pA22[uiLinkIndex2];
                pGamma23[uiSiteIndex] = _cuCaddf(pGamma23[uiSiteIndex], pA23[uiLinkIndex2]);
            }
        }
    }
}

/**
 * g(x)=exp(-i a Delta A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateG(
    deviceSU3* pG,
    const Real* __restrict__ pGamma11,
    const CLGComplex* __restrict__ pGamma12,
    const CLGComplex* __restrict__ pGamma13,
    const Real* __restrict__ pGamma22,
    const CLGComplex* __restrict__ pGamma23,
    Real fAlpha)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pG[uiSiteIndex] = deviceSU3::makeSU3Id();
    }
    else
    {
        deviceSU3 pA = deviceSU3::makeSU3TA(
            pGamma12[uiSiteIndex], pGamma13[uiSiteIndex], pGamma23[uiSiteIndex],
            pGamma11[uiSiteIndex], pGamma22[uiSiteIndex]);
        pG[uiSiteIndex] = (0 == _DC_ExpPrecision)
            ? pA.QuickExp(fAlpha)
            : pA.ExpReal(fAlpha, _DC_ExpPrecision);
    }
}

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const deviceSU3 left(pGx[uiSiteIndex]);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU3 res(pGauge[uiLinkDir]);

            const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
            const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_p_mu_site)];

            if (!site_p_mu.IsDirichlet())
            {
                res.MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateTrAGradientSq(
    Real* pDeviceRes,
    const Real* __restrict__ pDeltaA11,
    const CLGComplex* __restrict__ pDeltaA12,
    const CLGComplex* __restrict__ pDeltaA13,
    const Real* __restrict__ pDeltaA22,
    const CLGComplex* __restrict__ pDeltaA23)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex] = F(0.0);
        return;
    }

    Real fAbs1 = _cuCabsf(pDeltaA12[uiSiteIndex]);
    Real fAbs2 = _cuCabsf(pDeltaA13[uiSiteIndex]);
    Real fAbs3 = _cuCabsf(pDeltaA23[uiSiteIndex]);
    Real fM1122 = pDeltaA11[uiSiteIndex] + pDeltaA22[uiSiteIndex];
    pDeviceRes[uiSiteIndex] = F(2.0) * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}


#pragma endregion

#pragma region FFT accelaration

__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable(Real* pP, UINT uiV)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    Real fDenorm = static_cast<Real>(uiDir);
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        fDenorm -= _cos(F(2.0) * PI * sSite4.m_byData4[dir] / static_cast<Real>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = F(0.5);
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 4 / 2(Nd - sum cos)
    pP[uiSiteIndex] = F(8.0) / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC(const Real* __restrict__ realBuffer, CLGComplex* complexBuffer)
{
    intokernal;
    complexBuffer[uiSiteIndex] = _make_cuComplex(realBuffer[uiSiteIndex], F(0.0));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR(const CLGComplex* __restrict__ complexBuffer, Real* realBuffer)
{
    intokernal;
    realBuffer[uiSiteIndex] = complexBuffer[uiSiteIndex].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale(const Real* __restrict__ pP, CLGComplex* fftRes)
{
    intokernal;
    fftRes[uiSiteIndex] = cuCmulf_cr(fftRes[uiSiteIndex], pP[uiSiteIndex]);
}

#pragma endregion

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingLandauCornell)

void CGaugeFixingLandauCornell::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    //========== Initial Settings ==============
    if (!params.FetchValueReal(_T("Alpha"), m_fAlpha))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: Alpha not set, set to 0.08 by defualt."));
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: Accuracy not set, set to 0.00000000001 by defualt."));
    }

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: MaxIterate not set, set to 100000 by defualt."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = 1000;
    if (!params.FetchValueINT(_T("ShowErrorStep"), iValue))
    {
        appParanoiac(_T("CGaugeFixingLandauCornell: ShowErrorStep not set, set to 1000 by defualt."));
    }
    m_iShowErrorStep = iValue;

    iValue = 1;
    if (!params.FetchValueINT(_T("FFT"), iValue))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: FFT not set, set to 1 by defualt."));
    }
    m_bFA = (0 != iValue);

    //========== Initial Buffers ==============
    checkCudaErrors(cudaMalloc((void**)& m_pA11, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA12, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA13, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA22, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA23, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));

    checkCudaErrors(cudaMalloc((void**)& m_pGamma11, _HC_Volume * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma12, _HC_Volume * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma13, _HC_Volume * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma22, _HC_Volume * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma23, _HC_Volume * sizeof(CLGComplex)));

    checkCudaErrors(cudaMalloc((void**)& m_pG, _HC_Volume * sizeof(deviceSU3)));
    if (m_bFA)
    {
        checkCudaErrors(cudaMalloc((void**)& m_pMomentumTable, _HC_Volume * sizeof(Real)));
        checkCudaErrors(cudaMalloc((void**)& m_pTempFFTBuffer, _HC_Volume * sizeof(CLGComplex)));

        preparethread;
        _kernelBakeMomentumTable << <block, threads >> > (m_pMomentumTable, _HC_Volume);
    }
}

void CGaugeFixingLandauCornell::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge || EFT_GaugeSU3 != pResGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
    deviceSU3* pDeviceBufferPointer = pGaugeSU3->m_pDeviceData;

    preparethread;
    m_iIterate = 0;
    Real fTheta = F(0.0);

    while (m_iIterate < m_iMaxIterate)
    {
        //======= 1. Calculate Gamma    =========
        if (0 == _HC_ALog)
        {
            _kernelCalculateA << <block, threads >> > (
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }
        else
        {
            _kernelCalculateALog << <block, threads >> > (
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }

        _kernelCalculateAGradient << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);

        //======= 2. Calculate Theta    =========
        _kernelCalculateTrAGradientSq << <block, threads >> > (
            _D_RealThreadBuffer,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23);
        fTheta = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer) / (3 * _HC_Volume);
        if (m_iShowErrorStep > 0 && 0 == m_iIterate % m_iShowErrorStep)
        {
            appParanoiac(_T("Theta%d = %2.12f\n"), m_iIterate, fTheta);
        }
        
        if (fTheta < m_fAccuracy)
        {
            return;
        }

        //======= 3. FFT                =========
        if (m_bFA)
        {
            CCLGFFTHelper::FFT4D(m_pGamma12, TRUE);
            _kernelFFTScale << <block, threads >> > (m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::FFT4D(m_pGamma12, FALSE);

            CCLGFFTHelper::FFT4D(m_pGamma13, TRUE);
            _kernelFFTScale << <block, threads >> > (m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::FFT4D(m_pGamma13, FALSE);

            CCLGFFTHelper::FFT4D(m_pGamma23, TRUE);
            _kernelFFTScale << <block, threads >> > (m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::FFT4D(m_pGamma23, FALSE);

            _kernelFFTRtoC << <block, threads >> > (m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, TRUE);
            _kernelFFTScale << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, FALSE);
            _kernelFFTCtoR << <block, threads >> > (m_pTempFFTBuffer, m_pGamma11);

            _kernelFFTRtoC << <block, threads >> > (m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, TRUE);
            _kernelFFTScale << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, FALSE);
            _kernelFFTCtoR << <block, threads >> > (m_pTempFFTBuffer, m_pGamma22);
        }

        //======= 4. Gauge Transform    =========
        //Be careful not to use strictLog when A is really small
        _kernelCalculateG << <block, threads >> > (
            m_pG,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_fAlpha);

        _kernelGaugeTransform << <block, threads >> > (pResGauge->m_byFieldId, m_pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("Gauge fixing failed with last error = %f\n"), fTheta);
}

Real CGaugeFixingLandauCornell::CheckRes(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return F(0.0);
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelCalculateA << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);
    }
    else
    {
        _kernelCalculateALog << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);
    }

    _kernelCalculateAGradient << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        m_pGamma11,
        m_pGamma12,
        m_pGamma13,
        m_pGamma22,
        m_pGamma23,
        m_pA11,
        m_pA12,
        m_pA13,
        m_pA22,
        m_pA23);

    _kernelCalculateTrAGradientSq << <block, threads >> > (
        _D_RealThreadBuffer,
        m_pGamma11,
        m_pGamma12,
        m_pGamma13,
        m_pGamma22,
        m_pGamma23);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer) / (3 * _HC_Volume);
}

CCString CGaugeFixingLandauCornell::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingLandauCornell\n");
    sRet = sRet + tab + _T("accuray : ") + appFloatToString(m_fAccuracy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================