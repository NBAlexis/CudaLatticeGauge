//=============================================================================
// FILENAME : CGaugeFixingCoulombCornell.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/20/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CGaugeFixingCoulombCornell.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Cornell Steepest Descend

/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateA3D(
    SCHAR uiT,
    const deviceSU3* __restrict__ pU,
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23,
    BYTE byFieldId)
{
    intokernalInt4_S(uiT);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex3D] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex3D] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex3D] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex3D] = 0.0;
            pA12[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex3D] = 0.0;
            pA23[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelCalculateA3DLog(
    SCHAR uiT,
    const deviceSU3* __restrict__ pU,
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23,
    BYTE byFieldId)
{
    intokernalInt4_S(uiT);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
#pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A = su3A.Log();
            pA11[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex3D] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex3D] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex3D] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex3D] = 0.0;
            pA12[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex3D] = 0.0;
            pA23[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

/**
 * Gamma(n) = Delta _{-mu} A(n) = \sum _mu (A_mu(n - mu) - A_mu(n))
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAGradient3D(
    BYTE byFieldId,
    SCHAR uiT,
    DOUBLE* pGamma11,
    cuDoubleComplex* pGamma12,
    cuDoubleComplex* pGamma13,
    DOUBLE* pGamma22,
    cuDoubleComplex* pGamma23,
    const DOUBLE* __restrict__ pA11,
    const cuDoubleComplex* __restrict__ pA12,
    const cuDoubleComplex* __restrict__ pA13,
    const DOUBLE* __restrict__ pA22,
    const cuDoubleComplex* __restrict__ pA23)
{
    intokernalInt4_S_Only3D(uiT);

    pGamma11[uiSiteIndex3D] = 0.0;
    pGamma12[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);
    pGamma13[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);
    pGamma22[uiSiteIndex3D] = 0.0;
    pGamma23[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet())
    {
        return;
    }

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex3D];
            pGamma12[uiSiteIndex3D] = cuCsub(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex3D]);
            pGamma13[uiSiteIndex3D] = cuCsub(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex3D]);
            pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex3D];
            pGamma23[uiSiteIndex3D] = cuCsub(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex3D]);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, __bck(dir));
        const SIndex& p_m_mu_dir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(p_m_mu_site) + dir];
        const UINT uiLinkIndex2_3D = (p_m_mu_dir.m_uiSiteIndex / _DC_Lt) * (uiDir - 1) + p_m_mu_dir.m_byDir;

        //if (!__idx->_deviceIsBondOnSurface(p_m_mu, dir))
        if (!p_m_mu_dir.IsDirichlet())
        {
            if (p_m_mu_dir.NeedToDagger())
            {
                //dagger means A -> -A
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = cuCsub(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = cuCsub(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = cuCsub(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
            else
            {
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] + pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = cuCadd(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = cuCadd(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] + pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = cuCadd(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
        }
    }
}

/**
 * g(x)=exp(-i a Delta A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateG3D(
    BYTE byFieldId,
    SCHAR uiT,
    deviceSU3* pG,
    const DOUBLE* __restrict__ pGamma11,
    const cuDoubleComplex* __restrict__ pGamma12,
    const cuDoubleComplex* __restrict__ pGamma13,
    const DOUBLE* __restrict__ pGamma22,
    const cuDoubleComplex* __restrict__ pGamma23,
    DOUBLE fAlpha)
{
    intokernalInt4_S_Only3D(uiT);

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pG[uiSiteIndex3D] = deviceSU3::makeSU3Id();
    }
    else
    {
        const deviceSU3 pA = deviceSU3::makeSU3TA(
            _cToRealC(pGamma12[uiSiteIndex3D]), _cToRealC(pGamma13[uiSiteIndex3D]), _cToRealC(pGamma23[uiSiteIndex3D]),
            static_cast<Real>(pGamma11[uiSiteIndex3D]), static_cast<Real>(pGamma22[uiSiteIndex3D]));
        pG[uiSiteIndex3D] = _expreal(pA, fAlpha);
    }
}

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 * U_4(n) is left unchanged
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3D(
    BYTE byFieldId,
    SCHAR uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S(uiT);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const deviceSU3 left(pGx[uiSiteIndex3D]);

    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU3 res(pGauge[uiLinkDir]);

            const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, __fwd(dir));
            const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(p_p_mu_site)];

            if (!site_p_mu.IsDirichlet())
            {
                res.MulDagger(pGx[(site_p_mu.m_uiSiteIndex / _DC_Lt)]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

/**
 * g(n) U_t(n)
 * U_t(n -t) g(n)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3DT(
    BYTE byFieldId,
    SCHAR uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S(uiT);

    //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 3))
    {
        const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
        if (!site.IsDirichlet())
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            pGauge[uiLinkIndex] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex]);
        }
    }

    const SSmallInt4 p_m_t_site = _deviceSmallInt4OffsetC(sSite4, -4);
    const UINT p_m_t_bi = __idx->_deviceGetBigIndex(p_m_t_site);
    const SIndex& site_m_t_point = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][p_m_t_bi];
    const SIndex& site_m_t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][p_m_t_bi * _DC_Dir + 3];
    
    if (!site_m_t.IsDirichlet() && !site_m_t_point.IsDirichlet()) 
    {
        const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_t.m_uiSiteIndex, 3);
        if (site_m_t.NeedToDagger())
        {
            //never here
            printf("ever here???\n");
            pGauge[uiLinkIndex2] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex2]);
        }
        else
        {
            pGauge[uiLinkIndex2].MulDagger(pGx[uiSiteIndex3D]);
        }
    }
}


/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateTrAGradientSq3D(
    BYTE byFieldId,
    SCHAR uiT,
    DOUBLE* pDeviceRes,
    const DOUBLE* __restrict__ pDeltaA11,
    const cuDoubleComplex* __restrict__ pDeltaA12,
    const cuDoubleComplex* __restrict__ pDeltaA13,
    const DOUBLE* __restrict__ pDeltaA22,
    const cuDoubleComplex* __restrict__ pDeltaA23)
{
    intokernalInt4_S_Only3D(uiT);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex3D] = 0.0;
    }

    const DOUBLE fAbs1 = cuCabs(pDeltaA12[uiSiteIndex3D]);
    const DOUBLE fAbs2 = cuCabs(pDeltaA13[uiSiteIndex3D]);
    const DOUBLE fAbs3 = cuCabs(pDeltaA23[uiSiteIndex3D]);
    const DOUBLE fM1122 = pDeltaA11[uiSiteIndex3D] + pDeltaA22[uiSiteIndex3D];
    pDeviceRes[uiSiteIndex3D] = 2.0 * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}


#pragma endregion

#pragma region FFT accelaration

#if !_CLG_DOUBLEFLOAT
__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable3D(DOUBLE* pP, UINT uiV)
{
    intokernalInt4_S_Only3D(0);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE fDenorm = static_cast<DOUBLE>(uiDir - 1);
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        fDenorm -= cos(2.0 * PI * sSite4.m_byData4[dir] / static_cast<DOUBLE>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = 0.5;
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 3 / 2(Nd - sum cos)
    pP[uiSiteIndex3D] = 6.0 / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC3D(const DOUBLE* __restrict__ realBuffer, cuDoubleComplex* complexBuffer)
{
    intokernalInt4_S_Only3D(0);
    complexBuffer[uiSiteIndex3D] = make_cuDoubleComplex(realBuffer[uiSiteIndex3D], 0.0);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR3D(const cuDoubleComplex* __restrict__ complexBuffer, DOUBLE* realBuffer)
{
    intokernalInt4_S_Only3D(0);
    realBuffer[uiSiteIndex3D] = complexBuffer[uiSiteIndex3D].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale3D(const DOUBLE* __restrict__ pP, cuDoubleComplex* fftRes)
{
    intokernalInt4_S_Only3D(0);
    fftRes[uiSiteIndex3D] = cuCmulf_cd(fftRes[uiSiteIndex3D], pP[uiSiteIndex3D]);
}
#else
__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable3D(Real* pP, UINT uiV)
{
    intokernalInt4_S_Only3D(0);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    Real fDenorm = static_cast<Real>(uiDir - 1);
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        fDenorm -= _cos(F(2.0) * PI * sSite4.m_byData4[dir] / static_cast<Real>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = F(0.5);
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 3 / 2(Nd - sum cos)
    pP[uiSiteIndex3D] = F(6.0) / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC3D(const Real* __restrict__ realBuffer, CLGComplex* complexBuffer)
{
    intokernalInt4_S_Only3D(0);
    complexBuffer[uiSiteIndex3D] = _make_cuComplex(realBuffer[uiSiteIndex3D], F(0.0));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR3D(const CLGComplex* __restrict__ complexBuffer, Real* realBuffer)
{
    intokernalInt4_S_Only3D(0);
    realBuffer[uiSiteIndex3D] = complexBuffer[uiSiteIndex3D].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale3D(const Real* __restrict__ pP, CLGComplex* fftRes)
{
    intokernalInt4_S_Only3D(0);
    fftRes[uiSiteIndex3D] = cuCmulf_cr(fftRes[uiSiteIndex3D], pP[uiSiteIndex3D]);
}
#endif
#pragma endregion

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingCoulombCornell)

void CGaugeFixingCoulombCornell::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    //const TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
    //TArray<UINT> latticeDim;
    //latticeDim.AddItem(_HC_Lx);
    //latticeDim.AddItem(_HC_Ly);
    //latticeDim.AddItem(_HC_Lz);
    //TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    //m_pHDecomp[0] = decomp[0];
    //m_pHDecomp[1] = decomp[1];
    //m_pHDecomp[2] = decomp[2];
    //m_pHDecomp[3] = decomp[3];
    //m_pHDecomp[4] = decomp[4];
    //m_pHDecomp[5] = decomp[5];
    m_lstDims.RemoveAll();
    m_lstDims.AddItem(_HC_Lx);
    m_lstDims.AddItem(_HC_Ly);
    m_lstDims.AddItem(_HC_Lz);

    //checkCudaErrors(cudaMalloc((void**)& m_pDDecomp, sizeof(UINT) * 6));

    //========== Initial Settings ==============

#if !_CLG_DOUBLEFLOAT
    if (!params.FetchValueDOUBLE(_T("Alpha"), m_fAlpha))
#else
    if (!params.FetchValueReal(_T("Alpha"), m_fAlpha))
#endif
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: Alpha not set, set to 0.08 by defualt."));
    }
    
    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: Accuracy not set, set to 0.00000000001 by defualt."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
    }
    //m_fAccuracy = m_fAccuracy;

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: MaxIterate not set, set to 100000 by defualt."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = 1000;
    if (!params.FetchValueINT(_T("ShowErrorStep"), iValue))
    {
        appParanoiac(_T("CGaugeFixingCoulombCornell: ShowErrorStep not set, set to 1000 by defualt."));
    }
    m_iShowErrorStep = iValue;

    iValue = 1;
    if (!params.FetchValueINT(_T("FFT"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: FFT not set, set to 1 by defualt."));
    }
    m_bFA = (0 != iValue);

#if !_CLG_DOUBLEFLOAT
    if (m_bFA)
    {
        appCrucial(_T("Do not use FFT for single float point\n"));
        m_bFA = FALSE;
    }
#endif

    //========== Initial Buffers ==============
    checkCudaErrors(__cudaMalloc((void**)&m_pA11, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA12, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA13, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA22, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA23, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));

    checkCudaErrors(__cudaMalloc((void**)&m_pGamma11, _HC_Volume_xyz * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma12, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma13, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma22, _HC_Volume_xyz * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma23, _HC_Volume_xyz * sizeof(cuDoubleComplex)));

    checkCudaErrors(__cudaMalloc((void**)& m_pG, _HC_Volume_xyz * sizeof(deviceSU3)));
    if (m_bFA)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pMomentumTable, _HC_Volume_xyz * sizeof(DOUBLE)));
        checkCudaErrors(__cudaMalloc((void**)&m_pTempFFTBuffer, _HC_Volume_xyz * sizeof(cuDoubleComplex)));

        preparethread_S;
        _LAUNCH_KERNEL(_kernelBakeMomentumTable3D, block3d, threads3d, m_pMomentumTable, _HC_Volume_xyz);
    }
}

void CGaugeFixingCoulombCornell::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge || EFT_GaugeSU3 != pResGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
    deviceSU3* pDeviceBufferPointer = pGaugeSU3->m_pDeviceData;

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.3: same gather -> rank0 global-context fix -> scatter flow as
        //CGaugeFixingCoulombLosAlamos (P4-2.2). The cuFFT plans are re-created per
        //call inside CCLGFFTHelper with m_lstDims, so rebuilding m_lstDims to the
        //global dims (ResizeBuffersToGlobal) is enough to fix the whole lattice.
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(deviceSU3) * _HC_Dir);
        const UINT uiLocalBytes = static_cast<UINT>(sizeof(deviceSU3) * _HC_LinkCount);

        BYTE* pHostLocal = (BYTE*)malloc(uiLocalBytes);
        checkCudaErrors(cudaMemcpy(pHostLocal, pDeviceBufferPointer, uiLocalBytes, cudaMemcpyDeviceToHost));

        UINT uiGlobalBytes = 0;
        BYTE* pGlobal = appGetComm()->GatherFieldToRoot(pHostLocal, uiBytesPerSite, uiGlobalBytes);
        free(pHostLocal);

        if (appGetComm()->IsRoot())
        {
            deviceSU3* pDevGlobal = NULL;
            checkCudaErrors(__cudaMalloc((void**)&pDevGlobal, uiGlobalBytes));
            checkCudaErrors(cudaMemcpy(pDevGlobal, pGlobal, uiGlobalBytes, cudaMemcpyHostToDevice));
            free(pGlobal);
            pGlobal = NULL;

            MGEnterGlobalFixerContext();
            ResizeBuffersToGlobal();
            for (SCHAR uiT = 0; uiT < static_cast<SCHAR>(_HC_Lt); ++uiT)
            {
                GaugeFixingOneTimeSlice(pDevGlobal, uiT, pGaugeSU3->m_byFieldId);
            }
            RestoreLocalBuffers();
            MGExitGlobalFixerContext();

            pGlobal = (BYTE*)malloc(uiGlobalBytes);
            checkCudaErrors(cudaMemcpy(pGlobal, pDevGlobal, uiGlobalBytes, cudaMemcpyDeviceToHost));
            checkCudaErrors(__cudaFree(pDevGlobal));
        }

        BYTE* pLocalOut = (BYTE*)malloc(uiLocalBytes);
        appGetComm()->ScatterFieldFromRoot(pGlobal, uiBytesPerSite, pLocalOut);
        checkCudaErrors(cudaMemcpy(pDeviceBufferPointer, pLocalOut, uiLocalBytes, cudaMemcpyHostToDevice));
        free(pLocalOut);
        if (NULL != pGlobal)
        {
            free(pGlobal);
        }

        //Every local link was rewritten by the scatter memcpy, which bypasses
        //guarded launches; bump the owner handle so the next reader's guard
        //refills the halo (no-op for an unbound handle).
        pResGauge->NotifyWritten();
        return;
    }
#endif

    for (SCHAR uiT = 0; uiT < static_cast<SCHAR>(_HC_Lt); ++uiT)
    {
        GaugeFixingOneTimeSlice(pDeviceBufferPointer, uiT, pGaugeSU3->m_byFieldId);
    }
}

void CGaugeFixingCoulombCornell::GaugeFixingOneTimeSlice(deviceSU3* pDeviceBufferPointer, SCHAR uiT, BYTE byFieldId)
{
    preparethread_S;
    m_iIterate = 0;
    DOUBLE fTheta = 0.0;

    while (m_iIterate < m_iMaxIterate)
    {
        //======= 1. Calculate Gamma    =========
        if (0 == _HC_ALog)
        {
            _LAUNCH_KERNEL(_kernelCalculateA3D, block3d, threads3d, 
                uiT,
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23,
                byFieldId);
        }
        else
        {
            _LAUNCH_KERNEL(_kernelCalculateA3DLog, block3d, threads3d, 
                uiT,
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23,
                byFieldId);
        }

        _LAUNCH_KERNEL(_kernelCalculateAGradient3D, block3d, threads3d, 
            byFieldId,
            uiT,
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
        _LAUNCH_KERNEL(_kernelCalculateTrAGradientSq3D, block3d, threads3d, 
            byFieldId,
            uiT,
            _D_RealThreadBuffer,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23);
        fTheta = appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz);
        if (m_iShowErrorStep > 0 && 0 == m_iIterate % m_iShowErrorStep)
        {
            appParanoiac(_T("Theta%d = %2.12f\n"), m_iIterate, fTheta);
        }

        if (fTheta < m_fAccuracy)
        {
            return;
        }

#if _CLG_DOUBLEFLOAT
#define __FFT3DWithXYZ FFT3DWithXYZ
#else
#define __FFT3DWithXYZ FFT3DWithXYZDouble
#endif
        //======= 3. FFT =========
        if (m_bFA)
        {
            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma12, m_lstDims, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale3D, block3d, threads3d, m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma12, m_lstDims, FALSE);

            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma13, m_lstDims, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale3D, block3d, threads3d, m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma13, m_lstDims, FALSE);

            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma23, m_lstDims, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale3D, block3d, threads3d, m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pGamma23, m_lstDims, FALSE);

            _LAUNCH_KERNEL(_kernelFFTRtoC3D, block3d, threads3d, m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale3D, block3d, threads3d, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR3D, block3d, threads3d, m_pTempFFTBuffer, m_pGamma11);

            _LAUNCH_KERNEL(_kernelFFTRtoC3D, block3d, threads3d, m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale3D, block3d, threads3d, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::__FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR3D, block3d, threads3d, m_pTempFFTBuffer, m_pGamma22);
        }
#undef __FFT3DWithXYZ

        //======= 4. Gauge Transform    =========
        _LAUNCH_KERNEL(_kernelCalculateG3D, block3d, threads3d, 
            byFieldId,
            uiT,
            m_pG,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_fAlpha);
        _LAUNCH_KERNEL(_kernelGaugeTransform3D, block3d, threads3d, byFieldId, uiT, m_pG, pDeviceBufferPointer);
        _LAUNCH_KERNEL(_kernelGaugeTransform3DT, block3d, threads3d, byFieldId, uiT, m_pG, pDeviceBufferPointer);
        ++m_iIterate;
    }

    appGeneral(_T("Gauge fixing failed with last error = %2.15f\n"), fTheta);
}

DOUBLE CGaugeFixingCoulombCornell::CheckRes(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return F(0.0);
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.3: the deviation must be measured on the WHOLE lattice; same flow as
        //CGaugeFixingCoulombLosAlamos::CheckRes (P4-2.2). Non-root ranks report
        //their local deviation; the driver only compares rank 0.
        const UINT uiBytesPerSite = static_cast<UINT>(sizeof(deviceSU3) * _HC_Dir);
        const UINT uiLocalBytes = static_cast<UINT>(sizeof(deviceSU3) * _HC_LinkCount);

        BYTE* pHostLocal = (BYTE*)malloc(uiLocalBytes);
        checkCudaErrors(cudaMemcpy(pHostLocal, pGaugeSU3->m_pDeviceData, uiLocalBytes, cudaMemcpyDeviceToHost));

        UINT uiGlobalBytes = 0;
        BYTE* pGlobal = appGetComm()->GatherFieldToRoot(pHostLocal, uiBytesPerSite, uiGlobalBytes);
        free(pHostLocal);
        if (NULL == pGlobal)
        {
            return CheckResLocal(pGaugeSU3->m_pDeviceData, pGauge->m_byFieldId);
        }

        deviceSU3* pDevGlobal = NULL;
        checkCudaErrors(__cudaMalloc((void**)&pDevGlobal, uiGlobalBytes));
        checkCudaErrors(cudaMemcpy(pDevGlobal, pGlobal, uiGlobalBytes, cudaMemcpyHostToDevice));
        free(pGlobal);

        MGEnterGlobalFixerContext();
        ResizeBuffersToGlobal();
        const DOUBLE fRet = CheckResLocal(pDevGlobal, pGauge->m_byFieldId);
        RestoreLocalBuffers();
        MGExitGlobalFixerContext();

        checkCudaErrors(__cudaFree(pDevGlobal));
        return fRet;
    }
#endif
    return CheckResLocal(pGaugeSU3->m_pDeviceData, pGauge->m_byFieldId);
}

DOUBLE CGaugeFixingCoulombCornell::CheckResLocal(const deviceSU3* pGaugeData, BYTE byFieldId)
{
    DOUBLE fRet = 0.0;

    preparethread_S;
    for (SCHAR uiT = 0; uiT < static_cast<SCHAR>(_HC_Lt); ++uiT)
    {
        if (0 == _HC_ALog)
        {
            _LAUNCH_KERNEL(_kernelCalculateA3D, block3d, threads3d,
                uiT,
                pGaugeData,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23,
                byFieldId);
        }
        else
        {
            _LAUNCH_KERNEL(_kernelCalculateA3DLog, block3d, threads3d,
                uiT,
                pGaugeData,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23,
                byFieldId);
        }

        _LAUNCH_KERNEL(_kernelCalculateAGradient3D, block3d, threads3d,
            byFieldId,
            uiT,
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

        _LAUNCH_KERNEL(_kernelCalculateTrAGradientSq3D, block3d, threads3d,
            byFieldId,
            uiT,
            _D_RealThreadBuffer,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23);

        fRet += appAbs(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz));
    }
    return fRet / _HC_Lt;
}

#if _CLG_MULTI_GPU
void CGaugeFixingCoulombCornell::ResizeBuffersToGlobal()
{
    //Called under the temporary GLOBAL lattice context (rank 0 only), where
    //_HC_Volume_xyz / _HC_Lx/y/z are the GLOBAL values: save the local pointers
    //and dims, re-allocate the fixing buffers to the global 3D volume, rebuild
    //m_lstDims (feeds the cuFFT plans, re-created per call) and re-bake the
    //momentum table.
    m_pSavedA11 = m_pA11;
    m_pSavedA12 = m_pA12;
    m_pSavedA13 = m_pA13;
    m_pSavedA22 = m_pA22;
    m_pSavedA23 = m_pA23;
    m_pSavedGamma11 = m_pGamma11;
    m_pSavedGamma12 = m_pGamma12;
    m_pSavedGamma13 = m_pGamma13;
    m_pSavedGamma22 = m_pGamma22;
    m_pSavedGamma23 = m_pGamma23;
    m_pSavedG = m_pG;
    m_pSavedMomentumTable = m_pMomentumTable;
    m_pSavedTempFFTBuffer = m_pTempFFTBuffer;
    m_iSavedLstDims[0] = m_lstDims.Num() > 0 ? m_lstDims[0] : 0;
    m_iSavedLstDims[1] = m_lstDims.Num() > 1 ? m_lstDims[1] : 0;
    m_iSavedLstDims[2] = m_lstDims.Num() > 2 ? m_lstDims[2] : 0;

    const UINT uiVol = _HC_Volume_xyz;
    const UINT uiDirMinus1 = static_cast<UINT>(_HC_Dir - 1);
    checkCudaErrors(__cudaMalloc((void**)&m_pA11, uiVol * uiDirMinus1 * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA12, uiVol * uiDirMinus1 * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA13, uiVol * uiDirMinus1 * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA22, uiVol * uiDirMinus1 * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA23, uiVol * uiDirMinus1 * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma11, uiVol * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma12, uiVol * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma13, uiVol * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma22, uiVol * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma23, uiVol * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pG, uiVol * sizeof(deviceSU3)));
    if (m_bFA)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pMomentumTable, uiVol * sizeof(DOUBLE)));
        checkCudaErrors(__cudaMalloc((void**)&m_pTempFFTBuffer, uiVol * sizeof(cuDoubleComplex)));

        m_lstDims.RemoveAll();
        m_lstDims.AddItem(static_cast<INT>(_HC_Lx));
        m_lstDims.AddItem(static_cast<INT>(_HC_Ly));
        m_lstDims.AddItem(static_cast<INT>(_HC_Lz));

        preparethread_S;
        _LAUNCH_KERNEL(_kernelBakeMomentumTable3D, block3d, threads3d, m_pMomentumTable, uiVol);
    }
}

void CGaugeFixingCoulombCornell::RestoreLocalBuffers()
{
    //Free the temporary global buffers and put the local ones back (plus the
    //local dims). The destructor frees whatever is current, so it must always
    //see the local pointers afterwards.
    checkCudaErrors(__cudaFree(m_pA11));
    checkCudaErrors(__cudaFree(m_pA12));
    checkCudaErrors(__cudaFree(m_pA13));
    checkCudaErrors(__cudaFree(m_pA22));
    checkCudaErrors(__cudaFree(m_pA23));
    checkCudaErrors(__cudaFree(m_pGamma11));
    checkCudaErrors(__cudaFree(m_pGamma12));
    checkCudaErrors(__cudaFree(m_pGamma13));
    checkCudaErrors(__cudaFree(m_pGamma22));
    checkCudaErrors(__cudaFree(m_pGamma23));
    checkCudaErrors(__cudaFree(m_pG));
    checkCudaErrors(__cudaFree(m_pMomentumTable));
    checkCudaErrors(__cudaFree(m_pTempFFTBuffer));
    m_pA11 = m_pSavedA11;
    m_pA12 = m_pSavedA12;
    m_pA13 = m_pSavedA13;
    m_pA22 = m_pSavedA22;
    m_pA23 = m_pSavedA23;
    m_pGamma11 = m_pSavedGamma11;
    m_pGamma12 = m_pSavedGamma12;
    m_pGamma13 = m_pSavedGamma13;
    m_pGamma22 = m_pSavedGamma22;
    m_pGamma23 = m_pSavedGamma23;
    m_pG = m_pSavedG;
    m_pMomentumTable = m_pSavedMomentumTable;
    m_pTempFFTBuffer = m_pSavedTempFFTBuffer;
    m_lstDims.RemoveAll();
    m_lstDims.AddItem(m_iSavedLstDims[0]);
    m_lstDims.AddItem(m_iSavedLstDims[1]);
    m_lstDims.AddItem(m_iSavedLstDims[2]);
    m_pSavedA11 = NULL;
    m_pSavedA12 = NULL;
    m_pSavedA13 = NULL;
    m_pSavedA22 = NULL;
    m_pSavedA23 = NULL;
    m_pSavedGamma11 = NULL;
    m_pSavedGamma12 = NULL;
    m_pSavedGamma13 = NULL;
    m_pSavedGamma22 = NULL;
    m_pSavedGamma23 = NULL;
    m_pSavedG = NULL;
    m_pSavedMomentumTable = NULL;
    m_pSavedTempFFTBuffer = NULL;
}
#endif

CCString CGaugeFixingCoulombCornell::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingCoulombCornell\n");
    sRet = sRet + tab + _T("accuray : ") + appToString(m_fAccuracy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================