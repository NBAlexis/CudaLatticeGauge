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
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CGaugeFixingLandauCornell.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Cornell Steepest Descend

/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateA(
    const deviceSU3* __restrict__ pU,
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23,
    BYTE byFieldId)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex] = 0.0;
            pA12[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex] = 0.0;
            pA23[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
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
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23,
    BYTE byFieldId)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A = su3A.Log();
            pA11[uiLinkIndex] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex] = 0.0;
            pA12[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex] = 0.0;
            pA23[uiLinkIndex] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

/**
 * Gamma(n) = Delta _{-mu} A(n) = \sum _mu (A_mu(n - mu) - A_mu(n))
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAGradient(
    BYTE byFieldId,
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
    intokernalInt4;

    pGamma11[uiSiteIndex] = 0.0;
    pGamma12[uiSiteIndex] = make_cuDoubleComplex(0.0, 0.0);
    pGamma13[uiSiteIndex] = make_cuDoubleComplex(0.0, 0.0);
    pGamma22[uiSiteIndex] = 0.0;
    pGamma23[uiSiteIndex] = make_cuDoubleComplex(0.0, 0.0);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet())
    {
        return;
    }

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            pGamma11[uiSiteIndex] = pGamma11[uiSiteIndex] - pA11[uiLinkIndex];
            pGamma12[uiSiteIndex] = cuCsub(pGamma12[uiSiteIndex], pA12[uiLinkIndex]);
            pGamma13[uiSiteIndex] = cuCsub(pGamma13[uiSiteIndex], pA13[uiLinkIndex]);
            pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] - pA22[uiLinkIndex];
            pGamma23[uiSiteIndex] = cuCsub(pGamma23[uiSiteIndex], pA23[uiLinkIndex]);
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
                pGamma12[uiSiteIndex] = cuCsub(pGamma12[uiSiteIndex], pA12[uiLinkIndex2]);
                pGamma13[uiSiteIndex] = cuCsub(pGamma13[uiSiteIndex], pA13[uiLinkIndex2]);
                pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] - pA22[uiLinkIndex2];
                pGamma23[uiSiteIndex] = cuCsub(pGamma23[uiSiteIndex], pA23[uiLinkIndex2]);
            }
            else
            {
                pGamma11[uiSiteIndex] = pGamma11[uiSiteIndex] + pA11[uiLinkIndex2];
                pGamma12[uiSiteIndex] = cuCadd(pGamma12[uiSiteIndex], pA12[uiLinkIndex2]);
                pGamma13[uiSiteIndex] = cuCadd(pGamma13[uiSiteIndex], pA13[uiLinkIndex2]);
                pGamma22[uiSiteIndex] = pGamma22[uiSiteIndex] + pA22[uiLinkIndex2];
                pGamma23[uiSiteIndex] = cuCadd(pGamma23[uiSiteIndex], pA23[uiLinkIndex2]);
            }
        }
    }
}

/**
 * g(x)=exp(-i a Delta A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateG(
    BYTE byFieldId,
    deviceSU3* pG,
    const DOUBLE* __restrict__ pGamma11,
    const cuDoubleComplex* __restrict__ pGamma12,
    const cuDoubleComplex* __restrict__ pGamma13,
    const DOUBLE* __restrict__ pGamma22,
    const cuDoubleComplex* __restrict__ pGamma23,
    DOUBLE fAlpha
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pG[uiSiteIndex] = deviceSU3::makeSU3Id();
    }
    else
    {
        deviceSU3 pA = deviceSU3::makeSU3TA(
#if !_CLG_DOUBLEFLOAT
            _cToFloat(pGamma12[uiSiteIndex]), _cToFloat(pGamma13[uiSiteIndex]), _cToFloat(pGamma23[uiSiteIndex]),
            static_cast<Real>(pGamma11[uiSiteIndex]), static_cast<Real>(pGamma22[uiSiteIndex])
#else
            pGamma12[uiSiteIndex], pGamma13[uiSiteIndex], pGamma23[uiSiteIndex],
            pGamma11[uiSiteIndex], pGamma22[uiSiteIndex]
#endif
        );
        pG[uiSiteIndex] = _expreal(pA, fAlpha);
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
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
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
    BYTE byFieldId,
    DOUBLE* pDeviceRes,
    const DOUBLE* __restrict__ pDeltaA11,
    const cuDoubleComplex* __restrict__ pDeltaA12,
    const cuDoubleComplex* __restrict__ pDeltaA13,
    const DOUBLE* __restrict__ pDeltaA22,
    const cuDoubleComplex* __restrict__ pDeltaA23
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex] = 0.0;

        return;
    }

    DOUBLE fAbs1 = cuCabs(pDeltaA12[uiSiteIndex]);
    DOUBLE fAbs2 = cuCabs(pDeltaA13[uiSiteIndex]);
    DOUBLE fAbs3 = cuCabs(pDeltaA23[uiSiteIndex]);
    DOUBLE fM1122 = pDeltaA11[uiSiteIndex] + pDeltaA22[uiSiteIndex];
    pDeviceRes[uiSiteIndex] = 2.0 * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}


#pragma endregion

#pragma region FFT accelaration

__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable(
    DOUBLE* pP,
    UINT uiV)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);

    DOUBLE fDenorm = static_cast<DOUBLE>(uiDir);
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        fDenorm -= cos(2.0 * PI * sSite4.m_byData4[dir] / static_cast<DOUBLE>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = 0.5;
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 4 / 2(Nd - sum cos)
    pP[uiSiteIndex] = 8.0 / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC(const DOUBLE* __restrict__ realBuffer, cuDoubleComplex* complexBuffer)
{
    intokernal;
    complexBuffer[uiSiteIndex] = make_cuDoubleComplex(realBuffer[uiSiteIndex], 0.0);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR(const cuDoubleComplex* __restrict__ complexBuffer, DOUBLE* realBuffer)
{
    intokernal;
    realBuffer[uiSiteIndex] = complexBuffer[uiSiteIndex].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale(const DOUBLE* __restrict__ pP, cuDoubleComplex* fftRes)
{
    intokernal;
    fftRes[uiSiteIndex] = cuCmulf_cd(fftRes[uiSiteIndex], pP[uiSiteIndex]);
}

#pragma endregion

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingLandauCornell)

void CGaugeFixingLandauCornell::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    //========== Initial Settings ==============
#if !_CLG_DOUBLEFLOAT
    if (!params.FetchValueDOUBLE(_T("Alpha"), m_fAlpha))
#else
    if (!params.FetchValueReal(_T("Alpha"), m_fAlpha))
#endif
    {
        appGeneral(_T("CGaugeFixingLandauCornell: Alpha not set, set to 0.08 by defualt."));
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: Accuracy not set, set to 0.00000000001 by defualt."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
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
#if !_CLG_DOUBLEFLOAT
    checkCudaErrors(__cudaMalloc((void**)&m_pA11, _HC_Volume * _HC_Dir * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA12, _HC_Volume * _HC_Dir * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA13, _HC_Volume * _HC_Dir * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA22, _HC_Volume * _HC_Dir * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA23, _HC_Volume * _HC_Dir * sizeof(cuDoubleComplex)));

    checkCudaErrors(__cudaMalloc((void**)&m_pGamma11, _HC_Volume * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma12, _HC_Volume * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma13, _HC_Volume * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma22, _HC_Volume * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pGamma23, _HC_Volume * sizeof(cuDoubleComplex)));
#else
    checkCudaErrors(__cudaMalloc((void**)& m_pA11, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(__cudaMalloc((void**)& m_pA12, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(__cudaMalloc((void**)& m_pA13, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(__cudaMalloc((void**)& m_pA22, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(__cudaMalloc((void**)& m_pA23, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));

    checkCudaErrors(__cudaMalloc((void**)& m_pGamma11, _HC_Volume * sizeof(Real)));
    checkCudaErrors(__cudaMalloc((void**)& m_pGamma12, _HC_Volume * sizeof(CLGComplex)));
    checkCudaErrors(__cudaMalloc((void**)& m_pGamma13, _HC_Volume * sizeof(CLGComplex)));
    checkCudaErrors(__cudaMalloc((void**)& m_pGamma22, _HC_Volume * sizeof(Real)));
    checkCudaErrors(__cudaMalloc((void**)& m_pGamma23, _HC_Volume * sizeof(CLGComplex)));
#endif
    checkCudaErrors(__cudaMalloc((void**)& m_pG, _HC_Volume * sizeof(deviceSU3)));
    if (m_bFA)
    {
#if !_CLG_DOUBLEFLOAT
        checkCudaErrors(__cudaMalloc((void**)&m_pMomentumTable, _HC_Volume * sizeof(DOUBLE)));
        checkCudaErrors(__cudaMalloc((void**)&m_pTempFFTBuffer, _HC_Volume * sizeof(cuDoubleComplex)));
#else
        checkCudaErrors(__cudaMalloc((void**)& m_pMomentumTable, _HC_Volume * sizeof(Real)));
        checkCudaErrors(__cudaMalloc((void**)& m_pTempFFTBuffer, _HC_Volume * sizeof(CLGComplex)));
#endif

        preparethread;
        _LAUNCH_KERNEL(_kernelBakeMomentumTable, block, threads, m_pMomentumTable, _HC_Volume);
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

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.4: same gather -> rank0 global-context fix -> scatter flow as the
        //Coulomb fixers (P4-2.2/P4-2.3). The iteration loop runs on the gathered
        //global buffer under the temporary GLOBAL lattice context; the 4D cuFFT
        //plans are re-created per call with the global _HC_Lx/y/z/t constants.
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
            GaugeFixingLoop(pDevGlobal, pResGauge->m_byFieldId);
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

    GaugeFixingLoop(pDeviceBufferPointer, pResGauge->m_byFieldId);
}

void CGaugeFixingLandauCornell::GaugeFixingLoop(deviceSU3* pDeviceBufferPointer, BYTE byFieldId)
{
    preparethread;
    m_iIterate = 0;
#if !_CLG_DOUBLEFLOAT
    DOUBLE fTheta = 0.0;
#else
    Real fTheta = F(0.0);
#endif

    while (m_iIterate < m_iMaxIterate)
    {
        //======= 1. Calculate Gamma    =========
        if (0 == _HC_ALog)
        {
            _LAUNCH_KERNEL(_kernelCalculateA, block, threads, 
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
            _LAUNCH_KERNEL(_kernelCalculateALog, block, threads, 
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23,
                byFieldId);
        }

        _LAUNCH_KERNEL(_kernelCalculateAGradient, block, threads, 
            byFieldId,
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
        _LAUNCH_KERNEL(_kernelCalculateTrAGradientSq, block, threads, 
            byFieldId,
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
#if !_CLG_DOUBLEFLOAT

            CCLGFFTHelper::FFT4DDouble(m_pGamma12, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::FFT4DDouble(m_pGamma12, FALSE);

            CCLGFFTHelper::FFT4DDouble(m_pGamma13, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::FFT4DDouble(m_pGamma13, FALSE);

            CCLGFFTHelper::FFT4DDouble(m_pGamma23, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::FFT4DDouble(m_pGamma23, FALSE);

            _LAUNCH_KERNEL(_kernelFFTRtoC, block, threads, m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4DDouble(m_pTempFFTBuffer, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4DDouble(m_pTempFFTBuffer, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR, block, threads, m_pTempFFTBuffer, m_pGamma11);

            _LAUNCH_KERNEL(_kernelFFTRtoC, block, threads, m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4DDouble(m_pTempFFTBuffer, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4DDouble(m_pTempFFTBuffer, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR, block, threads, m_pTempFFTBuffer, m_pGamma22);

#else
            CCLGFFTHelper::FFT4D(m_pGamma12, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::FFT4D(m_pGamma12, FALSE);

            CCLGFFTHelper::FFT4D(m_pGamma13, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::FFT4D(m_pGamma13, FALSE);

            CCLGFFTHelper::FFT4D(m_pGamma23, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::FFT4D(m_pGamma23, FALSE);

            _LAUNCH_KERNEL(_kernelFFTRtoC, block, threads, m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR, block, threads, m_pTempFFTBuffer, m_pGamma11);

            _LAUNCH_KERNEL(_kernelFFTRtoC, block, threads, m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, TRUE);
            _LAUNCH_KERNEL(_kernelFFTScale, block, threads, m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT4D(m_pTempFFTBuffer, FALSE);
            _LAUNCH_KERNEL(_kernelFFTCtoR, block, threads, m_pTempFFTBuffer, m_pGamma22);
#endif
        }

        //======= 4. Gauge Transform    =========
        //Be careful not to use strictLog when A is really small
        _LAUNCH_KERNEL(_kernelCalculateG, block, threads, 
            byFieldId,
            m_pG,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_fAlpha);

        _LAUNCH_KERNEL(_kernelGaugeTransform, block, threads, byFieldId, m_pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("Gauge fixing failed with last error = %f\n"), fTheta);
}

#if _CLG_MULTI_GPU
void CGaugeFixingLandauCornell::ResizeBuffersToGlobal()
{
    //Called under the temporary GLOBAL lattice context (rank 0 only), where
    //_HC_Volume / _HC_Dir are the GLOBAL values: save the local pointers,
    //re-allocate the 4D-volume fixing buffers to the global volume and re-bake
    //the momentum table (feeds FFT4D, whose plans are re-created per call from
    //the runtime _HC_Lx/y/z/t constants).
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
    const UINT uiVol = _HC_Volume;
    const UINT uiDir = static_cast<UINT>(_HC_Dir);
    checkCudaErrors(__cudaMalloc((void**)&m_pA11, uiVol * uiDir * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA12, uiVol * uiDir * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA13, uiVol * uiDir * sizeof(cuDoubleComplex)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA22, uiVol * uiDir * sizeof(DOUBLE)));
    checkCudaErrors(__cudaMalloc((void**)&m_pA23, uiVol * uiDir * sizeof(cuDoubleComplex)));
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

        preparethread;
        _LAUNCH_KERNEL(_kernelBakeMomentumTable, block, threads, m_pMomentumTable, uiVol);
    }
}

void CGaugeFixingLandauCornell::RestoreLocalBuffers()
{
    //Free the temporary global buffers and put the local ones back. The
    //destructor frees whatever is current, so it must always see the local
    //pointers afterwards.
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

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingLandauCornell::CheckRes(const CFieldGauge* pGauge)
#else
Real CGaugeFixingLandauCornell::CheckRes(const CFieldGauge* pGauge)
#endif
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
#if !_CLG_DOUBLEFLOAT
        return 0.0;
#else
        return F(0.0);
#endif
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.4: the deviation must be measured on the WHOLE lattice; same flow as
        //the Coulomb fixers (P4-2.2/P4-2.3). Non-root ranks report their local
        //deviation; the driver only compares rank 0.
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
#if !_CLG_DOUBLEFLOAT
        const DOUBLE fRet = CheckResLocal(pDevGlobal, pGauge->m_byFieldId);
#else
        const Real fRet = CheckResLocal(pDevGlobal, pGauge->m_byFieldId);
#endif
        RestoreLocalBuffers();
        MGExitGlobalFixerContext();

        checkCudaErrors(__cudaFree(pDevGlobal));
        return fRet;
    }
#endif
    return CheckResLocal(pGaugeSU3->m_pDeviceData, pGauge->m_byFieldId);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingLandauCornell::CheckResLocal(const deviceSU3* pGaugeData, BYTE byFieldId)
#else
Real CGaugeFixingLandauCornell::CheckResLocal(const deviceSU3* pGaugeData, BYTE byFieldId)
#endif
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _LAUNCH_KERNEL(_kernelCalculateA, block, threads, 
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
        _LAUNCH_KERNEL(_kernelCalculateALog, block, threads, 
            pGaugeData,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23,
            byFieldId);
    }

    _LAUNCH_KERNEL(_kernelCalculateAGradient, block, threads, 
        byFieldId,
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

    _LAUNCH_KERNEL(_kernelCalculateTrAGradientSq, block, threads, 
        byFieldId,
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
    sRet = sRet + tab + _T("accuray : ") + appToString(m_fAccuracy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================