//=============================================================================
// FILENAME : CGaugeFixingLandauLosAlamos.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/23/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGOdd_S(
    SBYTE uiT,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const deviceSU3 su3Id = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || !sSite4.IsOdd())
    {
        return;
    }

    pG[uiSiteIndex3D] = deviceSU3::makeSU3Zero();
    
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        const UINT p_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + dir];
        const SIndex site_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_mu];
        const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_mu.m_uiSiteIndex, dir);

        pG[uiSiteIndex3D].Add(__idx->_deviceIsBondOnSurface(uiBigIdx, dir) ? su3Id : pU[uiLinkIndex]);
        pG[uiSiteIndex3D].AddDagger(__idx->_deviceIsBondOnSurface(p_m_mu, dir) ? su3Id : pU[uiLinkIndex2]);
    }

    pG[uiSiteIndex3D].MulReal(fOmega);
    pG[uiSiteIndex3D].Add(su3Id.MulRealC(F(1.0) - fOmega));
    pG[uiSiteIndex3D].CabbiboMarinariProj();
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGEven_S(
    SBYTE uiT,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const deviceSU3 su3Id = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || sSite4.IsOdd())
    {
        return;
    }

    pG[uiSiteIndex3D] = deviceSU3::makeSU3Zero();

    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        const UINT p_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + dir];
        const SIndex site_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_mu];
        const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_mu.m_uiSiteIndex, dir);

        pG[uiSiteIndex3D].Add(__idx->_deviceIsBondOnSurface(uiBigIdx, dir) ? su3Id : pU[uiLinkIndex]);
        pG[uiSiteIndex3D].AddDagger(__idx->_deviceIsBondOnSurface(p_m_mu, dir) ? su3Id : pU[uiLinkIndex2]);
    }

    pG[uiSiteIndex3D].MulReal(fOmega);
    pG[uiSiteIndex3D].Add(su3Id.MulRealC(F(1.0) - fOmega));
    pG[uiSiteIndex3D].CabbiboMarinariProj();
}

/**
 * U_mu(n odd) = g(n) U_mu(n odd)
 * U_mu(n even) =U_mu(n even) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformOdd_S(
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir - 1; ++dir)
            {
                const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                pGauge[uiLinkDir] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkDir]);
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
            {
                const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const UINT p_p_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + dir + uiDir];
                const SIndex site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_p_mu];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex / _DC_Lt]);
                }
            }
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformEven_S(
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (!sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir - 1; ++dir)
            {
                const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                pGauge[uiLinkDir] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkDir]);
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
            {
                const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const UINT p_p_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + dir + uiDir];
                const SIndex site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_p_mu];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex / _DC_Lt]);
                }
            }
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3DTOdd(
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    if (!sSite4.IsOdd())
    {
        return;
    }

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
        if (!site.IsDirichlet())
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            pGauge[uiLinkIndex] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex]);
        }
    }

    const UINT p_m_t = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 3];

    if (!__idx->_deviceIsBondOnSurface(p_m_t, 3))
    {
        const SIndex site_m_t = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_t];
        if (!site_m_t.IsDirichlet())
        {
            const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_t.m_uiSiteIndex, 3);
            pGauge[uiLinkIndex2].MulDagger(pGx[uiSiteIndex3D]);
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3DTEven(
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    if (sSite4.IsOdd())
    {
        return;
    }

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
        if (!site.IsDirichlet())
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            pGauge[uiLinkIndex] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex]);
        }
    }

    const UINT p_m_t = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + 3];

    if (!__idx->_deviceIsBondOnSurface(p_m_t, 3))
    {
        const SIndex site_m_t = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_t];
        if (!site_m_t.IsDirichlet())
        {
            const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_t.m_uiSiteIndex, 3);
            pGauge[uiLinkIndex2].MulDagger(pGx[uiSiteIndex3D]);
        }
    }
}

#pragma region Commen

/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateASpace_S(
        SBYTE uiT,
        const deviceSU3* __restrict__ pU,
        Real* pA11,
        CLGComplex* pA12,
        CLGComplex* pA13,
        Real* pA22,
        CLGComplex* pA23)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex3D] = su3A.m_me[0].y;
            pA12[uiLinkIndex3D] = su3A.m_me[1];
            pA13[uiLinkIndex3D] = su3A.m_me[2];
            pA22[uiLinkIndex3D] = su3A.m_me[4].y;
            pA23[uiLinkIndex3D] = su3A.m_me[5];
        }
        else
        {
            pA11[uiLinkIndex3D] = F(0.0);
            pA12[uiLinkIndex3D] = _zeroc;
            pA13[uiLinkIndex3D] = _zeroc;
            pA22[uiLinkIndex3D] = F(0.0);
            pA23[uiLinkIndex3D] = _zeroc;
        }
    }
}

/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateCoulombDivation_S(
    SBYTE uiT,
    Real* pDeviceRes,
    const Real* __restrict__ pA11,
    const CLGComplex* __restrict__ pA12,
    const CLGComplex* __restrict__ pA13,
    const Real* __restrict__ pA22,
    const CLGComplex* __restrict__ pA23)
{
    intokernalInt4_S_Only3D;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex3D] = F(0.0);
        return;
    }

    Real g11 = F(0.0);
    CLGComplex g12 = _zeroc;
    CLGComplex g13 = _zeroc;
    Real g22 = F(0.0);
    CLGComplex g23 = _zeroc;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const BYTE uiDir2 = uiDir * 2;

    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            g11 = g11 - pA11[uiLinkIndex3D];
            g12 = _cuCsubf(g12, pA12[uiLinkIndex3D]);
            g13 = _cuCsubf(g13, pA13[uiLinkIndex3D]);
            g22 = g22 - pA22[uiLinkIndex3D];
            g23 = _cuCsubf(g23, pA23[uiLinkIndex3D]);
        }

        const UINT p_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + dir];
        const SIndex site_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][p_m_mu];
        const UINT uiLinkIndex2_3D = (site_m_mu.m_uiSiteIndex / _DC_Lt) * (uiDir - 1) + dir;

        if (!__idx->_deviceIsBondOnSurface(p_m_mu, dir))
        {
            g11 = g11 + pA11[uiLinkIndex2_3D];
            g12 = _cuCaddf(g12, pA12[uiLinkIndex2_3D]);
            g13 = _cuCaddf(g13, pA13[uiLinkIndex2_3D]);
            g22 = g22 + pA22[uiLinkIndex2_3D];
            g23 = _cuCaddf(g23, pA23[uiLinkIndex2_3D]);
        }
    }

    //gamma now is sum _{mu} Delta _{mu} A _{mu} (or, partial A)
    //then calculate Tr[g g^dagger]
    const Real fAbs1 = _cuCabsf(g12);
    const Real fAbs2 = _cuCabsf(g13);
    const Real fAbs3 = _cuCabsf(g23);
    const Real fM1122 = g11 + g22;
    pDeviceRes[uiSiteIndex3D] = F(2.0) * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}


#pragma endregion


#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingCoulombLosAlamos)

void CGaugeFixingCoulombLosAlamos::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    const TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx);
    latticeDim.AddItem(_HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    m_pHDecomp[0] = decomp[0];
    m_pHDecomp[1] = decomp[1];
    m_pHDecomp[2] = decomp[2];
    m_pHDecomp[3] = decomp[3];
    m_pHDecomp[4] = decomp[4];
    m_pHDecomp[5] = decomp[5];
    m_lstDims.RemoveAll();
    m_lstDims.AddItem(_HC_Lx);
    m_lstDims.AddItem(_HC_Ly);
    m_lstDims.AddItem(_HC_Lz);

    checkCudaErrors(cudaMalloc((void**)& m_pDDecomp, sizeof(UINT) * 6));

    //========== Initial Settings ==============
    if (!params.FetchValueReal(_T("Omega"), m_fOmega))
    {
        appGeneral(_T("CGaugeFixingCoulombLosAlamos: Omega not set, set to 1.0 by defualt."));
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingCoulombLosAlamos: Accuracy not set, set to 0.00000000001 by defualt."));
    }

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombLosAlamos: MaxIterate not set, set to 100000 by defualt."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = static_cast<INT>(m_iCheckErrorStep);
    if (!params.FetchValueINT(_T("CheckErrorStep"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombLosAlamos: CheckErrorStep not set, set to 1000 by defualt."));
    }
    m_iCheckErrorStep = static_cast<UINT>(iValue);

    //========== Initial Buffers ==============
    checkCudaErrors(cudaMalloc((void**)& m_pG, _HC_Volume_xyz * sizeof(deviceSU3)));
    checkCudaErrors(cudaMalloc((void**)& m_pA11, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA12, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA13, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA22, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA23, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));
}

Real CGaugeFixingCoulombLosAlamos::CheckRes(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return F(0.0);
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    return CheckResDeviceBuffer(pGaugeSU3->m_pDeviceData);
}

Real CGaugeFixingCoulombLosAlamos::CheckResDeviceBuffer(const deviceSU3* __restrict__ pGauge)
{
    Real fRes = F(0.0);
    preparethread_S;
    for (SBYTE uiT = 0; uiT < static_cast<SBYTE>(_HC_Lt); ++uiT)
    {
        _kernelCalculateASpace_S << <block, threads >> > (
            uiT,
            pGauge,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);

        _kernelCalculateCoulombDivation_S << <block, threads >> > (
            uiT,
            _D_RealThreadBuffer,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);
        fRes += appAbs(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz));
    }
    return fRes;
}

Real CGaugeFixingCoulombLosAlamos::CheckResDeviceBufferOnlyT(const deviceSU3* __restrict__ pGauge, SBYTE uiT)
{
    preparethread_S;

    _kernelCalculateASpace_S << <block, threads >> > (
        uiT,
        pGauge,
        m_pA11,
        m_pA12,
        m_pA13,
        m_pA22,
        m_pA23);

    _kernelCalculateCoulombDivation_S << <block, threads >> > (
        uiT,
        _D_RealThreadBuffer,
        m_pA11,
        m_pA12,
        m_pA13,
        m_pA22,
        m_pA23);
    return appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz);
}

void CGaugeFixingCoulombLosAlamos::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge || EFT_GaugeSU3 != pResGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingCoulombLosAlamos only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
    deviceSU3* pDeviceBufferPointer = pGaugeSU3->m_pDeviceData;

    for (SBYTE uiT = 0; uiT < static_cast<SBYTE>(_HC_Lt); ++uiT)
    {
        GaugeFixingForT(pDeviceBufferPointer, uiT);
    }

    //appGeneral(_T("Gauge fixing failed with last error = %f\n"), fTheta);
}

void CGaugeFixingCoulombLosAlamos::GaugeFixingForT(deviceSU3* pDeviceBufferPointer, SBYTE uiT)
{
    preparethread_S;
    m_iIterate = 0;

    while (m_iIterate < m_iMaxIterate)
    {
        //check res
        if (0 == m_iIterate % m_iCheckErrorStep)
        {
            const Real fTheta = CheckResDeviceBufferOnlyT(pDeviceBufferPointer, uiT);
            appDetailed(_T("Iterate : %d, error = %2.12f\n"), m_iIterate, fTheta);
            if (fTheta < m_fAccuracy)
            {
                return;
            }
        }

        _kernelCalculateGOdd_S << <block, threads >> > (uiT, pDeviceBufferPointer, m_fOmega, m_pG);
        _kernelGaugeTransformOdd_S << <block, threads >> > (uiT, m_pG, pDeviceBufferPointer);
        _kernelGaugeTransform3DTOdd << <block, threads >> > (uiT, m_pG, pDeviceBufferPointer);

        _kernelCalculateGEven_S << <block, threads >> > (uiT, pDeviceBufferPointer, m_fOmega, m_pG);
        _kernelGaugeTransformEven_S << <block, threads >> > (uiT, m_pG, pDeviceBufferPointer);
        _kernelGaugeTransform3DTEven << <block, threads >> > (uiT, m_pG, pDeviceBufferPointer);

        ++m_iIterate;
    }
}

CCString CGaugeFixingCoulombLosAlamos::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingCoulombLosAlamos\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================