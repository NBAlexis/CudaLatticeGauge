//=============================================================================
// FILENAME : CGaugeFixingLandauLosAlamos.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/21/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGAL(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const deviceSU3 su3Id = deviceSU3::makeSU3Id();
    if (site.IsDirichlet())
    {
        pG[uiSiteIndex] = su3Id;
        return;
    }

    pG[uiSiteIndex] = deviceSU3::makeSU3Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];

        pG[uiSiteIndex].Add(__idx->_deviceIsBondOnSurface(uiBigIdx, dir) ? su3Id : pU[uiLinkIndex]);
        pG[uiSiteIndex].AddDagger(_deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu));
    }

    pG[uiSiteIndex].MulReal(fOmega);
    pG[uiSiteIndex].Add(su3Id.MulRealC(F(1.0) - fOmega));
    pG[uiSiteIndex].CabbiboMarinariProj();
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGOdd(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const deviceSU3 su3Id = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || !sSite4.IsOdd())
    {
        //pG[uiSiteIndex] = su3Id;
        return;
    }

    pG[uiSiteIndex] = deviceSU3::makeSU3Zero();
    
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];

        pG[uiSiteIndex].Add(__idx->_deviceIsBondOnSurface(uiBigIdx, dir) ? su3Id : pU[uiLinkIndex]);
        pG[uiSiteIndex].AddDagger(_deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu));
    }

    pG[uiSiteIndex].MulReal(fOmega);
    pG[uiSiteIndex].Add(su3Id.MulRealC(F(1.0) - fOmega));
    pG[uiSiteIndex].CabbiboMarinariProj();
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGEven(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const deviceSU3 su3Id = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || sSite4.IsOdd())
    {
        //pG[uiSiteIndex] = su3Id;
        return;
    }

    pG[uiSiteIndex] = deviceSU3::makeSU3Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];

        pG[uiSiteIndex].Add(__idx->_deviceIsBondOnSurface(uiBigIdx, dir) ? su3Id : pU[uiLinkIndex]);
        pG[uiSiteIndex].AddDagger(_deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu));
    }

    pG[uiSiteIndex].MulReal(fOmega);
    pG[uiSiteIndex].Add(su3Id.MulRealC(F(1.0) - fOmega));
    pG[uiSiteIndex].CabbiboMarinariProj();
}

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformAL(
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
                pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

/**
 * U_mu(n odd) = g(n) U_mu(n odd)
 * U_mu(n even) =U_mu(n even) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformOdd(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir; ++dir)
            {
                //If site is not Dirichlet, the bound should not be on the surface
                //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
                //{
                    UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                    pGauge[uiLinkDir] = pGx[uiSiteIndex].MulC(pGauge[uiLinkDir]);
                //}
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
                const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_p_mu_site)];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
                }
            }
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformEven(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (!sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir; ++dir)
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                pGauge[uiLinkDir] = pGx[uiSiteIndex].MulC(pGauge[uiLinkDir]);
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
                const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_p_mu_site)];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
                }
            }
        }
    }
}

#pragma region Commen

/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAAll(
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

/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateLandauDivation(
    BYTE byFieldId,
    Real* pDeviceRes,
    const Real* __restrict__ pA11,
    const CLGComplex* __restrict__ pA12,
    const CLGComplex* __restrict__ pA13,
    const Real* __restrict__ pA22,
    const CLGComplex* __restrict__ pA23)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex] = F(0.0);
        return;
    }

    Real g11 = F(0.0);
    CLGComplex g12 = _zeroc;
    CLGComplex g13 = _zeroc;
    Real g22 = F(0.0);
    CLGComplex g23 = _zeroc;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            g11 = g11 - pA11[uiLinkIndex];
            g12 = _cuCsubf(g12, pA12[uiLinkIndex]);
            g13 = _cuCsubf(g13, pA13[uiLinkIndex]);
            g22 = g22 - pA22[uiLinkIndex];
            g23 = _cuCsubf(g23, pA23[uiLinkIndex]);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_mu.m_uiSiteIndex, dir);

        if (!site_m_mu.IsDirichlet())
        {
            if (site_m_mu.NeedToDagger())
            {
                g11 = g11 - pA11[uiLinkIndex2];
                g12 = _cuCsubf(g12, pA12[uiLinkIndex2]);
                g13 = _cuCsubf(g13, pA13[uiLinkIndex2]);
                g22 = g22 - pA22[uiLinkIndex2];
                g23 = _cuCsubf(g23, pA23[uiLinkIndex2]);
            }
            else
            {
                g11 = g11 + pA11[uiLinkIndex2];
                g12 = _cuCaddf(g12, pA12[uiLinkIndex2]);
                g13 = _cuCaddf(g13, pA13[uiLinkIndex2]);
                g22 = g22 + pA22[uiLinkIndex2];
                g23 = _cuCaddf(g23, pA23[uiLinkIndex2]);
            }
        }
    }

    //gamma now is sum _{mu} Delta _{mu} A _{mu} (or, partial A)
    //then calculate Tr[g g^dagger]
    const Real fAbs1 = _cuCabsf(g12);
    const Real fAbs2 = _cuCabsf(g13);
    const Real fAbs3 = _cuCabsf(g23);
    const Real fM1122 = g11 + g22;
    pDeviceRes[uiSiteIndex] = F(2.0) * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}


#pragma endregion


#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingLandauLosAlamos)

void CGaugeFixingLandauLosAlamos::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    //========== Initial Settings ==============
    if (!params.FetchValueReal(_T("Omega"), m_fOmega))
    {
        appGeneral(_T("CGaugeFixingLandauLosAlamos: Omega not set, set to 1.0 by defualt."));
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

    iValue = static_cast<INT>(m_iCheckErrorStep);
    if (!params.FetchValueINT(_T("CheckErrorStep"), iValue))
    {
        appGeneral(_T("CGaugeFixingLandauCornell: CheckErrorStep not set, set to 1000 by defualt."));
    }
    m_iCheckErrorStep = static_cast<UINT>(iValue);

    //========== Initial Buffers ==============
    checkCudaErrors(cudaMalloc((void**)& m_pG, _HC_Volume * sizeof(deviceSU3)));
    checkCudaErrors(cudaMalloc((void**)& m_pA11, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA12, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA13, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA22, _HC_Volume * _HC_Dir * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA23, _HC_Volume * _HC_Dir * sizeof(CLGComplex)));
}

Real CGaugeFixingLandauLosAlamos::CheckRes(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return F(0.0);
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    preparethread;
    _kernelCalculateAAll << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        m_pA11,
        m_pA12,
        m_pA13,
        m_pA22,
        m_pA23);

    _kernelCalculateLandauDivation << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        _D_RealThreadBuffer,
        m_pA11,
        m_pA12,
        m_pA13,
        m_pA22,
        m_pA23);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer) / (3 * _HC_Volume);
}

//void CGaugeFixingLandauLosAlamos::TestMeasureFU(const CFieldGauge* pGauge)
//{
//    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
//    {
//        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
//        return;
//    }
//    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
//
//    preparethread;
//    _kernelMeasureF << <block, threads >> > (_D_RealThreadBuffer, pGaugeSU3->m_pDeviceData);
//    const Real fTrace = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer) / (3 * _HC_Volume);
//    appParanoiac(_T("F(U) = %f\n"), fTrace);
//
//}

void CGaugeFixingLandauLosAlamos::GaugeFixing(CFieldGauge* pResGauge)
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
        //check res
        if (0 == m_iIterate % m_iCheckErrorStep)
        {
            fTheta = CheckRes(pGaugeSU3);
            appParanoiac(_T("Iterate : %d, error = %2.12f\n"), m_iIterate, fTheta);
            //TestMeasureFU(pGaugeSU3);
            if (fTheta < m_fAccuracy)
            {
                return;
            }
        }

        _kernelCalculateGOdd << <block, threads >> > (pResGauge->m_byFieldId, pDeviceBufferPointer, m_fOmega, m_pG);
        _kernelGaugeTransformOdd << <block, threads >> > (pResGauge->m_byFieldId, m_pG, pDeviceBufferPointer);
        _kernelCalculateGEven << <block, threads >> > (pResGauge->m_byFieldId, pDeviceBufferPointer, m_fOmega, m_pG);
        _kernelGaugeTransformEven << <block, threads >> > (pResGauge->m_byFieldId, m_pG, pDeviceBufferPointer);

        //_kernelCalculateGAL << <block, threads >> > (pDeviceBufferPointer, m_fOmega, m_pG);
        //_kernelGaugeTransformAL << <block, threads >> > (m_pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("Gauge fixing failed with last error = %f\n"), fTheta);
}

CCString CGaugeFixingLandauLosAlamos::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingLandauLosAlamos\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================