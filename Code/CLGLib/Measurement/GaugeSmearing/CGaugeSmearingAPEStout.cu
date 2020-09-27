//=============================================================================
// FILENAME : CGaugeSmearingAPEStout.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleWithoutT(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3* staple)
{
    intokernal;

    const UINT plaqLengthm1 = plaqLength - 1;     //3
    UINT plaqCountAll = plaqCount * plaqLengthm1; //18

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        //plaqCount = 6, they are:
        //j=0,1,2,3 / mu, for example, when linkIndex is for y, 1, it caches 10,10-,12,12-,13,13-
        
        //If only spatial, we only use the first 4
        #pragma unroll
        for (int i = 0; i < 4; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }

        staple[linkIndex] = res;
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelAPEStoutSU3WithoutT(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fRho)
{
    intokernal;

    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulDagger(gauge[uiLinkIndex]);
        //Note, q = -i (staple*u).TA(), and exp(i q) = exp((staple*u).TA())
        omega.Ta();
        omega = (0 == _DC_ExpPrecision) ?
            omega.QuickExp(fRho)
            : omega.ExpReal(fRho, _DC_ExpPrecision);
        omega.Mul(gauge[uiLinkIndex]);
        gauge[uiLinkIndex] = omega;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAPEStoutSU3(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fRho)
{
    intokernal;

    for (BYTE dir = 0; dir < _DC_Dir; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulDagger(gauge[uiLinkIndex]);
        omega.Ta();
        omega = (0 == _DC_ExpPrecision) ? 
            omega.QuickExp(fRho)
          : omega.ExpReal(fRho, _DC_ExpPrecision);
        omega.Mul(gauge[uiLinkIndex]);
        gauge[uiLinkIndex] = omega;
    }
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeSmearingAPEStout)


void CGaugeSmearingAPEStout::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    if (!params.FetchValueReal(_T("Rho"), m_fRho))
    {
        appCrucial(_T("CGaugeSmearingAPEStout: rhos not set, set to 0.1 by defualt."));
    }

    INT iValue = 0;
    params.FetchValueINT(_T("HasT"), iValue);
    m_bHasT = (0 != iValue);

    iValue = 1;
    if (!params.FetchValueINT(_T("Iterate"), iValue))
    {
        appCrucial(_T("CGaugeSmearingAPEStout: Iterate not set, set to 1 by defualt."));
    }
    else if (iValue > 0)
    {
        m_uiIterate = static_cast<UINT>(iValue);
    }
    else
    {
        appCrucial(_T("CGaugeSmearingAPEStout: Iterate not correct, set to 1 by defualt."));
    }

    appGeneral(_T(" ---- CGaugeSmearingAPEStout created with rho = %f, HasT = %d, iterate = %d\n"), m_fRho, m_bHasT, m_uiIterate);
}

void CGaugeSmearingAPEStout::GaugeSmearing(CFieldGauge* pGauge, CFieldGauge* pStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    if (NULL == pStaple || EFT_GaugeSU3 != pStaple->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pStapleSU3 = dynamic_cast<CFieldGaugeSU3*>(pStaple);
    appParanoiac(_T("GaugeSmearing : CGaugeSmearingAPEStout\n"));

    preparethread;
    for (UINT i = 0; i < m_uiIterate; ++i)
    {
        if (m_bHasT)
        {
            _kernelAPEStoutSU3 << <block, threads >> > (
                pGaugeSU3->m_pDeviceData,
                pStapleSU3->m_pDeviceData,
                m_fRho);

            pGaugeSU3->CalculateOnlyStaple(pStapleSU3);
        }
        else
        {
            _kernelStapleWithoutT << <block, threads >> > (
                pGaugeSU3->m_pDeviceData,
                appGetLattice()->m_pIndexCache->m_pStappleCache,
                appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
                appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
                pStapleSU3->m_pDeviceData);

            _kernelAPEStoutSU3WithoutT << <block, threads >> > (
                pGaugeSU3->m_pDeviceData,
                pStapleSU3->m_pDeviceData,
                m_fRho);
        }
    }

    if (!m_bHasT)
    {
        pGaugeSU3->CalculateOnlyStaple(pStapleSU3);
    }
}

CCString CGaugeSmearingAPEStout::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeSmearingAPEStout\n");
    sRet = sRet + tab + _T("Rhos : ") + appFloatToString(m_fRho) + _T("\n");
    sRet = sRet + tab + _T("HasT : ") + (m_bHasT ? _T("1\n") : _T("0\n"));
    sRet = sRet + tab + _T("Iterate : ") + appIntToString(static_cast<INT>(m_uiIterate)) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================