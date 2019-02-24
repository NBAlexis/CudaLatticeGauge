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

__global__ void
_CLG_LAUNCH_BOUND
_kernelAPEStoutSU3(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fRhoS,
    Real fRhoT,
    BYTE maxDir
)
{
    intokernal;

    for (BYTE dir = 0; dir < maxDir; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulDagger(gauge[uiLinkIndex]);
        omega.Ta();
        omega = 0 == _DC_ExpPrecision ? 
            omega.QuickExp(3 == dir ? fRhoT : fRhoS)
          : omega.ExpReal(3 == dir ? fRhoT : fRhoS, _DC_ExpPrecision);
        omega.Mul(gauge[uiLinkIndex]);
        gauge[uiLinkIndex] = omega;
    }
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeSmearingAPEStout)


void CGaugeSmearingAPEStout::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    if (!params.FetchValueReal(_T("RhoS"), m_fRhoS))
    {
        appCrucial(_T("CGaugeSmearingAPEStout: rhos not set, set to 0.1 by defualt."));
    }

    if (!params.FetchValueReal(_T("RhoT"), m_fRhoT))
    {
        appCrucial(_T("CGaugeSmearingAPEStout: rhot not set, set to 0.0 by defualt."));
    }
    m_bHasT = appAbs(m_fRhoT) > F(0.00001);
    INT iValue;
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

    appGeneral(_T(" ---- CGaugeSmearingAPEStout created with rhoS = %f, rhoT = %f, iterate = %d\n"), m_fRhoS, m_fRhoT, m_uiIterate);
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

    preparethread;
    for (UINT i = 0; i < m_uiIterate; ++i)
    {
        _kernelAPEStoutSU3 << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            pStapleSU3->m_pDeviceData,
            m_fRhoS,
            m_fRhoT,
            m_bHasT ? 4 : 3
            );

        pGaugeSU3->CalculateOnlyStaple(pStapleSU3);
    }
}

CCString CGaugeSmearingAPEStout::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeSmearingAPEStout\n");
    sRet = sRet + tab + _T("RhosRhot : [") + appFloatToString(m_fRhoS) + _T(", ") + appFloatToString(m_fRhoT) +_T("]\n");
    sRet = sRet + tab + _T("Iterate : ") + appIntToString(static_cast<INT>(m_uiIterate)) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================