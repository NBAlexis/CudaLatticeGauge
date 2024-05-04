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
_kernelAPEStoutSU3(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    BYTE maxlink,
    Real fRho)
{
    intokernal;

    for (BYTE dir = 0; dir < maxlink; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulDagger(gauge[uiLinkIndex]);
        //Note, q = -i (staple.u).TA(), and exp(i q) = exp((staple.u).TA())
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
    CGaugeSmearing::Initial(pOwner, params);

    if (!params.FetchValueReal(_T("Rho"), m_fRho))
    {
        appCrucial(_T("CGaugeSmearingAPEStout: rhos not set, set to 0.1 by defualt."));
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
    CFieldGaugeSU3* pStapleSU3 = NULL;
    if (NULL == pStaple)
    {
        if (NULL == m_pTmpStaple)
        {
            m_pTmpStaple = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        }
        pStapleSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pTmpStaple);
        if (m_bHasT)
        {
            pGauge->CalculateOnlyStaple(pStapleSU3);
        }
    }
    else if (EFT_GaugeSU3 != pStaple->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    else
    {
        pStapleSU3 = dynamic_cast<CFieldGaugeSU3*>(pStaple);
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);
    appParanoiac(_T("GaugeSmearing : CGaugeSmearingAPEStout\n"));
    const BYTE maxlink = static_cast<BYTE>(m_bHasT ? _HC_Dim : _HC_Dim - 1);
    preparethread;
    for (UINT i = 0; i < m_uiIterate; ++i)
    {
        if (!m_bHasT)
        {
            CalculateSpatialFatLink(pGauge, pStapleSU3);
        }

        _kernelAPEStoutSU3 << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            pStapleSU3->m_pDeviceData,
            maxlink,
            m_fRho);

        if (m_bHasT)
        {
            pGaugeSU3->CalculateOnlyStaple(pStapleSU3);
        }
    }

    if (NULL != pStaple && !m_bHasT)
    {
        pGaugeSU3->CalculateOnlyStaple(pStaple);
    }
}

CCString CGaugeSmearingAPEStout::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeSmearingAPEStout\n");
    sRet = sRet + tab + _T("Rhos : ") + appToString(m_fRho) + _T("\n");
    sRet = sRet + tab + _T("HasT : ") + (m_bHasT ? _T("1\n") : _T("0\n"));
    sRet = sRet + tab + _T("Iterate : ") + appToString(static_cast<INT>(m_uiIterate)) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================