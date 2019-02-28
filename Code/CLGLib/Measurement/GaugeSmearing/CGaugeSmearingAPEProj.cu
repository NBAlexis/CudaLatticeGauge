//=============================================================================
// FILENAME : CGaugeSmearingAPEProj.cu
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
_kernelAPEProjSU3(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fAlpha,
    BYTE byProjItera
)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulReal(fAlpha * F(0.16666666666666666666666666666667));
        gauge[uiLinkIndex].MulReal(F(1.0) - fAlpha);
        gauge[uiLinkIndex].Add(omega);
        gauge[uiLinkIndex].Proj(byProjItera);
    }
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeSmearingAPEProj)


void CGaugeSmearingAPEProj::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    if (!params.FetchValueReal(_T("Alpha"), m_fAlpha))
    {
        appCrucial(_T("CGaugeSmearingAPEProj: alpha not set, set to 0.1 by defualt."));
    }

    INT iValue;
    if (!params.FetchValueINT(_T("Iterate"), iValue))
    {
        appCrucial(_T("CGaugeSmearingAPEProj: Iterate not set, set to 1 by defualt."));
    }
    else if (iValue > 0)
    {
        m_uiIterate = static_cast<UINT>(iValue);
    }
    else
    {
        appCrucial(_T("CGaugeSmearingAPEProj: Iterate not correct, set to 1 by defualt."));
    }

    if (!params.FetchValueINT(_T("ProjIterate"), iValue))
    {
        appCrucial(_T("CGaugeSmearingAPEProj: ProjIterate not set, set to 6 by defualt."));
    }
    else if (iValue >= 4 && iValue < 255)
    {
        m_byProjIterate = static_cast<BYTE>(iValue);
    }
    else
    {
        appCrucial(_T("CGaugeSmearingAPEProj: ProjIterate not correct, set to 6 by defualt."));
    }

    appGeneral(_T(" ---- CGaugeSmearingAPEProj created with alpha = %f, iterate = %d, ProjIterate = %d\n"), m_fAlpha, m_uiIterate, m_byProjIterate);
}

void CGaugeSmearingAPEProj::GaugeSmearing(CFieldGauge* pGauge, CFieldGauge* pStaple)
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
        _kernelAPEProjSU3 << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            pStapleSU3->m_pDeviceData,
            m_fAlpha,
            m_byProjIterate
            );

        pGaugeSU3->CalculateOnlyStaple(pStapleSU3);
    }
}

CCString CGaugeSmearingAPEProj::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeSmearingAPEProj\n");
    sRet = sRet + tab + _T("alpha : ") + appFloatToString(m_fAlpha) + _T("\n");
    sRet = sRet + tab + _T("Iterate : ") + appIntToString(static_cast<INT>(m_uiIterate)) + _T("\n");
    sRet = sRet + tab + _T("ProjIterate : ") + appIntToString(static_cast<INT>(m_byProjIterate)) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================