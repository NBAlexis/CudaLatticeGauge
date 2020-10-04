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

__global__ void _CLG_LAUNCH_BOUND
_kernelAPEProjSU3(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fAlphaLeft,
    Real fAlphaRight,
    BYTE maxlink,
    BYTE byProjItera)
{
    intokernal;

    for (BYTE dir = 0; dir < maxlink; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulReal(fAlphaRight);
        gauge[uiLinkIndex].MulReal(fAlphaLeft);
        gauge[uiLinkIndex].Add(omega);
        gauge[uiLinkIndex].Proj(byProjItera);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAPEProjSU3CM(
    deviceSU3* gauge,
    const deviceSU3* __restrict__ staple,
    Real fAlphaLeft,
    Real fAlphaRight,
    BYTE maxlink,
    BYTE byIteration)
{
    intokernal;

    for (BYTE dir = 0; dir < maxlink; ++dir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        deviceSU3 omega = staple[uiLinkIndex];
        omega.MulReal(fAlphaRight);
        gauge[uiLinkIndex].MulReal(fAlphaLeft);
        gauge[uiLinkIndex].Add(omega);

        //omega is not used, so use omega as m
        omega = gauge[uiLinkIndex];
        deviceSU3 u = omega;
        u.CabbiboMarinariProj();
        gauge[uiLinkIndex] = u;
        for (BYTE byIt = 1; byIt < byIteration; ++byIt)
        {
            //when it=1, this is u1.m
            //when it=2, this is u2.u1.m ...
            omega = u.MulC(omega);
            u = omega;
            u.CabbiboMarinariProj();

            //when it=1, this is u2.u1
            //when it=2, this is u3.u2.u1
            gauge[uiLinkIndex] = u.MulC(gauge[uiLinkIndex]);
        }
    }
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeSmearingAPEProj)


void CGaugeSmearingAPEProj::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    CGaugeSmearing::Initial(pOwner, params);
    if (!params.FetchValueReal(_T("AlphaLeft"), m_fAlphaLeft))
    {
        appCrucial(_T("CGaugeSmearingAPEProj: AlphaLeft not set, set to 1.0 by defualt."));
    }
    if (!params.FetchValueReal(_T("AlphaRight"), m_fAlphaRight))
    {
        appCrucial(_T("CGaugeSmearingAPEProj: AlphaRight not set, set to 0.1 by defualt."));
    }

    INT iValue;
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

    INT iCMProj = 0;
    params.FetchValueINT(_T("Cabibbo"), iCMProj);
    m_bCMProj = (0 != iCMProj);
    appGeneral(_T(" ---- CGaugeSmearingAPEProj created with u * %f + staple * %f, projectype = %d iterate = %d ProjIterate = %d\n"), m_fAlphaLeft, m_fAlphaRight, m_bCMProj, m_uiIterate, m_byProjIterate);
}

void CGaugeSmearingAPEProj::GaugeSmearing(CFieldGauge* pGauge, CFieldGauge* pStaple)
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
    appParanoiac(_T("GaugeSmearing : CGaugeSmearingAPEProj\n"));
    const BYTE byMaxLink = static_cast<BYTE>(m_bHasT ? _HC_Dim : _HC_Dim - 1);

    preparethread;
    for (UINT i = 0; i < m_uiIterate; ++i)
    {
        if (!m_bHasT)
        {
            CalculateSpatialFatLink(pGauge, pStapleSU3);
        }

        if (m_bCMProj)
        {
            _kernelAPEProjSU3CM << <block, threads >> > (
                pGaugeSU3->m_pDeviceData,
                pStapleSU3->m_pDeviceData,
                m_fAlphaLeft, m_fAlphaRight,
                byMaxLink,
                m_byProjIterate
                );
        }
        else
        {
            _kernelAPEProjSU3 << <block, threads >> > (
                pGaugeSU3->m_pDeviceData,
                pStapleSU3->m_pDeviceData,
                m_fAlphaLeft, m_fAlphaRight, 
                byMaxLink,
                m_byProjIterate
                );
        }

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

CCString CGaugeSmearingAPEProj::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeSmearingAPEProj\n");
    sRet = sRet + tab + _T("# u=u*alphaleft + staple*alpharight\n");
    sRet = sRet + tab + _T("alpha left : ") + appFloatToString(m_fAlphaLeft) + _T("\n");
    sRet = sRet + tab + _T("alpha rihgt : ") + appFloatToString(m_fAlphaRight) + _T("\n");
    sRet = sRet + tab + _T("CabibboMarinari : ") + (m_bCMProj ? _T("1\n") : _T("0\n"));
    sRet = sRet + tab + _T("HasT : ") + (m_bHasT ? _T("1\n") : _T("0\n"));
    sRet = sRet + tab + _T("Iterate : ") + appIntToString(static_cast<INT>(m_uiIterate)) + _T("\n");
    sRet = sRet + tab + _T("ProjIterate : ") + appIntToString(static_cast<INT>(m_byProjIterate)) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================