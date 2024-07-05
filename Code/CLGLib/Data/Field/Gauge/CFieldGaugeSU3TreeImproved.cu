//=============================================================================
// FILENAME : CFieldGaugeSU3TreeImproved.cu
// 
// DESCRIPTION:
// On plaq term and rect term
//
//
// REVISION:
//  [10/03/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldGaugeSU3TreeImproved.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_Normal(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE betaOverNTimesCRect,
    DOUBLE* results
#else
    const Real betaOverNTimesCRect,
    Real* results
#endif
)
{
    intokernalInt4;

#if !_CLG_DOUBLEFLOAT
    DOUBLE res = 0.0;
#else
    Real res = F(0.0);
#endif
    for (BYTE byMu = 0; byMu < _DC_Dir; ++byMu)
    {
        for (BYTE byNu = byMu + 1; byNu < _DC_Dir; ++byNu)
        {
            INT path[6] = { __fwd(byMu), __fwd(byMu), __fwd(byNu), __bck(byMu), __bck(byMu), __bck(byNu) };
            res += _deviceLink(pDeviceData, sSite4, 6, byFieldId, path).ReTr();

            path[0] = __fwd(byMu);
            path[1] = __fwd(byNu);
            path[2] = __fwd(byNu);
            path[3] = __bck(byMu);
            path[4] = __bck(byNu);
            path[5] = __bck(byNu);
            res += _deviceLink(pDeviceData, sSite4, 6, byFieldId, path).ReTr();
        }
    }

    results[uiSiteIndex] = (F(36.0) - res) * betaOverNTimesCRect;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_Clover(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE betaOverNTimesCRect,
    DOUBLE* results
#else
    const Real betaOverNTimesCRect,
    Real* results
#endif
)
{
    intokernalInt4;

#if !_CLG_DOUBLEFLOAT
    DOUBLE res = 0.0;
#else
    Real res = F(0.0);
#endif
    for (BYTE byMu = 0; byMu < _DC_Dir; ++byMu)
    {
        for (BYTE byNu = byMu + 1; byNu < _DC_Dir; ++byNu)
        {
            res += F(12.0) - _deviceCloverRectangleRetr(byFieldId, pDeviceData, sSite4, byMu, byNu);
        }
    }

    results[uiSiteIndex] = res * betaOverNTimesCRect * F(0.5);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_Force(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    Real betaOverNTimesCRect,
    deviceSU3* pForce)
{
    intokernalInt4;
    betaOverNTimesCRect = betaOverNTimesCRect * F(-0.5);

    for (BYTE byMu = 0; byMu < _DC_Dir; ++byMu)
    {
        deviceSU3 res = deviceSU3::makeSU3Zero();
        for (BYTE byNu = 0; byNu < _DC_Dir; ++byNu)
        {
            if (byNu == byMu)
            {
                continue;
            }

            //Calculate staples
            // ---------
            // |       |
            // x--------
            INT path[5] = { __fwd(byNu), __fwd(byMu), __fwd(byMu), __bck(byNu), __bck(byMu) };
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));

            path[0] = __bck(byNu);
            path[1] = __fwd(byMu);
            path[2] = __fwd(byMu);
            path[3] = __fwd(byNu);
            path[4] = __bck(byMu);
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));

            // ---------
            // |       |
            // ----x----
            path[0] = __bck(byMu);
            path[1] = __fwd(byNu);
            path[2] = __fwd(byMu);
            path[3] = __fwd(byMu);
            path[4] = __bck(byNu);
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));

            path[0] = __bck(byMu);
            path[1] = __bck(byNu);
            path[2] = __fwd(byMu);
            path[3] = __fwd(byMu);
            path[4] = __fwd(byNu);
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));

            // ----
            // |  |
            // |  |
            // x---
            path[0] = __fwd(byNu);
            path[1] = __fwd(byNu);
            path[2] = __fwd(byMu);
            path[3] = __bck(byNu);
            path[4] = __bck(byNu);
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));

            path[0] = __bck(byNu);
            path[1] = __bck(byNu);
            path[2] = __fwd(byMu);
            path[3] = __fwd(byNu);
            path[4] = __fwd(byNu);
            res.Add(_deviceLink(pDeviceData, sSite4, 5, byFieldId, path));
        }
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
        deviceSU3 force(pDeviceData[uiLinkIndex]);
        force.MulDagger(res);
        force.Ta();
        force.MulReal(betaOverNTimesCRect);

        //force is additive
        pForce[uiLinkIndex].Add(force);
    }
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3TreeImproved)

CFieldGaugeSU3TreeImproved::CFieldGaugeSU3TreeImproved()
    : CFieldGaugeSU3()
    , m_fRectOverPlaq(F(-0.05))
{
    
}

CFieldGaugeSU3TreeImproved::~CFieldGaugeSU3TreeImproved()
{
    
}

void CFieldGaugeSU3TreeImproved::InitialOtherParameters(CParameters& param)
{
    CFieldGaugeSU3::InitialOtherParameters(param);

    Real fValue = F(-0.05);
    if (param.FetchValueReal(_T("RectOverPlaq"), fValue))
    {
        m_fRectOverPlaq = fValue;
    }
}

void CFieldGaugeSU3TreeImproved::CopyTo(CField* U) const
{
    CFieldGaugeSU3::CopyTo(U);
    CFieldGaugeSU3TreeImproved* pFieldGauge = dynamic_cast<CFieldGaugeSU3TreeImproved*>(U);
    if (NULL == pFieldGauge)
    {
        return;
    }
    pFieldGauge->m_fRectOverPlaq = m_fRectOverPlaq;
}

void CFieldGaugeSU3TreeImproved::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const
{
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    CFieldGaugeSU3::CalculateForceAndStaple(pForce, pStaple, betaOverN);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    preparethread;
    _kernelPlaqutteEnergy_Rectangular_Force << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        m_fRectOverPlaq * betaOverN,
        pForceSU3->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
#else
Real CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergy(Real betaOverN) const
#endif
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergy(betaOverN);
#else
    Real fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergy(betaOverN);
#endif
    preparethread;
    _kernelPlaqutteEnergy_Rectangular_Normal << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        m_fRectOverPlaq * betaOverN,
        _D_RealThreadBuffer);

    fPlaqTerm += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return fPlaqTerm;
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
#else
Real CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergyUseClover(Real betaOverN) const
#endif
{
#if !_CLG_DOUBLEFLOAT
    DOUBLE fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergyUseClover(betaOverN);
#else
    Real fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergyUseClover(betaOverN);
#endif

    preparethread;
    _kernelPlaqutteEnergy_Rectangular_Clover << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        m_fRectOverPlaq * betaOverN,
        _D_RealThreadBuffer);

    fPlaqTerm += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    return fPlaqTerm;
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const
#else
Real CFieldGaugeSU3TreeImproved::CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge* pStaple) const
#endif
{
    appGeneral(_T("Calculate energy using stable is not supported...\n"));
    //We don't known whether to use clover or not, so use normal

#if !_CLG_DOUBLEFLOAT
    DOUBLE fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergyUsingStable(betaOverN, pStaple);
#else
    Real fPlaqTerm = CFieldGaugeSU3::CalculatePlaqutteEnergyUsingStable(betaOverN, pStaple);
#endif
    preparethread;
    _kernelPlaqutteEnergy_Rectangular_Normal << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        m_fRectOverPlaq * betaOverN,
        _D_RealThreadBuffer);

    fPlaqTerm += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return fPlaqTerm;
}

CCString CFieldGaugeSU3TreeImproved::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldGaugeSU3::GetInfos(tab);
    sRet = sRet + appToString(m_fRectOverPlaq) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================