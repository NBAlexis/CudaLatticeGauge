//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingU1_3D.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingU1_3D)



#pragma region kernels

#pragma region Clover

__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermU1_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    DOUBLE fXSq = (sSite4.x - _DC_Centerx + 0.5);
    fXSq = fXSq * fXSq;
    DOUBLE fYSq = (sSite4.y - _DC_Centery + 0.5);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const DOUBLE fU13 = fXSq * _device4PlaqutteTermU1(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU23 = fYSq * _device4PlaqutteTermU1(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

#else
    Real fXSq = (sSite4.x - _DC_Centerx + F(0.5));
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - _DC_Centery + F(0.5));
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const Real fU13 = fXSq * _device4PlaqutteTermU1(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const Real fU23 = fYSq * _device4PlaqutteTermU1(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

#endif

    results[uiSiteIndex] = (fU23 + fU13) * betaOverN * fOmegaSq;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermU1_XYZ_Shifted3D(
    BYTE byFieldId,
    UBOOL bTorus,
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq
#else
    Real betaOverN, Real fOmegaSq
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    //deviceSU3 plaqSum = deviceSU3::makeSU3Zero();
    BYTE idx[4] = { 2, 1, 2, 1 };
    BYTE byOtherDir[4] = { 2, 2, 0, 1 };

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //mu = idir, nu = 4, i = mu
        CLGComplex stap = _deviceStapleTermGfactorU1(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[idir],
            idx[idir],
            TRUE);
        if (2 == idir)
        {
            stap = _cuCaddf(stap, _deviceStapleTermGfactorU1(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
                idir,
                byOtherDir[idir + 1],
                idx[idir + 1],
                TRUE));
        }
        
        CLGComplex force = pDeviceData[linkIndex];
        force = _cuCmulf(force, _cuConjf(stap));
        pForceData[linkIndex].y = pForceData[linkIndex].y + force.y * betaOverN;
    }
}

#pragma endregion

#pragma region Chair energy

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermU1_Term1234_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
#else
    Real betaOverN, Real fOmega,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = 0.125 * betaOverN;
    const DOUBLE fXOmega = (sSite4.x - _DC_Centerx + 0.5) * fOmega;
    const DOUBLE fYOmega = (sSite4.y - _DC_Centery + 0.5) * fOmega;

    //===============
    //- x Omega V312
    const DOUBLE fV312 = -fXOmega * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);

    //===============
    //+ y Omega V321
    const DOUBLE fV321 = fYOmega * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

#else
    betaOverN = F(0.125) * betaOverN;
    const Real fXOmega = (sSite4.x - _DC_Centerx + F(0.5)) * fOmega;
    const Real fYOmega = (sSite4.y - _DC_Centery + F(0.5)) * fOmega;

    //===============
    //-x Omega V312
    const Real fV312 = -fXOmega * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);

    //===============
    //+y Omega V321
    const Real fV321 = fYOmega * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);
#endif

    results[uiSiteIndex] = (fV312 + fV321) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermU1_Term5_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = 0.125 * betaOverN;
    const DOUBLE fXYOmega2 = (sSite4.x - _DC_Centerx + 0.5) * (sSite4.y - _DC_Centery + 0.5) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#else
    betaOverN = F(0.125) * betaOverN;
    const Real fXYOmega2 = (sSite4.x - _DC_Centerx + F(0.5)) * (sSite4.y - _DC_Centery + F(0.5)) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    const Real fV132 = fXYOmega2 * _deviceChairTermU1(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#endif

    results[uiSiteIndex] = fV132 * betaOverN;
}

#pragma endregion

#pragma region Chair force

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermU1_Term1_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V312
    //add force for dir=3
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const CLGComplex staple_term1_4 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);
    CLGComplex force4 = pDeviceData[uiLink4];
    force4 = _cuCmulf(force4, _cuConjf(staple_term1_4));
    pForceData[uiLink4].y = pForceData[uiLink4].y + force4.y * betaOverN;
    //}


    //===============
    //+x Omega V412
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const CLGComplex staple_term1_2 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 2, _deviceHiShifted0);
    CLGComplex force2 = pDeviceData[uiLink2];
    force2 = _cuCmulf(force2, _cuConjf(staple_term1_2));
    pForceData[uiLink2].y = pForceData[uiLink2].y + force2.y * betaOverN;
   // }

    //===============
    //+x Omega V412
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const CLGComplex staple_term1_1 = _deviceStapleChairTerm2ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);
    CLGComplex force1 = pDeviceData[uiLink1];
    force1 = _cuCmulf(force1, _cuConjf(staple_term1_1));
    pForceData[uiLink1].y = pForceData[uiLink1].y + force1.y * betaOverN;
    //}
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermU1_Term3_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+ y Omega V321
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const CLGComplex staple_term3_4 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);
    CLGComplex force4 = pDeviceData[uiLink4];
    force4 = _cuCmulf(force4, _cuConjf(staple_term3_4));
    pForceData[uiLink4].y = pForceData[uiLink4].y + force4.y * betaOverN;
    //}

    //===============
    //+ y Omega V421
    //add force for mu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const CLGComplex staple_term3_1 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 2, _deviceHiShifted1);
    CLGComplex force1 = pDeviceData[uiLink1];
    force1 = _cuCmulf(force1, _cuConjf(staple_term3_1));
    pForceData[uiLink1].y = pForceData[uiLink1].y + force1.y * betaOverN;
    //}


    //===============
    //+ y Omega V421
    //add force for mu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const CLGComplex staple_term3_2 = _deviceStapleChairTerm2ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);
    CLGComplex force2 = pDeviceData[uiLink2];
    force2 = _cuCmulf(force2, _cuConjf(staple_term3_2));
    pForceData[uiLink2].y = pForceData[uiLink2].y + force2.y * betaOverN;
    //}

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermU1_Term5_Shifted3D(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq
#else
    Real betaOverN, Real fOmegaSq
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //- Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const CLGComplex staple_term5_1 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);
    CLGComplex force1 = pDeviceData[uiLink1];
    force1 = _cuCmulf(force1, _cuConjf(staple_term5_1));
    pForceData[uiLink1].y = pForceData[uiLink1].y + force1.y * betaOverN;
    //}

    //===============
    //- Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const CLGComplex staple_term5_2 = _deviceStapleChairTerm1ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 0, _deviceHiShifted2);
    CLGComplex force2 = pDeviceData[uiLink2];
    force2 = _cuCmulf(force2, _cuConjf(staple_term5_2));
    pForceData[uiLink2].y = pForceData[uiLink2].y + force2.y * betaOverN;
    //}

    //===============
    //- Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const CLGComplex staple_term5_3 = _deviceStapleChairTerm2ShiftedU1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);
    CLGComplex force3 = pDeviceData[uiLink3];
    force3 = _cuCmulf(force3, _cuConjf(staple_term5_3));
    pForceData[uiLink3].y = pForceData[uiLink3].y + force3.y * betaOverN;
    //}
}

#pragma endregion

#pragma endregion

CActionGaugePlaquetteRotatingU1_3D::CActionGaugePlaquetteRotatingU1_3D()
    : CActionGaugePlaquetteRotatingU1()
{
}

UBOOL CActionGaugePlaquetteRotatingU1_3D::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
#if !_CLG_DOUBLEFLOAT
    pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
#else
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);
#endif

    const CFieldGaugeU1* pGaugeU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldGaugeU1* pForceU1 = dynamic_cast<CFieldGaugeU1*>(pForce);
    if (NULL == pGaugeU1 || NULL == pForceU1)
    {
        appCrucial(_T("CActionGaugePlaquetteRotatingU1 only work with U1 now.\n"));
        return TRUE;
    }

    preparethread;


    if (!m_bShiftHalfCoord)
    {
        appCrucial(_T("Dirichlet not supported yet"));
    }
    else
    {

        _kernelAddForce4PlaqutteTermU1_XYZ_Shifted3D << <block, threads >> > (pGaugeU1->m_byFieldId, FALSE, pGaugeU1->m_pDeviceData,
            pForceU1->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        _kernelAddForceChairTermU1_Term1_Shifted3D << <block, threads >> > (pGaugeU1->m_byFieldId, pGaugeU1->m_pDeviceData,
            pForceU1->m_pDeviceData, m_fBetaOverN, m_fOmega);
        
        _kernelAddForceChairTermU1_Term3_Shifted3D << <block, threads >> > (pGaugeU1->m_byFieldId, pGaugeU1->m_pDeviceData,
            pForceU1->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermU1_Term5_Shifted3D << <block, threads >> > (pGaugeU1->m_byFieldId, pGaugeU1->m_pDeviceData,
            pForceU1->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquetteRotatingU1_3D::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    if (m_bCloverEnergy)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
    }
    else
    {
        //m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
        appCrucial(_T("Dirichlet not supported yet"));
    }
    
    const CFieldGaugeU1* pGaugeU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    if (NULL == pGaugeU1)
    {
        appCrucial(_T("CActionGaugePlaquetteRotatingU1 only work with U1 now.\n"));
        return m_fNewEnergy;
    }

    preparethread;

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    if (m_bShiftHalfCoord)
    {

        _kernelAdd4PlaqutteTermU1_Shifted3D << <block, threads >> > (
            pGaugeU1->m_byFieldId,
            pGaugeU1->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
        
        _kernelAddChairTermU1_Term1234_Shifted3D << <block, threads >> > (
            pGaugeU1->m_byFieldId,
            pGaugeU1->m_pDeviceData,
            m_fBetaOverN, 
            m_fOmega, 
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermU1_Term5_Shifted3D << <block, threads >> > (
            pGaugeU1->m_byFieldId,
            pGaugeU1->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }
    else
    {
        appCrucial(_T("Dirichlet not supported yet"));
    }


    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteRotatingU1_3D::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================