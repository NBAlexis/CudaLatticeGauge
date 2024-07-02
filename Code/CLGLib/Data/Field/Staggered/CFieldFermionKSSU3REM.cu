//=============================================================================
// FILENAME : CFieldFermionKSSU3R.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/23/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3REM)

#pragma region DOperator

#pragma region kernel

#if discard

/**
 * Explain:
 * qBz: u_y = exp(i qBz x), u_x(L_x) = exp(-i Lx y)
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSREM_Simple(
        const deviceSU3Vector* __restrict__ pDeviceData,
        const deviceSU3* __restrict__ pGauge,
        const SIndex* __restrict__ pGaugeMove,
        const SIndex* __restrict__ pFermionMove,
        const BYTE* __restrict__ pEtaTable,
        deviceSU3Vector* pResultData,
        Real fam,
        Real fqBz,
        BYTE byGaugeType,
        UBOOL bTwistedBoundary,
        SSmallInt4 sCenter,
        BYTE byFieldId,
        UBOOL bDDagger,
        EOperatorCoefficientType eCoeff,
        Real fCoeff,
        CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];

        if (0 == idir || 1 == idir)
        {
            const SSmallInt4 x_m_site = __deviceSiteIndexToInt4(x_m_mu_Gauge.m_uiSiteIndex);
            const Real forwardPhase = __deviceGetMagneticPhase(sSite4, sCenter, static_cast<BYTE>(idir), byGaugeType, bTwistedBoundary) * fqBz;
            const Real backwardPhase = __deviceGetMagneticPhase(x_m_site, sCenter, static_cast<BYTE>(idir), byGaugeType, bTwistedBoundary) * fqBz;

            x_Gauge_element.MulComp(_make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(backwardPhase), _sin(backwardPhase)));
        }

        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0) * eta_mu);
        }
        else
        {
            u_phi_x_p_m.MulReal(eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        //u_dagger_phi_x_m_m.MulReal(bShiftCenter ? eta_mu2 : eta_mu);
        u_dagger_phi_x_m_m.MulReal(eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        //u_phi_x_p_m.MulReal(eta_mu);
        result.Add(u_phi_x_p_m);
    }

    pResultData[uiSiteIndex].MulReal(fam);
    if (bDDagger)
    {
        pResultData[uiSiteIndex].Sub(result);
    }
    else
    {
        pResultData[uiSiteIndex].Add(result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].MulReal(fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex].MulComp(cCoeff);
        break;
    }
}
#endif

/**
* When link n and n+mu, the coordinate is stick with n
* When link n and n-mu, the coordinate is stick with n-mu
* Irrelavent with tau
* Optimization: bXorY removed, block.x *= 2 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PREM_XYTerm(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const BYTE * __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    Real fCharge,
    SSmallInt4 sCenter,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu  = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
                         + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {            
            this_eta_tau = this_eta_tau + 1;
        }

        deviceSU3Vector right = _deviceVXXTauOptimizedEM(pGauge, pU1, sSite4, fCharge,
            byGaugeFieldId, bXorY, bPlusMu, bPlusTau).MulVector(
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        right.MulReal(sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + F(0.5));

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            result.Sub(right);
        }
        else
        {
            result.Add(right);
        }
    }

    if (bDDagger)
    {
        result.MulReal(F(-0.25) * fOmega);
    }
    else
    {
        result.MulReal(F(0.25) * fOmega);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pResultData[uiSiteIndex].Add(result);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PREM_XYTau_Term(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    Real fCharge,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        
        const deviceSU3Vector right = _deviceVXYTOptimizedEM(
            pGauge, pU1, sSite4, fCharge,
            byGaugeFieldId, bPlusX, bPlusY, bPlusT)
        .MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? (sSite4.y + sSite4.z) : (site_target.y + site_target.z + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            result.Sub(right);
        }
        else
        {
            result.Add(right);
        }
    }

    if (bDDagger)
    {
        result.MulReal(-F(0.125) * fOmega);
    }
    else
    {
        result.MulReal(F(0.125) * fOmega);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pResultData[uiSiteIndex].Add(result);
}

#pragma endregion

#pragma region Derivate

#if discard

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceREM_Simple(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    Real fqBz,
    BYTE byGaugeType,
    UBOOL bTwistedBoundary,
    SSmallInt4 sCenter,
    BYTE byFieldId)
{
    intokernalInt4;

    const Real fPhaseX = __deviceGetMagneticPhase(sSite4, sCenter, 0, byGaugeType, bTwistedBoundary) * fqBz;
    const Real fPhaseY = __deviceGetMagneticPhase(sSite4, sCenter, 1, byGaugeType, bTwistedBoundary) * fqBz;

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        deviceSU3 gaugeElement = pGauge[linkIndex];
        if (0 == idir)
        {
            gaugeElement.MulComp(_make_cuComplex(_cos(fPhaseX), _sin(fPhaseX)));
        }
        else if (1 == idir)
        {
            gaugeElement.MulComp(_make_cuComplex(_cos(fPhaseY), _sin(fPhaseY)));
        }

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[uiR];
            const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];

            deviceSU3Vector toContract = gaugeElement.MulVector(phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = gaugeElement.MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            thisTerm.Add(deviceSU3::makeSU3ContractV(toContract, phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR]);
            }

            thisTerm.Ta();

            pForce[linkIndex].Sub(thisTerm);
        }
    }

}

#endif

/**
 * Have n, n->n1, n->n2,
 * 1. we need to obtain V_(n, n1) , V_(n, n2)
 * 2. we need phi(n1), phi(n2), phid(n1), phid(n2)
 *
 * byContribution: 0 for mu, 1 for tau, 2 for both mu and tau
 *
 * iTau = 1 for +t, -1 for -t
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PREM_XYTerm( 
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    deviceSU3* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    Real fCharge,
    BYTE byMu, INT iTau,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    //const SSmallInt4 middleSite = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    //From now on, site_n1 is smiddle
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    
    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + F(0.5))
        : static_cast<Real>(_DC_Centerx - site_n1.x - F(0.5));

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLinkEM(pGauge, pU1, fCharge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLinkEM(pGauge, pU1, fCharge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];
        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
        deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            phi1.MulReal(F(-1.0));
            phi3.MulReal(F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            phi2.MulReal(F(-1.0));
            phi4.MulReal(F(-1.0));
        }
        deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        const Real eta_tau = (iTau > 0 ? 
            ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1) 
            : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1) )
            ? F(-1.0) : F(1.0);
        res.MulReal(OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);

        //For mu
        if (0 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
            pForce[linkIndex].Sub(res);
        }

        //For tau
        if (1 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            if (iTau > 0)
            {
                pForce[linkIndex].Sub(res);
            }
            else
            {
                pForce[linkIndex].Add(res);
            }
        }
    }
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PREM_XYTau_Term(
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
#if _CLG_DOUBLEFLOAT
    Real fOmega,
#else
    DOUBLE fOmega,
#endif
    Real fCharge,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

    //Why use sn2? shouldn't it be sn1?
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLinkEM(pGauge, pU1, fCharge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLinkEM(pGauge, pU1, fCharge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
        deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            phi1.MulReal(F(-1.0));
            phi3.MulReal(F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            phi2.MulReal(F(-1.0));
            phi4.MulReal(F(-1.0));
        }
        deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
        //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
        //The sign of the second term is because of 'dagger'
        //However, eta124(n1) = -eta124(n2), so use Add directly.
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        res.MulReal(OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);

        //Use eta124 of n2 so Add left Sub right
        //Change to use eta124 of n1, Sub left and Add right
        if (pathLdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);
            pForce[linkIndex].Add(res);
        }

        if (pathRdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);
            pForce[linkIndex].Sub(res);
        }
    }

}

/*
__global__ void _CLG_LAUNCH_BOUND
_giveupkernelDFermionKSForce_PR_XYTau_Term2(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    INT pathLdir1, INT pathLdir2, INT pathLdir3)
{
    intokernalInt4;
    const INT full[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Ldirs[3];
    INT Rdirs[3];
    BYTE Llength, Rlength;

    #pragma unroll
    for (BYTE sep = 0; sep < 4; ++sep)
    {
        _deviceSeperate(full, sep, 3, Ldirs, Rdirs, Llength, Rlength);

        const UBOOL bHasLeft = Llength > 0 && Ldirs[0] > 0;
        const UBOOL bHasRight = Rlength > 0 && Rdirs[0] > 0;
        if (!bHasLeft && !bHasRight)
        {
            continue;
        }

        //=================================
        // 1. Find n1, n2
        const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
        const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
        const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
        const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

        //Why use sn2? shouldn't it be sn1?
        const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex));
        //=================================
        // 2. Find V(n,n1), V(n,n2)
        const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
        const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
            const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
            deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
            deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
            deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
            if (sn1.NeedToOpposite())
            {
                phi1.MulReal(F(-1.0));
                phi3.MulReal(F(-1.0));
            }
            if (sn2.NeedToOpposite())
            {
                phi2.MulReal(F(-1.0));
                phi4.MulReal(F(-1.0));
            }
            deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
            //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
            //The sign of the second term is because of 'dagger'
            //However, eta124(n1) = -eta124(n2), so use Add directly.
            res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
            res.Ta();
            res.MulReal(OneOver48 * fOmega * pNumerators[rfieldId] * eta124);

            if (bHasLeft)
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, Ldirs[0] - 1);
                pForce[linkIndex].Add(res);
            }

            if (bHasRight)
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, Rdirs[0] - 1);
                pForce[linkIndex].Sub(res);
            }
        }
    }
}

*/


#if discard
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaREM(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byDir,
    Real fQBz,
    BYTE byGaugeType,
    UBOOL bTwistedBoundary,
    SSmallInt4 sCenter)
{
    intokernalInt4;

    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    BYTE eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1));
    BYTE eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir) & 1));

    deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];

    //============= add phase ===============
    //Note that, the coordinate should be used as same as the gauge
    if (0 == byDir || 1 == byDir)
    {
        const SSmallInt4 x_m_site = __deviceSiteIndexToInt4(x_m_mu_Gauge.m_uiSiteIndex);
        const Real forwardPhase = __deviceGetMagneticPhase(sSite4, sCenter, byDir, byGaugeType, bTwistedBoundary) * fQBz;
        const Real backwardPhase = __deviceGetMagneticPhase(x_m_site, sCenter, byDir, byGaugeType, bTwistedBoundary) * fQBz;
        x_Gauge_element.MulComp(_make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
        x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(backwardPhase), _sin(backwardPhase)));
    }

    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    pMe[uiSiteIndex] = x_Gauge_element.MulVector(pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        pMe[uiSiteIndex].MulReal(F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        pMe[uiSiteIndex].Sub(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        pMe[uiSiteIndex].Add(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    pMe[uiSiteIndex].MulReal(F(0.5));
}
#endif

#pragma endregion


#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3REM::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byEMFieldID));
    if (NULL == pU1)
    {
        appCrucial(_T("Background field cannot be found!\n"));
        return;
    }

    CFieldFermionKSSU3GammaEM::DOperatorEM(pTargetBuffer, pBuffer, pGaugeBuffer, pU1->m_pDeviceData, f2am, m_fQ,
        m_bEachSiteEta, bDagger, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);

    preparethread;

#if 1
    _kernelDFermionKS_PREM_XYTerm << <block, threads >> > (
        pSource,
        pGauge,
        pU1->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        byGaugeFieldId,
        CCommonData::m_fOmega,
        m_fQ,
        _HC_Center,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

#endif

#if 1

    _kernelDFermionKS_PREM_XYTau_Term << <block, threads >> > (
        pSource,
        pGauge,
        pU1->m_pDeviceData,
        pTarget,
        m_byFieldId,
        byGaugeFieldId,
        CCommonData::m_fOmega,
        m_fQ,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

#endif
}

void CFieldFermionKSSU3REM::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byEMFieldID));
    if (NULL == pU1)
    {
        appCrucial(_T("Background field cannot be found!\n"));
        return;
    }

    CFieldFermionKSSU3GammaEM::KSForceEM(pForce, pGaugeBuffer, pU1->m_pDeviceData, m_fQ,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        m_byFieldId);

    preparethread;

#if 1

    #pragma region X Y Term

    INT mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        INT dirs[6][3] =
        {
            {4, mu[imu] + 1, mu[imu] + 1},
            {mu[imu] + 1, 4, mu[imu] + 1},
            {mu[imu] + 1, mu[imu] + 1, 4},
            //{4, -mu[imu] - 1, -mu[imu] - 1},
            //{-mu[imu] - 1, 4, -mu[imu] - 1},
            //{-mu[imu] - 1, -mu[imu] - 1, 4},
            {mu[imu] + 1, mu[imu] + 1, -4},
            {mu[imu] + 1, -4, mu[imu] + 1},
            {-4, mu[imu] + 1, mu[imu] + 1},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };
        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 3},
            {0, 1, 0, 3},
            {0, 0, 1, 3},
            //{1, 3, 0, 0},
            //{3, 2, 3, 0},
            //{3, 0, 2, 3},
            {0, 0, 3, 1},
            {0, 3, 2, 3},
            {3, 2, 0, 3},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (3 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _kernelDFermionKSForce_PREM_XYTerm << <block, threads >> > (
                    (const deviceSU3*)pGaugeBuffer,
                    pU1->m_pDeviceData,
                    (deviceSU3*)pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    m_pRationalFieldPointers,
                    m_pMDNumerator,
                    m_rMD.m_uiDegree,
                    m_byFieldId,
                    CCommonData::m_fOmega,
                    m_fQ,
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

    #pragma endregion
#endif

#if 1

    #pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    //INT linkTypes[4][3] =
    //{
    //    {1, 2, 4},
    //    {1, 2, -4},
    //    {-1, 2, 4},
    //    {-1, 2, -4}
    //};
    INT linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        INT sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            //bearly no change of time, because force calculation is not frequent
            /*
            _giveupkernelDFermionKSForce_PR_XYTau_Term2 << <block, threads >> > (
                (const deviceSU3*)pGaugeBuffer,
                (deviceSU3*)pForce,
                m_pRationalFieldPointers,
                m_pMDNumerator,
                m_rMD.m_uiDegree,
                m_byFieldId,
                CCommonData::m_fOmega,
                sixlinks[isixtype][0], sixlinks[isixtype][1], sixlinks[isixtype][2]
                );
            */

            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _kernelDFermionKSForce_PREM_XYTau_Term << <block, threads >> > (
                        (const deviceSU3*)pGaugeBuffer,
                        pU1->m_pDeviceData,
                        (deviceSU3*)pForce,
                        m_pRationalFieldPointers,
                        m_pMDNumerator,
                        m_rMD.m_uiDegree,
                        m_byFieldId,
                        CCommonData::m_fOmega,
                        m_fQ,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }
    
    #pragma endregion
#endif
}

#pragma endregion

CFieldFermionKSSU3REM::CFieldFermionKSSU3REM()
    : CFieldFermionKSSU3()
    , m_byEMFieldID(0)
    , m_fQ(F(0.0))
    , m_pDevicePathBuffer(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePathBuffer, sizeof(INT) * 4));
}

CFieldFermionKSSU3REM::~CFieldFermionKSSU3REM()
{
    checkCudaErrors(cudaFree(m_pDevicePathBuffer));
}

void CFieldFermionKSSU3REM::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);
    m_bEachSiteEta = TRUE;

    Real fValue = F(0.0);
    if (params.FetchValueReal(_T("Charge"), fValue))
    {
        m_fQ = fValue;
    }
    INT iValue = 0;
    if (params.FetchValueINT(_T("EMFieldID"), iValue))
    {
        m_byEMFieldID = static_cast<BYTE>(iValue);
    }
}

void CFieldFermionKSSU3REM::ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
{
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byEMFieldID));
    if (NULL == pU1)
    {
        appCrucial(_T("Background field cannot be found!\n"));
        return;
    }
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    InitialField(EFIT_Zero);

    //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
    UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
               || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);

    //appApplyGammaKSEM is used for simulation, where (2a) is multiplied, so coefficient should be 0.5
    //note that, finally the 2a is multiplied again, we decide to leave it to the Mathematica

    CFieldFermionKSSU3GammaEM::appApplyGammaKSEM(
        m_pDeviceData,
        pPooled->m_pDeviceData, 
        pFieldSU3->m_pDeviceData, 
        pU1->m_pDeviceData, 
        m_fQ,
        eGamma, 
        m_bEachSiteEta, 
        FALSE, 
        F(0.5), 
        bImag ? EOCT_Complex : EOCT_None,
        F(1.0), 
        bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,
        m_byFieldId, 
        1);

    pPooled->Return();
}

void CFieldFermionKSSU3REM::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);

    CFieldFermionKSSU3REM* pTarget = dynamic_cast<CFieldFermionKSSU3REM*>(U);
    pTarget->m_byEMFieldID = m_byEMFieldID;
    pTarget->m_fQ = m_fQ;
}

CCString CFieldFermionKSSU3REM::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKSSU3::GetInfos(tab);
    sRet = sRet + tab + _T("Omega : ") + appToString(CCommonData::m_fOmega) + _T("\n");
    sRet = sRet + tab + _T("Background Field ID : ") + appToString(m_byEMFieldID) + _T("\n");
    sRet = sRet + tab + _T("fQ : ") + appToString(m_fQ) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================