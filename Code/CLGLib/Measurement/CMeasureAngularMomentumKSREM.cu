//=============================================================================
// FILENAME : CMeasureAngularMomentumKSREM.cpp
// 
// DESCRIPTION:
// almost copy from CMeasureChiralCondensate.cpp, but with Wilson SU3 vector to SU3 vector
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3REM.h"
#include "CMeasureAngularMomentumKSREM.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAngularMomentumKSREM)

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermCopyREM(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    SSmallInt4 sCenter,
    Real fCharge)
{
    intokernalInt4;

    pResultData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
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

        deviceSU3Vector right = _deviceVXXTauOptimizedEMT(pGauge, pU1, sSite4, fCharge,
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
            pResultData[uiSiteIndex].Sub(right);
        }
        else
        {
            pResultData[uiSiteIndex].Add(right);
        }
    }

    pResultData[uiSiteIndex].MulReal(F(0.25));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermCopyREM(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCharge)
{
    intokernalInt4;

    pResultData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();

#ifndef _CLG_DTK
    #pragma unroll
#endif
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

        const deviceSU3Vector right = _deviceVXYTOptimizedEMT(
            pGauge, pU1, sSite4, fCharge,
            byGaugeFieldId, bPlusX, bPlusY, bPlusT)
            .MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            pResultData[uiSiteIndex].Sub(right);
        }
        else
        {
            pResultData[uiSiteIndex].Add(right);
        }
    }

    pResultData[uiSiteIndex].MulReal(F(0.125));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEtaCopyREM(
    BYTE byGaugeFieldId,
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3* __restrict__ pAphys,
    Real fCharge)
{
    intokernalInt4;

    BYTE byDir = 3;
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    BYTE eta_mu = (pEtaTable[uiSiteIndex] >> byDir) & 1;
    BYTE eta_mu2 = (pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir) & 1;

    deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const Real forwardPhase = pU1[linkIndex] * fCharge;
    x_Gauge_element.MulComp(_make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
    const UINT x_m_mu_link = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, x_m_mu_Gauge.m_byDir);
    const Real backPhase = pU1[x_m_mu_link] * fCharge;
    deviceSU3 x_m_mu_Gauge_element = pGauge[x_m_mu_link];
    x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(backPhase), _sin(backPhase)));

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
    //pMe[uiSiteIndex].MulReal(F(0.5));

    //Here it is gamma _4 psi, we still need r x Aphys times it
    //Multi-GPU (P4-1.6): the angular-momentum lever arms x, y must be GLOBAL;
    //identity on single-GPU (offsets 0). The site stays local for the
    //index-table lookups below.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const INT iXg = sSite4G.x;
    const INT iYg = sSite4G.y;
    const Real fY = static_cast<Real>(iYg - _DC_Centery + F(0.5));
    const Real fX = static_cast<Real>(iXg - _DC_Centerx + F(0.5));
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //x ay - y ax
    deviceSU3 midY = _deviceGetGaugeBCDirZeroT(byGaugeFieldId, pAphys, uiBigIdx, 1);
    deviceSU3 midX = _deviceGetGaugeBCDirZeroT(byGaugeFieldId, pAphys, uiBigIdx, 0);
    midY.MulReal(fX);
    midX.MulReal(fY);
    midY.Sub(midX);

    pMe[uiSiteIndex] = midY.MulVector(pMe[uiSiteIndex]);
}

#pragma endregion

void CMeasureAngularMomentumKSREM::ApplyOrbitalMatrix(
    deviceSU3Vector* pAppliedBuffer, 
    const deviceSU3Vector* pInverseZ4,
    const deviceSU3* pGauge, BYTE byGaugeFieldId) const
{
    const CFieldFermionKSSU3REM* pFieldREM = dynamic_cast<const CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(GetFermionFieldId()));
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(pFieldREM->m_byEMFieldID));

    if (NULL == pFieldREM)
    {
        appCrucial(_T("CMeasureAngularMomentumKSREM must work with CFieldFermionKSSU3REM!!!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermCopyREM, block, threads, 
        pInverseZ4,
        pGauge,
        pU1->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pAppliedBuffer,
        GetFermionFieldId(),
        byGaugeFieldId,
        _HC_Center,
        pFieldREM->m_fQ);
}

void CMeasureAngularMomentumKSREM::ApplySpinMatrix(
    deviceSU3Vector* pAppliedBuffer, 
    const deviceSU3Vector* pInverseZ4, 
    const deviceSU3* pGauge, BYTE byGaugeFieldId) const
{
    const CFieldFermionKSSU3REM* pFieldREM = dynamic_cast<const CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(GetFermionFieldId()));
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(pFieldREM->m_byEMFieldID));

    if (NULL == pFieldREM)
    {
        appCrucial(_T("CMeasureAngularMomentumKSREM must work with CFieldFermionKSSU3REM!!!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermCopyREM, block, threads, 
        pInverseZ4,
        pGauge,
        pU1->m_pDeviceData,
        pAppliedBuffer,
        GetFermionFieldId(),
        byGaugeFieldId,
        pFieldREM->m_fQ);
}

void CMeasureAngularMomentumKSREM::ApplyPotentialMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge, BYTE byGaugeFieldId) const
{
    const CFieldFermionKSSU3REM* pFieldREM = dynamic_cast<const CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(GetFermionFieldId()));
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(pFieldREM->m_byEMFieldID));

    const CFieldGaugeSU3* pAphys = dynamic_cast<const CFieldGaugeSU3*>(appGetLattice()->m_pAphys);
    if (NULL == pAphys)
    {
        appCrucial(_T("CMeasureAMomentumStochastic: A phys undefined.\n"));
    }
    preparethread;
    _LAUNCH_KERNEL(_kernelKSApplyGammaEtaCopyREM, block, threads, 
        byGaugeFieldId,
        pAppliedBuffer,
        pInverseZ4,
        pGauge,
        pU1->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[GetFermionFieldId()],
        appGetLattice()->m_pIndexCache->m_pMoveCache[GetFermionFieldId()],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pAphys->m_pDeviceData,
        pFieldREM->m_fQ);
}

void CMeasureAngularMomentumKSREM::OnConfigurationAcceptedZ4SingleField(
    const class CFieldGauge* pAcceptGauge, 
    const class CFieldGauge* pCorrespondingStaple, 
    const class CFieldFermion* pZ4, 
    const class CFieldFermion* pInverseZ4, 
    UBOOL bStart, 
    UBOOL bEnd)
{
    CMeasureAngularMomentumKS::OnConfigurationAcceptedZ4SingleField(pAcceptGauge, pCorrespondingStaple, pZ4, pInverseZ4, bStart, bEnd);

    if (bEnd && NULL != m_pOwner)
    {
        const UINT uiConfigIdx = m_uiConfigurationCount - 1;

        m_pOwner->AddOneConfigurationResult(this, _T("OrbitalCondAll"), m_lstCondAll[OrbitalKS][uiConfigIdx]);
        m_pOwner->AddOneConfigurationResult(this, _T("SpinCondAll"), m_lstCondAll[SpinKS][uiConfigIdx]);
        m_pOwner->AddOneConfigurationResult(this, _T("PotentialCondAll"), m_lstCondAll[PotentialKS][uiConfigIdx]);

        m_pOwner->AddOneConfigurationResult(this, _T("OrbitalCondIn"), m_lstCondIn[OrbitalKS][uiConfigIdx]);
        m_pOwner->AddOneConfigurationResult(this, _T("SpinCondIn"), m_lstCondIn[SpinKS][uiConfigIdx]);
        m_pOwner->AddOneConfigurationResult(this, _T("PotentialCondIn"), m_lstCondIn[PotentialKS][uiConfigIdx]);

        const UINT uiRCount = m_lstR.Num();
        if (uiRCount > 0)
        {
            const UINT uiRStart = uiConfigIdx * uiRCount;
            TArray<CLGComplex> orbitalR;
            TArray<CLGComplex> spinR;
            TArray<CLGComplex> potentialR;
            for (UINT j = 0; j < uiRCount; ++j)
            {
                orbitalR.AddItem(m_lstCond[OrbitalKS][uiRStart + j]);
                spinR.AddItem(m_lstCond[SpinKS][uiRStart + j]);
                potentialR.AddItem(m_lstCond[PotentialKS][uiRStart + j]);
            }
            m_pOwner->AddOneConfigurationResult(this, _T("OrbitalCondR"), orbitalR);
            m_pOwner->AddOneConfigurationResult(this, _T("SpinCondR"), spinR);
            m_pOwner->AddOneConfigurationResult(this, _T("PotentialCondR"), potentialR);
        }

        if (m_bMeasureZSlice)
        {
            const UINT uiZSliceCount = _HC_Lz;
            const UINT uiZSliceStart = uiConfigIdx * uiZSliceCount;
            TArray<CLGComplex> orbitalZSlice;
            TArray<CLGComplex> spinZSlice;
            TArray<CLGComplex> potentialZSlice;
            for (UINT j = 0; j < uiZSliceCount; ++j)
            {
                orbitalZSlice.AddItem(m_lstCondZSlice[OrbitalKS][uiZSliceStart + j]);
                spinZSlice.AddItem(m_lstCondZSlice[SpinKS][uiZSliceStart + j]);
                potentialZSlice.AddItem(m_lstCondZSlice[PotentialKS][uiZSliceStart + j]);
            }
            m_pOwner->AddOneConfigurationResult(this, _T("OrbitalZSlice"), orbitalZSlice);
            m_pOwner->AddOneConfigurationResult(this, _T("SpinZSlice"), spinZSlice);
            m_pOwner->AddOneConfigurationResult(this, _T("PotentialZSlice"), potentialZSlice);
        }
    }
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================