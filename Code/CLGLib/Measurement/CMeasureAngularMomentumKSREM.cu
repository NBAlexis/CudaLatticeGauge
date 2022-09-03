//=============================================================================
// FILENAME : CMeasureAngularMomentumKS.cpp
// 
// DESCRIPTION:
// almost copy from CMeasureChiralCondensate.cpp, but with Wilson SU3 vector to SU3 vector
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAngularMomentumKSREM)

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermCopyREM(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    SSmallInt4 sCenter,
    Real fQBz,
    BYTE byGaugeType,
    UBOOL bTwisted)
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

        deviceSU3Vector right = _deviceVXXTauOptimizedEM(pGauge, sSite4, sCenter, fQBz, byGaugeType, bTwisted,
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
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    SSmallInt4 sCenter,
    Real fQBz,
    BYTE byGaugeType,
    UBOOL bTwisted)
{
    intokernalInt4;

    pResultData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();

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
            pGauge, sSite4, sCenter, fQBz, byGaugeType, bTwisted,
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
            pResultData[uiSiteIndex].Sub(right);
        }
        else
        {
            pResultData[uiSiteIndex].Add(right);
        }
    }

    pResultData[uiSiteIndex].MulReal(F(0.125));
}

#pragma endregion

void CMeasureAngularMomentumKSREM::ApplyOrbitalMatrix(
    deviceSU3Vector* pAppliedBuffer, 
    const deviceSU3Vector* pInverseZ4,
    const deviceSU3* pGauge) const
{
    const CFieldFermionKSSU3REM* pFieldREM = dynamic_cast<const CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(m_byFieldId));

    if (NULL == pFieldREM)
    {
        appCrucial(_T("CMeasureAngularMomentumKSREM must work with CFieldFermionKSSU3REM!!!\n"));
        return;
    }

    preparethread;
    _kernelDFermionKS_PR_XYTermCopyREM << <block, threads >> > (
        pInverseZ4,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pAppliedBuffer,
        m_byFieldId,
        1,
        CCommonData::m_sCenter,
        CCommonData::m_fBz * pFieldREM->GetQ(),
        m_byGaugeType,
        m_bTwistedBoundary);
}

void CMeasureAngularMomentumKSREM::ApplySpinMatrix(
    deviceSU3Vector* pAppliedBuffer, 
    const deviceSU3Vector* pInverseZ4, 
    const deviceSU3* pGauge) const
{
    const CFieldFermionKSSU3REM* pFieldREM = dynamic_cast<const CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(m_byFieldId));

    if (NULL == pFieldREM)
    {
        appCrucial(_T("CMeasureAngularMomentumKSREM must work with CFieldFermionKSSU3REM!!!\n"));
        return;
    }

    preparethread;
    _kernelDFermionKS_PR_XYTau_TermCopyREM << <block, threads >> > (
        pInverseZ4,
        pGauge,
        pAppliedBuffer,
        m_byFieldId,
        1,
        CCommonData::m_sCenter,
        CCommonData::m_fBz * pFieldREM->GetQ(),
        m_byGaugeType,
        m_bTwistedBoundary);
}

void CMeasureAngularMomentumKSREM::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& params, BYTE byId)
{
    CMeasureAngularMomentumKS::Initial(pOwner, pLatticeData, params, byId);

    INT iValue = 0;
    params.FetchValueINT(_T("MagneticType"), iValue);
    m_byGaugeType = static_cast<BYTE>(iValue);

    iValue = 1;
    params.FetchValueINT(_T("TwistedBoundary"), iValue);
    m_bTwistedBoundary = (0 != iValue);
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================