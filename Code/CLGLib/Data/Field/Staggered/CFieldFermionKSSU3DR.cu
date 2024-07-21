//=============================================================================
// FILENAME : CFieldFermionKSSU3DR.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/05/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineSU3.h"
#include "CFieldFermionKSSU3DR.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3DR)

#pragma region DOperator

#pragma region kernel

/**
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_R_XYTerm(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const BYTE * __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bShiftCenter,
    UBOOL bXorY,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIndex = __idx->_deviceGetBigIndex(sSite4);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        pResultData[uiSiteIndex] = result;
        return;
    }
    const Real shift = bShiftCenter ? F(0.5) : F(0.0);
    #pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlusMu  = (0 != (idx & 1));
        const UBOOL bPlusTau = (0 != (idx & 2));
        SSmallInt4 sTargetSite = sSite4;
        if (bXorY)
        {
            sTargetSite.x = sTargetSite.x + (bPlusMu ? 2 : -2);
        }
        else
        {
            sTargetSite.y = sTargetSite.y + (bPlusMu ? 2 : -2);
        }
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        Real eta_tau = F(1.0);
        if (bPlusTau)
        {
            if (1 == ((pEtaTable[uiSiteIndex] >> 3) & 1))
            {
                eta_tau = F(-1.0);
            }
            if (bXorY)
            {
                eta_tau = eta_tau * (sSite4.y - _DC_Centery + shift);
            }
            else
            {
                eta_tau = eta_tau * (_DC_Centerx - sSite4.x - shift);
            }
        }
        else
        {
            if (1 == ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
            {
                eta_tau = F(-1.0);
            }
            //sTargetSite no longer use
            //Why to use target set?????
            sTargetSite = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);
            if (bXorY)
            {
                eta_tau = eta_tau * (sTargetSite.y - _DC_Centery + shift);
            }
            else
            {
                eta_tau = eta_tau * (_DC_Centerx - sTargetSite.x - shift);
            }
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta_tau = eta_tau * F(-1.0);
        }

        deviceSU3Vector right = _deviceVXXTau(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau).MulVector(
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        right.MulReal(eta_tau);

        if (bPlusMu)
        {
            result.Add(right);
        }
        else
        {
            result.Sub(right);
        }
    }

    if (bDDagger)
    {
        fOmega = -fOmega;
    }
    result.MulReal(F(0.25) * fOmega);

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
_kernelDFermionKS_R_XYTau_Term(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    const UINT uiBigIndex = __idx->_deviceGetBigIndex(sSite4);
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        pResultData[uiSiteIndex] = result;
        return;
    }
    const Real eta124 = _deviceEta124(sSite4);

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
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        
        deviceSU3Vector right = _deviceVXYT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT)
        .MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        if (sTargetBigIndex.NeedToOpposite())
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
        result.MulReal(-F(0.125) * fOmega * eta124);
    }
    else
    {
        result.MulReal(F(0.125) * fOmega * eta124);
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

/**
 * Have n, n->n1, n->n2,
 * 1. we need to obtain V_(n, n1) , V_(n, n2)
 * 2. we need phi(n1), phi(n2), phid(n1), phid(n2)
 *
 * byContribution: 0 for mu, 1 for tau, 2 for both mu and tau
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_R_XYTerm( 
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bShiftCenter,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byMu, INT iMu,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Ldirs, Llength))];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    {
        return;
    }
    const SSmallInt4 site_n1 = __deviceSiteIndexToInt4(sn1.m_uiSiteIndex);
    const Real shift = bShiftCenter ? F(0.5) : F(0.0);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + shift)
        : static_cast<Real>(_DC_Centerx - site_n1.x - shift);

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
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        const Real eta_tau = (1 == ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)) ? F(-1.0) : F(1.0);
        res.MulReal(OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);

        //For mu
        if (0 == byContribution || 2 == byContribution)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, byMu))
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
                pForce[linkIndex].Sub(res);
            }
        }

        //For tau
        if (1 == byContribution || 2 == byContribution)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, 3))
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
                if (iMu > 0)
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
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_R_XYTau_Term(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    {
        return;
    }
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

        if (pathLdir1 > 0)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, pathLdir1 - 1))
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);
                pForce[linkIndex].Sub(res);
            }
        }

        if (pathRdir1 > 0)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, pathRdir1 - 1))
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);
                pForce[linkIndex].Add(res);
            }
        }
    }

}

#pragma endregion


#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3DR::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionKSSU3D::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKS_R_XYTerm << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        byGaugeFieldId,
        m_bShiftCenter,
        TRUE,
        m_fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionKS_R_XYTerm << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        byGaugeFieldId,
        m_bShiftCenter,
        FALSE,
        m_fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionKS_R_XYTau_Term << <block, threads >> > (
        pSource,
        pGauge,
        pTarget,
        m_byFieldId,
        byGaugeFieldId,
        m_fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionKSSU3DR::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer,
    BYTE byGaugeFieldId) const
{
    CFieldFermionKSSU3D::DerivateD0(pForce, pGaugeBuffer, byGaugeFieldId);

    preparethread;

    #pragma region X Y Term

    INT mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        INT dirs[6][3] =
        {
            {4, mu[imu] + 1, mu[imu] + 1},
            {mu[imu] + 1, 4, mu[imu] + 1},
            {mu[imu] + 1, mu[imu] + 1, 4},
            {4, -mu[imu] - 1, -mu[imu] - 1},
            {-mu[imu] - 1, 4, -mu[imu] - 1},
            {-mu[imu] - 1, -mu[imu] - 1, 4},
        };

        INT iMu[6] = { 1, 1, 1, -1, -1, -1 };
        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 3},
            {0, 1, 0, 3},
            {0, 0, 1, 3},
            {1, 3, 0, 0},
            {3, 2, 3, 0},
            {3, 0, 2, 3},
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

                _kernelDFermionKSForce_R_XYTerm << <block, threads >> > (
                    (const deviceSU3*)pGaugeBuffer,
                    (deviceSU3*)pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    m_pRationalFieldPointers,
                    m_pMDNumerator,
                    m_rMD.m_uiDegree,
                    m_byFieldId,
                    byGaugeFieldId,
                    m_bShiftCenter,
                    m_fOmega,
                    static_cast<BYTE>(imu), iMu[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

    #pragma endregion

    #pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    INT linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, 2, -4},
        {-1, 2, 4},
        {-1, 2, -4}
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
                    _kernelDFermionKSForce_R_XYTau_Term << <block, threads >> > (
                        (const deviceSU3*)pGaugeBuffer,
                        (deviceSU3*)pForce,
                        m_pRationalFieldPointers,
                        m_pMDNumerator,
                        m_rMD.m_uiDegree,
                        m_byFieldId,
                        byGaugeFieldId,
                        m_fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

    #pragma endregion
}

#pragma endregion

void CFieldFermionKSSU3DR::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3D::InitialOtherParameters(params);

    DOUBLE fValue = 0.1;
    if (params.FetchValueDOUBLE(_T("Omega"), fValue))
    {
        m_fOmega = fValue;
    }

    INT iShiftCenter = 0;
    if (params.FetchValueINT(_T("ShiftCenter"), iShiftCenter))
    {
        m_bShiftCenter = (0 != iShiftCenter);
    }
}

void CFieldFermionKSSU3DR::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
    CFieldFermionKSSU3DR* pOther = dynamic_cast<CFieldFermionKSSU3DR*>(U);
    if (NULL != pOther)
    {
        pOther->m_bShiftCenter = m_bShiftCenter;
        pOther->m_fOmega = m_fOmega;
    }
}

CCString CFieldFermionKSSU3DR::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKSSU3D::GetInfos(tab);
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    sRet = sRet + tab + _T("ShiftCenter : ") + appToString(m_bShiftCenter) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================