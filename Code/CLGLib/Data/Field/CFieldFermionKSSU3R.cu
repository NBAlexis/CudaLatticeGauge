//=============================================================================
// FILENAME : CFieldFermionKSSU3R.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/05/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3R)

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
    UBOOL bXorY,
    Real fOmega,
    SSmallInt4 sCenter,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    Real eta_tau = (1 == ((pEtaTable[uiSiteIndex] >> 3) & 1)) ? F(-1.0) : F(1.0);
    const Real fNvOmega = (bXorY
        ? static_cast<Real>(sSite4.y - sCenter.y)
        : static_cast<Real>(sCenter.x - sSite4.x)) * fOmega;


    #pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlusMu  = (0 != (idx & 1));
        const UBOOL bPlusTau = (0 != (idx & 2));

        deviceSU3Vector right = _deviceOffsetXXTau(pDeviceData, sSite4, byFieldId, bXorY, bPlusMu, bPlusTau);
        right = _deviceVXXTau(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau).MulVector(right);
        
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
        eta_tau = eta_tau * F(-1.0);
    }
    result.MulReal(eta_tau * F(0.25) * fNvOmega);

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
    Real fOmega,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    Real eta124 = _deviceEta124(sSite4);

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        deviceSU3Vector right = _deviceOffsetXYTau(pDeviceData, sSite4, byFieldId, bPlusX, bPlusY, bPlusT);
        right = _deviceVXYT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT).MulVector(right);

        result.Add(right);
    }

    if (bDDagger)
    {
        eta124 = eta124 * F(-1.0);
    }
    result.MulReal(eta124 * F(0.125) * fOmega);

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
    Real fOmega, SSmallInt4 sCenter, BYTE byMu, INT iTauSign,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength,
    INT iEtaShift, BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
    ? static_cast<Real>(sSite4.y - sCenter.y)
    : static_cast<Real>(sCenter.x - sSite4.x);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SIndex sn1 = __idx->_deviceGetMappingIndex(
        _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength), byFieldId);
    const SIndex sn2 = __idx->_deviceGetMappingIndex(
        _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength), byFieldId);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE byFieldId = 0; byFieldId < uiRational; ++byFieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[byFieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[byFieldId + uiRational];
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
        const Real eta_tau = (1 == ((pEtaTable[uiSiteIndex] >> 3) & 1)) ? F(-1.0) : F(1.0);
        res.MulReal(OneOver12 * fOmega * fNv * pNumerators[byFieldId] * eta_tau * iEtaShift);

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
            if (iTauSign > 0)
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
_kernelDFermionKSForce_R_XYTau_Term(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    Real fOmega, 
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex sn1 = __idx->_deviceGetMappingIndex(siten1, byFieldId);
    const SIndex sn2 = __idx->_deviceGetMappingIndex(siten2, byFieldId);
    const Real eta124 = _deviceEta124(siten1);
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE byFieldId = 0; byFieldId < uiRational; ++byFieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[byFieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[byFieldId + uiRational];

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
        res.MulReal(OneOver48 * fOmega * pNumerators[byFieldId] * eta124);

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

#pragma endregion


#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3R::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionKSSU3::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

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
        1,
        TRUE,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
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
        1,
        FALSE,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);


    _kernelDFermionKS_R_XYTau_Term << <block, threads >> > (
        pSource,
        pGauge,
        pTarget,
        m_byFieldId,
        1,
        CCommonData::m_fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

}

void CFieldFermionKSSU3R::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{
    CFieldFermionKSSU3::DerivateD0(pForce, pGaugeBuffer);

    preparethread;

    INT mu[2] = { 0, 1 };
    //Test only x, not y for now
    for (INT imu = 0; imu < 2; ++imu)
    {
        INT dirs[6][3] =
        {
            {4, mu[imu] + 1, mu[imu] + 1},
            {mu[imu] + 1, 4, mu[imu] + 1},
            {mu[imu] + 1, mu[imu] + 1, 4},
            {mu[imu] + 1, mu[imu] + 1, -4},
            {mu[imu] + 1, -4, mu[imu] + 1},
            {-4, mu[imu] + 1, mu[imu] + 1}
        };

        INT iTaus[6] = {1, 1, 1, -1, -1, -1};

        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 3},
            {0, 1, 0, 3},
            {0, 0, 1, 3},
            {0, 0, 3, 1},
            {0, 3, 2, 3},
            {3, 2, 0, 3}
        };

        INT etashift[6][4] =
        {
            {1, 1, -1, 1},
            {1, -1, -1, 1},
            {1, -1, 1, 1},
            {1, -1, 1, 1},
            {1, -1, -1, 1},
            {1, 1, -1, 1}
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            INT iTau = iTaus[pathidx];
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
                    CCommonData::m_fOmega, CCommonData::m_sCenter,
                    static_cast<BYTE>(imu), iTau,
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    etashift[pathidx][iSeperation],
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }


    //===========================
    //polarization terms
    
    //ilinkType is +-x +-y +t,
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
                        CCommonData::m_fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }
}

#pragma endregion

void CFieldFermionKSSU3R::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

CCString CFieldFermionKSSU3R::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3R\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(CCommonData::m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================