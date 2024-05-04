//=============================================================================
// FILENAME : CFieldFermionKSSU3DAcc.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/19/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3Acc)

#pragma region DOperator

#pragma region kernel

/**
* gz * gamma3 (pt + iAt)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSACC_GZPT(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fG,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    const INT eta_z = pEtaTable[uiSiteIndex] >> 2;

    #pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlusZ = idx & 1;
        const UBOOL bPlusT = idx & 2;
        SSmallInt4 sTargetSite = sSite4;

        //gamma_z partial_t
        sTargetSite.z = sTargetSite.z + (bPlusZ ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusT ? 2 : -2);
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        sTargetSite = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        INT this_eta_z = eta_z;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_z = this_eta_z + 1;
        }

        //gamma_z partial_t
        //deviceSU3Vector right = _deviceGMUPNU(pGauge, sSite4, byGaugeFieldId, 2, 3, bPlusZ, bPlusT).MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        deviceSU3Vector right = _deviceGMUPNUOptimized(pGauge, sSite4, byGaugeFieldId, 2, 3, bPlusZ, bPlusT).MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //right.MulComp(_make_cuComplex(F(0.0), fG * sSite4.z));
        right.MulReal(F(0.5) * (sSite4.w + sTargetSite.w));

        if (!bPlusT)
        {
            //for -2t terms, there is another minus sign
            this_eta_z = this_eta_z + 1;
        }

        if (this_eta_z & 1)
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
        result.MulReal(fG * F(-0.25));
    }
    else
    {
        result.MulReal(fG * F(0.25));
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
_kernelDFermionKSForce_Acc(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    Real fG,
    INT iZ,
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
    SSmallInt4 site_n2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n2)];

    site_n1 = __deviceSiteIndexToInt4(sn1.m_uiSiteIndex);
    site_n2 = __deviceSiteIndexToInt4(sn2.m_uiSiteIndex);
    const Real fNv = fG * F(0.5) * (site_n1.w + site_n2.w);

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

        //if (sn1.NeedToOpposite())
        //{
        //    phi1.MulReal(F(-1.0));
        //    phi3.MulReal(F(-1.0));
        //}
        //if (sn2.NeedToOpposite())
        //{
        //    phi2.MulReal(F(-1.0));
        //    phi4.MulReal(F(-1.0));
        //}
        deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        INT etaTau = (sn1.NeedToOpposite() ^ sn2.NeedToOpposite());
        etaTau += (pEtaTable[sn1.m_uiSiteIndex] >> 2);

        if (etaTau & 1)
        {
            res.MulReal(-OneOver12 * fNv * pNumerators[rfieldId]);
        }
        else
        {
            res.MulReal(OneOver12 * fNv * pNumerators[rfieldId]);
        }

        //partial term
        if (0 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            pForce[linkIndex].Sub(res);
        }

        //gamma term
        if (1 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);
            if (iZ > 0)
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


#pragma endregion


void CFieldFermionKSSU3Acc::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionKSSU3::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKSACC_GZPT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        1,
        CCommonData::m_fG,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
void CFieldFermionKSSU3Acc::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{

    CFieldFermionKSSU3::DerivateD0(pForce, pGaugeBuffer);

    preparethread;

    INT dirs[6][3] =
    {
        {3, 4, 4},
        {4, 3, 4},
        {4, 4, 3},
        {4, 4, -3},
        {4, -3, 4},
        {-3, 4, 4},
    };

    INT iZ[6] = { 1, 1, 1, -1, -1, -1 };
    BYTE contributionOf[6][4] =
    {
        {1, 0, 0, 3},
        {0, 1, 0, 3},
        {0, 0, 1, 3},
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

            _kernelDFermionKSForce_Acc << <block, threads >> > (
                (const deviceSU3*)pGaugeBuffer,
                (deviceSU3*)pForce,
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                m_pRationalFieldPointers,
                m_pMDNumerator,
                m_rMD.m_uiDegree,
                m_byFieldId,
                CCommonData::m_fG,
                iZ[pathidx],
                L[0], L[1], L[2], LLength,
                R[0], R[1], R[2], RLength,
                contributionOf[pathidx][iSeperation]
                );
        }
    }
}

#pragma endregion

void CFieldFermionKSSU3Acc::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

CCString CFieldFermionKSSU3Acc::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3Acc\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("G : ") + appToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================