//=============================================================================
// FILENAME : CFieldFermionKSSU3Gamma.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/19/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CFieldFermionKSSU3Gamma.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3Gamma)

#pragma region gamma kernels

/**
 * gamma_i ^dagger = gamma_i for i = 1,2,3,4,5, so 'bDDagger' is used for only coefficient
 *
 * 1/2a eta_mu(x) sum_{mu=+-} psibar(x) U_mu psi(x+mu)
 * the 1/2a is absorbed (be sure to add when measuring)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernal;

    const Real eta_mu = ((pEtaTable[uiSiteIndex] >> byDir) & 1) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceSU3Vector result = x_Gauge_element.MulVector(pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        result.MulReal(F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        result.Add(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    //result.MulReal(eta_mu);

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    result.MulComp(_make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pMe[uiSiteIndex].Add(result);
}

/**
* This function applys i Gamma_i on the field
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEta1234(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernal;

    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    BYTE eta_mu = pEtaTable[uiSiteIndex] >> byDir;
    BYTE eta_mu2 = pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir;

    const deviceSU3& x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, x_m_mu_Gauge.m_byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceSU3Vector result = x_Gauge_element.MulVector(pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        result.MulReal(F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        result.Sub(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        result.Add(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    result.MulComp(_make_cuComplex(F(0.0), fGammCoefficient));

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pMe[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSForceGamma1234(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE idir,
    Real fCoeff,
    BYTE byFieldId)
{
    intokernal;
    
    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];

    for (UINT uiR = 0; uiR < uiRational; ++uiR)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[uiR];
        const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];

        deviceSU3Vector toContract = pGauge[linkIndex].MulVector(phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

        toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
        thisTerm.Sub(deviceSU3::makeSU3ContractV(toContract, phi_i[uiSiteIndex]));

        if (x_p_mu_Fermion.NeedToOpposite())
        {
            thisTerm.MulReal(eta_mu * pNumerators[uiR] * F(-1.0));
        }
        else
        {
            thisTerm.MulReal(eta_mu * pNumerators[uiR]);
        }

        thisTerm.MulComp(_make_cuComplex(F(0.0), fCoeff));
        thisTerm.Ta();


        pForce[linkIndex].Sub(thisTerm);
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaSigmaIJ(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    BYTE byEtaShift,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    SBYTE byDir1,
    SBYTE byDir2,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlus12[2] = { (0 != (idx & 1)), (0 != (idx & 2)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim12[2] = {
            bPlus12[0] ? static_cast<SBYTE>(byDir1 + 1) : static_cast<SBYTE>(-byDir1 - 1),
            bPlus12[1] ? static_cast<SBYTE>(byDir2 + 1) : static_cast<SBYTE>(-byDir2 - 1)
        };
        sOffset.m_byData4[byDir1] = sOffset.m_byData4[byDir1] + (bPlus12[0] ? 1 : -1);
        sOffset.m_byData4[byDir2] = sOffset.m_byData4[byDir2] + (bPlus12[1] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceSU3Vector right = _devicePlaneDiagonal(pGauge, sSite4, byGaugeFieldId, dim12[0], dim12[1])
            .MulVector(pOther[sTargetBigIndex.m_uiSiteIndex]);

        //eta12 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta2 = _deviceEta2(pEtaTable[uiSiteIndex], byDir1, byDir2);
        if (1 == byEtaShift || 2 == byEtaShift)
        {
            //type 1, it is not the corner, and if it cross the Y-boundary
            if (!bPlus12[byEtaShift - 1])
            {
                eta2 = _deviceEta2(pEtaTable[sTargetBigIndex.m_uiSiteIndex], byDir1, byDir2) + 1;
            }
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta2 = eta2 + 1;
        }

        if (eta2 & 1)
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
        result.MulReal(-F(0.5) * fGammCoefficient);
    }
    else
    {
        result.MulReal(F(0.5) * fGammCoefficient);
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

    pMe[uiSiteIndex].Add(result);
}

/**
 * gamma5i corresponds to diagonal links of cubic in the other 3 dimensions.
 * For each diagonal link, there are 6 different ways to add the gauge field
 * We simply use the average of all 6 links
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma51234(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    UBOOL bDDagger,
    BYTE byMissingDir,
    SBYTE byEtaShift,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll

    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlus123[3] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim123[3] = {0, 0, 0};
        BYTE byDimIndex = 0;
        for (SBYTE byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
        {
            if (byCubeDir != static_cast<SBYTE>(byMissingDir))
            {
                sOffset.m_byData4[byCubeDir] = sOffset.m_byData4[byCubeDir] + (bPlus123[byDimIndex] ? 1 : -1);
                dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                byDimIndex++;
            }
        }

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceSU3Vector right = _deviceCubicDiagonal(pGauge, sSite4, byGaugeFieldId, dim123[0], dim123[1], dim123[2])
            .MulVector(pOther[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta3 = _deviceEta3(sSite4, byMissingDir);
        if (byEtaShift >= 0 && byEtaShift <= 2 && !bPlus123[byEtaShift])
        {
            eta3 = _deviceEta3(site_target, byMissingDir) + 1;
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta3 = eta3 + 1;
        }

        if (eta3 & 1)
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
        result.MulReal(-F(0.25) * fGammCoefficient);
    }
    else
    {
        result.MulReal(F(0.25) * fGammCoefficient);
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

    pMe[uiSiteIndex].Add(result);
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma5(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    UBOOL bDDagger,
    UBOOL bEtaShift,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 16; ++idx)
    {
        const UBOOL bPlus1234[4] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)), (0 != (idx & 8)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim1234[4] = 
        {
            bPlus1234[0] ? static_cast<SBYTE>(1) : static_cast<SBYTE>(-1),
            bPlus1234[1] ? static_cast<SBYTE>(2) : static_cast<SBYTE>(-2),
            bPlus1234[2] ? static_cast<SBYTE>(3) : static_cast<SBYTE>(-3),
            bPlus1234[3] ? static_cast<SBYTE>(4) : static_cast<SBYTE>(-4)
        };
        sOffset.m_byData4[0] = sOffset.m_byData4[0] + (bPlus1234[0] ? 1 : -1);
        sOffset.m_byData4[1] = sOffset.m_byData4[1] + (bPlus1234[1] ? 1 : -1);
        sOffset.m_byData4[2] = sOffset.m_byData4[2] + (bPlus1234[2] ? 1 : -1);
        sOffset.m_byData4[3] = sOffset.m_byData4[3] + (bPlus1234[3] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceSU3Vector right = _deviceHyperCubicDiagonal(pGauge, sSite4, byGaugeFieldId, dim1234[0], dim1234[1], dim1234[2], dim1234[3])
            .MulVector(pOther[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta51 is gamma5 (x+z)
        INT eta4 = _deviceEta3(sSite4, 0);
        if (bEtaShift && !bPlus1234[3])
        {
            //target is almost always site4 except for boundaries
            eta4 = _deviceEta3(site_target, 0);
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta4 = eta4 + 1;
        }

        if (eta4 & 1)
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
        result.MulComp(_make_cuComplex(F(0.0), -F(0.125) * fGammCoefficient));
    }
    else
    {
        result.MulComp(_make_cuComplex(F(0.0), F(0.125) * fGammCoefficient));
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

    pMe[uiSiteIndex].Add(result);
}

/**
 * similar as _kernelDFermionKSForce_WithLink.
 * but _kernelDFermionKSForce_WithLink use gamma mu as the gamma matrix
 * we use sigma12
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_SigmaIJ(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE bDir1,
    BYTE bDir2,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[2];
    INT pathRight[2];
    for (BYTE iSeperation = 0; iSeperation <= 2; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 2, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta2(pEtaTable[sn1.m_uiSiteIndex], bDir1, bDir2);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
                const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
                deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
                //deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                //deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);

                deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
                phi2 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                res.Add(deviceSU3::makeSU3ContractV(phi1, phi2));
                //res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
                res.Ta();
                res.MulReal(fGammCoefficient * pNumerators[rfieldId]);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Sub(res);
                    }
                    else
                    {
                        pForce[linkIndex].Add(res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Add(res);
                    }
                    else
                    {
                        pForce[linkIndex].Sub(res);
                    }
                }
            }
        }
    }
}

/**
 * similar as _kernelDFermionKSForce_WithLink.
 * but _kernelDFermionKSForce_WithLink use gamma mu as the gamma matrix
 * we use gamma5i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma51234(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE byMissingDir,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[3];
    INT pathRight[3];
    for (BYTE iSeperation = 0; iSeperation <= 3; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 3, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta3(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex), byMissingDir);
            //INT iEtaMu2 = _deviceEta3(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex), byMissingDir);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
                const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
                deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
                //deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                //deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);

                deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
                phi2 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                res.Add(deviceSU3::makeSU3ContractV(phi1, phi2));
                //res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
                res.Ta();
                res.MulReal(fGammCoefficient * pNumerators[rfieldId]);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Sub(res);
                    }
                    else
                    {
                        pForce[linkIndex].Add(res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Add(res);
                    }
                    else
                    {
                        pForce[linkIndex].Sub(res);
                    }
                }
            }
        }
    }
}

//similar as above
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma5(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[4];
    INT pathRight[4];

    //if (0 == uiSiteIndex)
    //{
    //    printf("%f, %d %d %d %d\n", fGammCoefficient, path[0], path[1], path[2], path[3]);
    //}

    for (BYTE iSeperation = 0; iSeperation <= 4; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 4, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta3(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex), 0);
            //INT iEtaMu2 = _deviceEta3(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex), byMissingDir);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
                const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
                deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
                //deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                //deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);

                deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
                phi2 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                res.Sub(deviceSU3::makeSU3ContractV(phi1, phi2));
                //res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));

                res.MulComp(_make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));
                //res.MulReal(fGammCoefficient * pNumerators[rfieldId]);
                res.Ta();

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Sub(res);
                    }
                    else
                    {
                        pForce[linkIndex].Add(res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].Add(res);
                    }
                    else
                    {
                        pForce[linkIndex].Sub(res);
                    }
                }
            }
        }
    }
}


#pragma endregion



void CFieldFermionKSSU3Gamma::appApplyGammaKS(
    void* pTargetBuffer,
    const void* pBuffer,
    const void* pGaugeBuffer,
    EGammaMatrix eGamma,
    UBOOL bShiftCenter,
    UBOOL bDagger,
    Real fGammaCoeff,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    preparethread;

    switch (eGamma)
    {
    case GAMMA1:
    case GAMMA2:
    case GAMMA3:
    case GAMMA4:
        {
            INT iDir = static_cast<INT>(eGamma) - 1;
            
            if (bShiftCenter)
            {
                _kernelKSApplyGammaEta1234 << <block, threads >> > (
                    pTarget,
                    pSource,
                    pGauge,
                    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    bDagger,
                    fGammaCoeff,
                    eOCT,
                    fRealCoeff,
                    cCmpCoeff,
                    static_cast<BYTE>(iDir));
            }
            else
            {
                _kernelKSApplyGamma1234 << <block, threads >> > (
                    pTarget,
                    pSource,
                    pGauge,
                    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    bDagger,
                    fGammaCoeff,
                    eOCT,
                    fRealCoeff,
                    cCmpCoeff,
                    static_cast<BYTE>(iDir));
            }
        }
        break;
    case SIGMA12:
        if (bShiftCenter)
        {
            appCrucial(_T("Sigma 12 in projective plane boundary condition is not supported!\n"));
        }
        else
        {
            _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
                pTarget,
                pSource,
                pGauge,
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fGammaCoeff,
                bShiftCenter ? 3 : 0,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                0,
                1,
                byFieldID,
                byGaugeFieldID);
        }
        break;
    case SIGMA31:
        //this is sigma 13
        _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 1 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            0,
            2,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA41:
        //this is sigma 14
        _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 1 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            0,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA23:
        _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 2 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            1,
            2,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA42:
        //this is sigma 24
        _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 2 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            1,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA43:
        //this is sigma 34
        _kernelKSApplyGammaSigmaIJ << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            2,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case GAMMA51:
    case GAMMA52:
    case GAMMA53:
    case GAMMA54:
        {
            //eta shift is:
            //x->y
            //y->x
            //z->t
            //t->z

            const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
            SBYTE etaShift = -1;
            if (bShiftCenter)
            {
                if (byMissingDir < 2)
                {
                    etaShift = 0;
                }
                else
                {
                    etaShift = 2;
                }
            }

            _kernelKSApplyGamma51234 << <block, threads >> > (
                pTarget,
                pSource,
                pGauge,
                bDagger,
                byMissingDir,
                etaShift,
                fGammaCoeff,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                byFieldID,
                byGaugeFieldID);
        }
        break;
    case GAMMA5:
        _kernelKSApplyGamma5 << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            bDagger,
            bShiftCenter,
            fGammaCoeff,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            byFieldID,
            byGaugeFieldID);
        break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}

void CFieldFermionKSSU3Gamma::GammaKSForce(
    void* pForce,
    const void* pGaugeBuffer,
    const deviceSU3Vector* const* pRationalFields,
    const Real* pRationalNumerator,
    UINT uiRationalDegree,
    Real fCoeff,
    EGammaMatrix eGamma,
    INT* devicePathBuffer,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    preparethread;

    switch (eGamma)
    {
    case GAMMA1:
    case GAMMA2:
    case GAMMA3:
    case GAMMA4:
    {
        BYTE byDir = static_cast<BYTE>(eGamma) - 1;
        _kernelFermionKSForceGamma1234 << <block, threads >> > (
            (const deviceSU3*)pGaugeBuffer,
            (deviceSU3*)pForce,
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pRationalFields,
            pRationalNumerator,
            uiRationalDegree,
            byDir,
            fCoeff,
            byFieldID);
    }
    break;
    case SIGMA12:
    case SIGMA31:
    case SIGMA41:
    case SIGMA23:
    case SIGMA42:
    case SIGMA43:
        {
            SBYTE byDirs[2] = { 0, 1 };
            if (SIGMA31 == eGamma)
            {
                byDirs[0] = 0; byDirs[1] = 2;
            }
            else if (SIGMA41 == eGamma)
            {
                byDirs[0] = 0; byDirs[1] = 3;
            }
            else if (SIGMA23 == eGamma)
            {
                byDirs[0] = 1; byDirs[1] = 2;
            }
            else if (SIGMA42 == eGamma)
            {
                byDirs[0] = 1; byDirs[1] = 3;
            }
            else if (SIGMA43 == eGamma)
            {
                byDirs[0] = 2; byDirs[1] = 3;
            }
            for (INT idx = 0; idx < 2; ++idx)
            {
                UBOOL bPlus2 = (0 == (idx & 1));

                INT dimordered[2];
                for (INT order = 0; order < 2; ++order)
                {
                    dimordered[0] = byDirs[0] + 1;
                    dimordered[1] = bPlus2 ? static_cast<SBYTE>(byDirs[1] + 1) : static_cast<SBYTE>(-byDirs[1] - 1);
                    if (1 == order)
                    {
                        dimordered[0] = bPlus2 ? static_cast<SBYTE>(byDirs[1] + 1) : static_cast<SBYTE>(-byDirs[1] - 1);
                        dimordered[1] = byDirs[0] + 1;
                    }
                    checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 2, cudaMemcpyHostToDevice));

                    _kernelDFermionKSForce_WithLink_SigmaIJ << <block, threads >> > (
                        (const deviceSU3*)pGaugeBuffer,
                        (deviceSU3*)pForce,
                        pRationalFields,
                        pRationalNumerator,
                        uiRationalDegree,
                        appGetLattice()->m_pIndexCache->m_pEtaMu,
                        byFieldID,
                        byGaugeFieldID,
                        fCoeff * F(0.25),
                        byDirs[0],
                        byDirs[1],
                        devicePathBuffer
                        );
                }
            }
        }
        break;
    case GAMMA51:
    case GAMMA52:
    case GAMMA53:
    case GAMMA54:
    {
        const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        for (INT idx = 0; idx < 4; ++idx)
        {
            UBOOL bPlus123[3] = { (0 == (idx & 1)), (0 == (idx & 2)), TRUE };
            if (byMissingDir < 2)
            {
                bPlus123[0] = TRUE;
                bPlus123[1] = (0 == (idx & 1));
                bPlus123[2] = (0 == (idx & 2));
            }

            INT dim123[3];
            BYTE byDimIndex = 0;
            for (INT byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
            {
                if (byCubeDir != static_cast<INT>(byMissingDir))
                {
                    dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                    byDimIndex++;
                }
            }

            INT dimordered[3];
            for (INT order = 0; order < 6; ++order)
            {
                switch (order)
                {
                case 1:
                    dimordered[0] = dim123[0];
                    dimordered[1] = dim123[2];
                    dimordered[2] = dim123[1];
                    break;
                case 2:
                    dimordered[0] = dim123[1];
                    dimordered[1] = dim123[0];
                    dimordered[2] = dim123[2];
                    break;
                case 3:
                    dimordered[0] = dim123[1];
                    dimordered[1] = dim123[2];
                    dimordered[2] = dim123[0];
                    break;
                case 4:
                    dimordered[0] = dim123[2];
                    dimordered[1] = dim123[0];
                    dimordered[2] = dim123[1];
                    break;
                case 5:
                    dimordered[0] = dim123[2];
                    dimordered[1] = dim123[1];
                    dimordered[2] = dim123[0];
                    break;
                default:
                    dimordered[0] = dim123[0];
                    dimordered[1] = dim123[1];
                    dimordered[2] = dim123[2];
                    break;
                }

                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 3, cudaMemcpyHostToDevice));

                _kernelDFermionKSForce_WithLink_Gamma51234 << <block, threads >> > (
                    (const deviceSU3*)pGaugeBuffer,
                    (deviceSU3*)pForce,
                    pRationalFields,
                    pRationalNumerator,
                    uiRationalDegree,
                    byFieldID,
                    byGaugeFieldID,
                    fCoeff * OneOver24,
                    byMissingDir,
                    devicePathBuffer
                    );
            }
        }
    }
    break;
    case GAMMA5:
        {
            for (INT idx = 0; idx < 8; ++idx)
            {
                const UBOOL bPlus1234[4] =
                {
                    (0 == (idx & 1)),
                    (0 == (idx & 2)),
                    (0 == (idx & 4)),
                    TRUE
                };

                const INT dim1234[4] =
                {
                    bPlus1234[0] ? 1 : -1,
                    bPlus1234[1] ? 2 : -2,
                    bPlus1234[2] ? 3 : -3,
                    4
                };

                INT dimordered[4];
                INT dim234[3];
                for (BYTE k = 0; k < 4; ++k)
                {
                    dimordered[0] = dim1234[k];
                    for (BYTE k2 = 0; k2 < 3; ++k2)
                    {
                        BYTE idx2 = k2 + 1 + k;
                        idx2 = idx2 > 3 ? (idx2 - 4) : idx2;
                        dim234[k2] = dim1234[idx2];
                    }

                    for (BYTE order2 = 0; order2 < 6; ++order2)
                    {
                        switch (order2)
                        {
                        case 1:
                            dimordered[1] = dim234[0];
                            dimordered[2] = dim234[2];
                            dimordered[3] = dim234[1];
                            break;
                        case 2:
                            dimordered[1] = dim234[1];
                            dimordered[2] = dim234[0];
                            dimordered[3] = dim234[2];
                            break;
                        case 3:
                            dimordered[1] = dim234[1];
                            dimordered[2] = dim234[2];
                            dimordered[3] = dim234[0];
                            break;
                        case 4:
                            dimordered[1] = dim234[2];
                            dimordered[2] = dim234[0];
                            dimordered[3] = dim234[1];
                            break;
                        case 5:
                            dimordered[1] = dim234[2];
                            dimordered[2] = dim234[1];
                            dimordered[3] = dim234[0];
                            break;
                        default:
                            dimordered[1] = dim234[0];
                            dimordered[2] = dim234[1];
                            dimordered[3] = dim234[2];
                            break;
                        }
                        //appGeneral(_T("dimordered=%d %d %d %d\n"), dimordered[0], dimordered[1], dimordered[2], dimordered[3]);
                        checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 4, cudaMemcpyHostToDevice));
                        _kernelDFermionKSForce_WithLink_Gamma5 << <block, threads >> > (
                            (const deviceSU3*)pGaugeBuffer,
                            (deviceSU3*)pForce,
                            pRationalFields,
                            pRationalNumerator,
                            uiRationalDegree,
                            byFieldID,
                            byGaugeFieldID,
                            fCoeff * OneOver192,
                            devicePathBuffer
                            );
                    }
                }
            }
        }
        break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}

#pragma region DOperator

#pragma region kernel



#pragma endregion


void CFieldFermionKSSU3Gamma::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionKSSU3::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    if (abs(m_fCoeffGamma1) > _CLG_FLT_EPSILON)
    {
        //bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,

        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA1, m_bEachSiteEta, bDagger, m_fCoeffGamma1, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma2) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA2, m_bEachSiteEta, bDagger, m_fCoeffGamma2, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma3) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA3, m_bEachSiteEta, bDagger, m_fCoeffGamma3, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma4) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA4, m_bEachSiteEta, bDagger, m_fCoeffGamma4, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma12) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA12, m_bEachSiteEta, bDagger, m_fCoeffSigma12, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma13) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA31, m_bEachSiteEta, bDagger, m_fCoeffSigma13, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma14) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA41, m_bEachSiteEta, bDagger, m_fCoeffSigma14, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma23) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA23, m_bEachSiteEta, bDagger, m_fCoeffSigma23, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma24) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA42, m_bEachSiteEta, bDagger, m_fCoeffSigma24, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma34) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, SIGMA43, m_bEachSiteEta, bDagger, m_fCoeffSigma34, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma51) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA51, m_bEachSiteEta, bDagger, m_fCoeffGamma51, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma52) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA52, m_bEachSiteEta, bDagger, m_fCoeffGamma52, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma53) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA53, m_bEachSiteEta, bDagger, m_fCoeffGamma53, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma54) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA54, m_bEachSiteEta, bDagger, m_fCoeffGamma54, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma5) > _CLG_FLT_EPSILON)
    {
        appApplyGammaKS(pTargetBuffer, pBuffer, pGaugeBuffer, GAMMA5, m_bEachSiteEta, bDagger, m_fCoeffGamma5, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
void CFieldFermionKSSU3Gamma::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    preparethread;
    CFieldFermionKSSU3::DerivateD0(pForce, pGaugeBuffer, byGaugeFieldId);
    if (abs(m_fCoeffGamma1) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce, 
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma1, 
            GAMMA1,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma2) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma2,
            GAMMA2,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma3) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma3,
            GAMMA3,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma4) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma4,
            GAMMA4,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma12) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma12,
            SIGMA12,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma13) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma13,
            SIGMA31,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma14) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma14,
            SIGMA41,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma23) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma23,
            SIGMA23,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma24) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma24,
            SIGMA42,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma34) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffSigma34,
            SIGMA43,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma51) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma51,
            GAMMA51,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma52) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma52,
            GAMMA52,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma53) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma53,
            GAMMA53,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma54) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma54,
            GAMMA54,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma5) > _CLG_FLT_EPSILON)
    {
        GammaKSForce(pForce,
            pGaugeBuffer,
            m_pRationalFieldPointers,
            m_pMDNumerator,
            m_rMD.m_uiDegree,
            m_fCoeffGamma5,
            GAMMA5,
            m_pDevicePathBuffer,
            m_byFieldId,
            byGaugeFieldId);
    }
}

#pragma endregion

#pragma region Other Kernel



#pragma endregion

CFieldFermionKSSU3Gamma::CFieldFermionKSSU3Gamma()
    : CFieldFermionKSSU3()
    , m_bImagine(TRUE)
    , m_fCoeffGamma1(F(0.0))
    , m_fCoeffGamma2(F(0.0))
    , m_fCoeffGamma3(F(0.0))
    , m_fCoeffGamma4(F(0.0))
    , m_fCoeffGamma5(F(0.0))
    , m_fCoeffGamma51(F(0.0))
    , m_fCoeffGamma52(F(0.0))
    , m_fCoeffGamma53(F(0.0))
    , m_fCoeffGamma54(F(0.0))
    , m_fCoeffSigma12(F(0.0))
    , m_fCoeffSigma13(F(0.0))
    , m_fCoeffSigma14(F(0.0))
    , m_fCoeffSigma23(F(0.0))
    , m_fCoeffSigma24(F(0.0))
    , m_fCoeffSigma34(F(0.0))
    , m_pDevicePathBuffer(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePathBuffer, sizeof(INT) * 4));
}

CFieldFermionKSSU3Gamma::~CFieldFermionKSSU3Gamma()
{
    checkCudaErrors(cudaFree(m_pDevicePathBuffer));
}

void CFieldFermionKSSU3Gamma::InitialOtherParameters(CParameters & params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);

    INT iValue = 1;
    if (params.FetchValueINT(_T("Imagine"), iValue))
    {
        m_bImagine = (0 != iValue);
        if (!m_bImagine)
        {
            appCrucial(_T("NOTE: Imagine is set to be FALSE, the M.C. has sign problem!!! Do this only in measure!!!\n"));
        }
    }

    Real fValue = F(0.0);
    if (params.FetchValueReal(_T("Gamma1"), fValue))
    {
        m_fCoeffGamma1 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma2"), fValue))
    {
        m_fCoeffGamma2 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma3"), fValue))
    {
        m_fCoeffGamma3 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma4"), fValue))
    {
        m_fCoeffGamma4 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma5"), fValue))
    {
        m_fCoeffGamma5 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma51"), fValue))
    {
        m_fCoeffGamma51 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma52"), fValue))
    {
        m_fCoeffGamma52 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma53"), fValue))
    {
        m_fCoeffGamma53 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma54"), fValue))
    {
        m_fCoeffGamma54 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma12"), fValue))
    {
        m_fCoeffSigma12 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma13"), fValue))
    {
        m_fCoeffSigma13 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma14"), fValue))
    {
        m_fCoeffSigma14 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma23"), fValue))
    {
        m_fCoeffSigma23 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma24"), fValue))
    {
        m_fCoeffSigma24 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma34"), fValue))
    {
        m_fCoeffSigma34 = fValue;
    }
}

void CFieldFermionKSSU3Gamma::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
    CFieldFermionKSSU3Gamma* pOther = dynamic_cast<CFieldFermionKSSU3Gamma*>(U);

    pOther->m_fCoeffGamma1 = m_fCoeffGamma1;
    pOther->m_fCoeffGamma2 = m_fCoeffGamma2;
    pOther->m_fCoeffGamma3 = m_fCoeffGamma3;
    pOther->m_fCoeffGamma4 = m_fCoeffGamma4;
    pOther->m_fCoeffGamma5 = m_fCoeffGamma5;
    pOther->m_fCoeffGamma51 = m_fCoeffGamma51;
    pOther->m_fCoeffGamma52 = m_fCoeffGamma52;
    pOther->m_fCoeffGamma53 = m_fCoeffGamma53;
    pOther->m_fCoeffGamma54 = m_fCoeffGamma54;
    pOther->m_fCoeffSigma12 = m_fCoeffSigma12;
    pOther->m_fCoeffSigma13 = m_fCoeffSigma13;
    pOther->m_fCoeffSigma14 = m_fCoeffSigma14;
    pOther->m_fCoeffSigma23 = m_fCoeffSigma23;
    pOther->m_fCoeffSigma24 = m_fCoeffSigma24;
    pOther->m_fCoeffSigma34 = m_fCoeffSigma34;
}

void CFieldFermionKSSU3Gamma::ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    InitialField(EFIT_Zero);

    //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
    UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
        || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);

    appApplyGammaKS(
        m_pDeviceData,
        pPooled->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        eGamma,
        m_bEachSiteEta,
        FALSE,
        F(0.5),
        bImag ? EOCT_Complex : EOCT_None,
        F(1.0),
        bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,
        m_byFieldId,
        pGauge->m_byFieldId
    );

    pPooled->Return();
}

CCString CFieldFermionKSSU3Gamma::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKSSU3::GetInfos(tab);

    sRet = sRet + tab + _T("Gamma1 : ") + appToString(m_fCoeffGamma1) + _T("\n");
    sRet = sRet + tab + _T("Gamma2 : ") + appToString(m_fCoeffGamma2) + _T("\n");
    sRet = sRet + tab + _T("Gamma3 : ") + appToString(m_fCoeffGamma3) + _T("\n");
    sRet = sRet + tab + _T("Gamma4 : ") + appToString(m_fCoeffGamma4) + _T("\n");
    sRet = sRet + tab + _T("Gamma5 : ") + appToString(m_fCoeffGamma5) + _T("\n");

    sRet = sRet + tab + _T("Gamma51 : ") + appToString(m_fCoeffGamma51) + _T("\n");
    sRet = sRet + tab + _T("Gamma52 : ") + appToString(m_fCoeffGamma52) + _T("\n");
    sRet = sRet + tab + _T("Gamma53 : ") + appToString(m_fCoeffGamma53) + _T("\n");
    sRet = sRet + tab + _T("Gamma54 : ") + appToString(m_fCoeffGamma54) + _T("\n");

    sRet = sRet + tab + _T("Sigma12 : ") + appToString(m_fCoeffSigma12) + _T("\n");
    sRet = sRet + tab + _T("Sigma13 : ") + appToString(m_fCoeffSigma13) + _T("\n");
    sRet = sRet + tab + _T("Sigma14 : ") + appToString(m_fCoeffSigma14) + _T("\n");
    sRet = sRet + tab + _T("Sigma23 : ") + appToString(m_fCoeffSigma23) + _T("\n");
    sRet = sRet + tab + _T("Sigma24 : ") + appToString(m_fCoeffSigma24) + _T("\n");
    sRet = sRet + tab + _T("Sigma34 : ") + appToString(m_fCoeffSigma34) + _T("\n");

    return sRet;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================