//=============================================================================
// FILENAME : CFieldFermionKSTKernel.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [mm/dd/yy]
//  [05/12/2025 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKST.h"

__BEGIN_NAMESPACE

#pragma region Gamma for KS

#pragma region gamma kernels

/**
 * gamma_i ^dagger = gamma_i for i = 1,2,3,4,5, so 'bDDagger' is used for only coefficient
 *
 * 1/2a eta_mu(x) sum_{mu=+-} psibar(x) U_mu psi(x+mu)
 * the 1/2a is absorbed (be sure to add when measuring)
 */
    template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
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

    const deviceGauge& x_Gauge_element = pGauge[linkIndex];
    deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        _mul(result, F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    //result.MulReal(eta_mu);

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
* This function applys i Gamma_i on the field
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEta1234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
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

    const deviceGauge& x_Gauge_element = pGauge[linkIndex];
    deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, x_m_mu_Gauge.m_byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        _mul(result, F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234EvenOddT(
    deviceVector* pMe,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernalInt4EO;

    const Real eta_mu = ((pEtaTable[uiSiteIndex] >> byDir) & 1) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    const deviceGauge& x_Gauge_element = pGauge[linkIndex];
    deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pMe[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        _mul(result, F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pMe[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pMe[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSForceGamma1234T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
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
        const deviceVector* phi_i = pFermionPointers[uiR];
        const deviceVector* phi_id = pFermionPointers[uiR + uiRational];

        //deviceVector toContract = _mulVec(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_i[x_p_mu_Fermion.m_uiSiteIndex], phi_id[uiSiteIndex]);

        //toContract = _mulVec(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
        _sub(thisTerm, _makeContract<deviceGauge, deviceVector>(phi_id[x_p_mu_Fermion.m_uiSiteIndex], phi_i[uiSiteIndex]));

        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(thisTerm, _make_cuComplex(F(0.0), eta_mu * pNumerators[uiR] * F(-1.0) * fCoeff));
        }
        else
        {
            _mul(thisTerm, _make_cuComplex(F(0.0), eta_mu * pNumerators[uiR] * fCoeff));
        }
        
        _add(pForce[linkIndex], thisTerm);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSForceGamma1234EvenOddT(
    deviceGauge* pForce,
    BYTE idir,
    Real fCoeff,
    BYTE byFieldId)
{
    intokernal;

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
    deviceGauge f0 = pForce[linkIndex];
    _mul(f0, _make_cuComplex(F(0.0), fCoeff));
    _sub(pForce[linkIndex], f0);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaSigmaIJT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    BYTE byEtaShift,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    SCHAR byDir1,
    SCHAR byDir2,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlus12[2] = { (0 != (idx & 1)), (0 != (idx & 2)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim12[2] = {
            bPlus12[0] ? static_cast<SCHAR>(byDir1 + 1) : static_cast<SCHAR>(-byDir1 - 1),
            bPlus12[1] ? static_cast<SCHAR>(byDir2 + 1) : static_cast<SCHAR>(-byDir2 - 1)
        };
        sOffset.m_byData4[byDir1] = sOffset.m_byData4[byDir1] + (bPlus12[0] ? 1 : -1);
        sOffset.m_byData4[byDir2] = sOffset.m_byData4[byDir2] + (bPlus12[1] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_devicePlaneDiagonalT(pGauge, sSite4, byGaugeFieldId, dim12[0], dim12[1]), pOther[sTargetBigIndex.m_uiSiteIndex]);

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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.5) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.5) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 * gamma5i corresponds to diagonal links of cubic in the other 3 dimensions.
 * For each diagonal link, there are 6 different ways to add the gauge field
 * We simply use the average of all 6 links
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma51234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    UBOOL bDDagger,
    BYTE byMissingDir,
    SCHAR byEtaShift,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlus123[3] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim123[3] = { 0, 0, 0 };
        BYTE byDimIndex = 0;
        for (SCHAR byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
        {
            if (byCubeDir != static_cast<SCHAR>(byMissingDir))
            {
                sOffset.m_byData4[byCubeDir] = sOffset.m_byData4[byCubeDir] + (bPlus123[byDimIndex] ? 1 : -1);
                dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                byDimIndex++;
            }
        }

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceCubicDiagonalT(pGauge, sSite4, byGaugeFieldId, dim123[0], dim123[1], dim123[2]), pOther[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.25) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.25) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 *
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma5T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
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

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 16; ++idx)
    {
        const UBOOL bPlus1234[4] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)), (0 != (idx & 8)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim1234[4] =
        {
            bPlus1234[0] ? static_cast<SCHAR>(1) : static_cast<SCHAR>(-1),
            bPlus1234[1] ? static_cast<SCHAR>(2) : static_cast<SCHAR>(-2),
            bPlus1234[2] ? static_cast<SCHAR>(3) : static_cast<SCHAR>(-3),
            bPlus1234[3] ? static_cast<SCHAR>(4) : static_cast<SCHAR>(-4)
        };
        sOffset.m_byData4[0] = sOffset.m_byData4[0] + (bPlus1234[0] ? 1 : -1);
        sOffset.m_byData4[1] = sOffset.m_byData4[1] + (bPlus1234[1] ? 1 : -1);
        sOffset.m_byData4[2] = sOffset.m_byData4[2] + (bPlus1234[2] ? 1 : -1);
        sOffset.m_byData4[3] = sOffset.m_byData4[3] + (bPlus1234[3] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceHyperCubicDiagonalT(pGauge, sSite4, byGaugeFieldId, dim1234[0], dim1234[1], dim1234[2], dim1234[3]), pOther[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, _make_cuComplex(F(0.0), -F(0.125) * fGammCoefficient));
    }
    else
    {
        _mul(result, _make_cuComplex(F(0.0), F(0.125) * fGammCoefficient));
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}


template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234EMEvenOddT(
    deviceVector* pMe,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCharge,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernalInt4EO;

    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1)) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    deviceGauge x_Gauge_element = pGauge[linkIndex];
    const Real forwardPhase = pU1[linkIndex] * fCharge;
    _mul(x_Gauge_element, _make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
    const UINT backwardLink = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir);
    deviceGauge x_m_mu_Gauge_element = pGauge[backwardLink];
    const Real backPhase = pU1[backwardLink] * fCharge;
    _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(backPhase), _sin(backPhase)));

    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pMe[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        _mul(result, F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pMe[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pMe[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    //result.MulReal(eta_mu);

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

#pragma region force

/**
 * similar as _kernelDFermionKSForce_WithLink.
 * but _kernelDFermionKSForce_WithLink use gamma mu as the gamma matrix
 * we use sigma12
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_SigmaIJT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE bDir1,
    BYTE bDir2,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[2];
    SCHAR pathRight[2];
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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
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
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma51234T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE byMissingDir,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[3];
    SCHAR pathRight[3];
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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

//similar as above
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma5T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[4];
    SCHAR pathRight[4];

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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    //_mul(res, fGammCoefficient * pNumerators[rfieldId]);
                    _mul(res, _make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    //_mul(res, fGammCoefficient * pNumerators[rfieldId]);
                    _mul(res, _make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));

                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_SigmaIJEMT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    Real fCharge,
    BYTE bDir1,
    BYTE bDir2,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[2];
    SCHAR pathRight[2];
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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
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
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma51234EMT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    Real fCharge,
    BYTE byMissingDir,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[3];
    SCHAR pathRight[3];
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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

//similar as above
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma5EMT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    Real fCharge,
    const SCHAR* __restrict__ path)
{
    intokernalInt4;
    SCHAR pathLeft[4];
    SCHAR pathRight[4];

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
            if (bHasLeft)
            {
                const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
                    //_mul(res, fGammCoefficient * pNumerators[rfieldId]);
                    _mul(res, _make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));

                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }

            if (bHasRight)
            {
                const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, RLength, byGaugeFieldId, pathRight);

                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);

                for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
                {
                    const deviceVector* phi_i = pFermionPointers[rfieldId];
                    const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                    //=================================
                    // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                    deviceVector phi1 = (LLength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                    deviceVector phi2 = (RLength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];

                    deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                    //This Add is required by partial(D^+D)
                    phi1 = (RLength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                    phi2 = (LLength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                    _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                    //_mul(res, fGammCoefficient * pNumerators[rfieldId]);
                    _mul(res, _make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));

                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

#pragma endregion

#pragma endregion


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::appApplyGammaKS(
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
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceVector* pSource = (const deviceVector*)pBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;
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
            _LAUNCH_KERNEL(_kernelKSApplyGammaEta1234T TMPARG(deviceVector, deviceGauge), block, threads,
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
            _LAUNCH_KERNEL(_kernelKSApplyGamma1234T TMPARG(deviceVector, deviceGauge), block, threads,
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
            _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
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
        SCHAR etaShift = -1;
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

        _LAUNCH_KERNEL(_kernelKSApplyGamma51234T TMPARG(deviceVector, deviceGauge), block, threads,
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
        _LAUNCH_KERNEL(_kernelKSApplyGamma5T TMPARG(deviceVector, deviceGauge), block, threads,
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

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::appApplyGammaKSEvenOdd(
    void* pTargetBuffer,
    const void* pGaugeBuffer,
    EGammaMatrix eGamma,
    UBOOL bShiftCenter,
    UBOOL bEven,
    UBOOL bDagger,
    Real fGammaCoeff,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;
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
            appCrucial(_T("not supported yet\n"));
        }
        else
        {
            _LAUNCH_KERNEL(_kernelKSApplyGamma1234EvenOddT TMPARG(deviceVector, deviceGauge), block, threads,
                pTarget,
                pGauge,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEven,
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
    case SIGMA31:
    case SIGMA41:
    case SIGMA23:
    case SIGMA42:
    case SIGMA43:
        appCrucial(_T("Sigma ij should not support even odd!\n"));
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

        //const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        //SCHAR etaShift = -1;
        //if (bShiftCenter)
        //{
        //    if (byMissingDir < 2)
        //    {
        //        etaShift = 0;
        //    }
        //    else
        //    {
        //        etaShift = 2;
        //    }
        //}

        appCrucial(_T("not supported yet\n"));
    }
    break;
    case GAMMA5:
        appCrucial(_T("not supported yet\n"));
        break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::GammaKSForce(
    void* pForce,
    const void* pGaugeBuffer,
    const deviceVector* const* pRationalFields,
    const Real* pRationalNumerator,
    UINT uiRationalDegree,
    Real fCoeff,
    EGammaMatrix eGamma,
    SCHAR* devicePathBuffer,
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
        _LAUNCH_KERNEL(_kernelFermionKSForceGamma1234T TMPARG(deviceVector, deviceGauge), block, threads,
            (const deviceGauge*)pGaugeBuffer,
            (deviceGauge*)pForce,
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
        SCHAR byDirs[2] = { 0, 1 };
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

            SCHAR dimordered[2];
            for (INT order = 0; order < 2; ++order)
            {
                dimordered[0] = byDirs[0] + 1;
                dimordered[1] = bPlus2 ? static_cast<SCHAR>(byDirs[1] + 1) : static_cast<SCHAR>(-byDirs[1] - 1);
                if (1 == order)
                {
                    dimordered[0] = bPlus2 ? static_cast<SCHAR>(byDirs[1] + 1) : static_cast<SCHAR>(-byDirs[1] - 1);
                    dimordered[1] = byDirs[0] + 1;
                }
                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 2, cudaMemcpyHostToDevice));

                _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_SigmaIJT TMPARG(deviceVector, deviceGauge), block, threads,
                    (const deviceGauge*)pGaugeBuffer,
                    (deviceGauge*)pForce,
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

            SCHAR dim123[3];
            BYTE byDimIndex = 0;
            for (SCHAR byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
            {
                if (byCubeDir != static_cast<INT>(byMissingDir))
                {
                    dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                    byDimIndex++;
                }
            }

            SCHAR dimordered[3];
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

                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 3, cudaMemcpyHostToDevice));

                _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_Gamma51234T TMPARG(deviceVector, deviceGauge), block, threads,
                    (const deviceGauge*)pGaugeBuffer,
                    (deviceGauge*)pForce,
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

            const SCHAR dim1234[4] =
            {
                bPlus1234[0] ? static_cast<SCHAR>(1) : static_cast<SCHAR>(-1),
                bPlus1234[1] ? static_cast<SCHAR>(2) : static_cast<SCHAR>(-2),
                bPlus1234[2] ? static_cast<SCHAR>(3) : static_cast<SCHAR>(-3),
                static_cast<SCHAR>(4)
            };

            SCHAR dimordered[4];
            SCHAR dim234[3];
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
                    checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 4, cudaMemcpyHostToDevice));
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_Gamma5T TMPARG(deviceVector, deviceGauge), block, threads,
                        (const deviceGauge*)pGaugeBuffer,
                        (deviceGauge*)pForce,
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

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::GammaKSForceEvenOdd(
    void* pForce,
    Real fCoeff,
    EGammaMatrix eGamma,
    BYTE byFieldID)
{
    deviceGauge* forcesu3 = (deviceGauge*)(pForce);
    preparethread;

    switch (eGamma)
    {
    case GAMMA1:
    case GAMMA2:
    case GAMMA3:
    case GAMMA4:
    {
        BYTE iDir = static_cast<BYTE>(eGamma) - 1;
        _LAUNCH_KERNEL(_kernelFermionKSForceGamma1234EvenOddT<deviceGauge>, block, threads, forcesu3, iDir, fCoeff, byFieldID);
    }
    break;
    case SIGMA12:
    case SIGMA31:
    case SIGMA41:
    case SIGMA23:
    case SIGMA42:
    case SIGMA43:
        appCrucial(_T("Sigma ij should not support even odd!\n"));
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

        //const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        //SCHAR etaShift = -1;
        //if (bShiftCenter)
        //{
        //    if (byMissingDir < 2)
        //    {
        //        etaShift = 0;
        //    }
        //    else
        //    {
        //        etaShift = 2;
        //    }
        //}

        appCrucial(_T("not supported yet\n"));
    }
    break;
    case GAMMA5:
        appCrucial(_T("not supported yet\n"));
        break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}



#pragma region gamma EM

#pragma region kernel

/**
 * gamma_i ^dagger = gamma_i for i = 1,2,3,4,5, so 'bDDagger' is used for only coefficient
 *
 * 1/2a eta_mu(x) sum_{mu=+-} psibar(x) U_mu psi(x+mu)
 * the 1/2a is absorbed (be sure to add when measuring)
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234EMT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fCharge,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernal;

    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1)) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    deviceGauge x_Gauge_element = pGauge[linkIndex];
    const Real forwardPhase = pU1[linkIndex] * fCharge;
    _mul(x_Gauge_element, _make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
    const UINT backwardLink = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir);
    deviceGauge x_m_mu_Gauge_element = pGauge[backwardLink];
    const Real backPhase = pU1[backwardLink] * fCharge;
    _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(backPhase), _sin(backPhase)));

    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        _mul(result, F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    //result.MulReal(eta_mu);

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEta1234EMT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fCharge,
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

    BYTE eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1));
    BYTE eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir) & 1));

    deviceGauge x_Gauge_element = pGauge[linkIndex];
    const Real forwardPhase = pU1[linkIndex] * fCharge;
    _mul(x_Gauge_element, _make_cuComplex(_cos(forwardPhase), _sin(forwardPhase)));
    const UINT backwardLink = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir);
    deviceGauge x_m_mu_Gauge_element = pGauge[backwardLink];
    const Real backPhase = pU1[backwardLink] * fCharge;
    _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(backPhase), _sin(backPhase)));

    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        _mul(result, F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaSigmaIJEMT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fCharge,
    Real fGammCoefficient,
    BYTE byEtaShift,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    SCHAR byDir1,
    SCHAR byDir2,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlus12[2] = { (0 != (idx & 1)), (0 != (idx & 2)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim12[2] = {
            bPlus12[0] ? static_cast<SCHAR>(byDir1 + 1) : static_cast<SCHAR>(-byDir1 - 1),
            bPlus12[1] ? static_cast<SCHAR>(byDir2 + 1) : static_cast<SCHAR>(-byDir2 - 1)
        };
        sOffset.m_byData4[byDir1] = sOffset.m_byData4[byDir1] + (bPlus12[0] ? 1 : -1);
        sOffset.m_byData4[byDir2] = sOffset.m_byData4[byDir2] + (bPlus12[1] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_devicePlaneDiagonalEMT(pGauge, pU1, fCharge, sSite4, byGaugeFieldId, dim12[0], dim12[1]), pOther[sTargetBigIndex.m_uiSiteIndex]);

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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.5) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.5) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 * gamma5i corresponds to diagonal links of cubic in the other 3 dimensions.
 * For each diagonal link, there are 6 different ways to add the gauge field
 * We simply use the average of all 6 links
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma51234EMT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    UBOOL bDDagger,
    BYTE byMissingDir,
    SCHAR byEtaShift,
    Real fCharge,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlus123[3] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim123[3] = { 0, 0, 0 };
        BYTE byDimIndex = 0;
        for (SCHAR byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
        {
            if (byCubeDir != static_cast<SCHAR>(byMissingDir))
            {
                sOffset.m_byData4[byCubeDir] = sOffset.m_byData4[byCubeDir] + (bPlus123[byDimIndex] ? 1 : -1);
                dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                byDimIndex++;
            }
        }

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceCubicDiagonalEMT(pGauge, pU1, fCharge, sSite4, byGaugeFieldId, dim123[0], dim123[1], dim123[2]), pOther[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.25) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.25) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 *
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma5EMT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    UBOOL bDDagger,
    UBOOL bEtaShift,
    Real fCharge,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

    #pragma unroll
    for (UINT idx = 0; idx < 16; ++idx)
    {
        const UBOOL bPlus1234[4] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)), (0 != (idx & 8)) };

        SSmallInt4 sOffset = sSite4;
        SCHAR dim1234[4] =
        {
            bPlus1234[0] ? static_cast<SCHAR>(1) : static_cast<SCHAR>(-1),
            bPlus1234[1] ? static_cast<SCHAR>(2) : static_cast<SCHAR>(-2),
            bPlus1234[2] ? static_cast<SCHAR>(3) : static_cast<SCHAR>(-3),
            bPlus1234[3] ? static_cast<SCHAR>(4) : static_cast<SCHAR>(-4)
        };
        sOffset.m_byData4[0] = sOffset.m_byData4[0] + (bPlus1234[0] ? 1 : -1);
        sOffset.m_byData4[1] = sOffset.m_byData4[1] + (bPlus1234[1] ? 1 : -1);
        sOffset.m_byData4[2] = sOffset.m_byData4[2] + (bPlus1234[2] ? 1 : -1);
        sOffset.m_byData4[3] = sOffset.m_byData4[3] + (bPlus1234[3] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceHyperCubicDiagonalEMT(pGauge, pU1, fCharge, sSite4, byGaugeFieldId, dim1234[0], dim1234[1], dim1234[2], dim1234[3]), pOther[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, _make_cuComplex(F(0.0), -F(0.125) * fGammCoefficient));
    }
    else
    {
        _mul(result, _make_cuComplex(F(0.0), F(0.125) * fGammCoefficient));
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::appApplyGammaKSEM(
    void* pTargetBuffer,
    const void* pBuffer,
    const void* pGaugeBuffer,
    const void* pEMFieldBuffer,
    Real fCharge,
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
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceVector* pSource = (const deviceVector*)pBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;
    const Real* pU1 = (const Real*)pEMFieldBuffer;

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
            _LAUNCH_KERNEL(_kernelKSApplyGammaEta1234EMT TMPARG(deviceVector, deviceGauge), block, threads,
                pTarget,
                pSource,
                pGauge,
                pU1,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fCharge,
                fGammaCoeff,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                static_cast<BYTE>(iDir));
        }
        else
        {
            _LAUNCH_KERNEL(_kernelKSApplyGamma1234EMT TMPARG(deviceVector, deviceGauge), block, threads,
                pTarget,
                pSource,
                pGauge,
                pU1,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fCharge,
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
            _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
                pTarget,
                pSource,
                pGauge,
                pU1,
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fCharge,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fCharge,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fCharge,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fCharge,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fCharge,
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
        _LAUNCH_KERNEL(_kernelKSApplyGammaSigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fCharge,
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
        SCHAR etaShift = -1;
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

        _LAUNCH_KERNEL(_kernelKSApplyGamma51234EMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            bDagger,
            byMissingDir,
            etaShift,
            fCharge,
            fGammaCoeff,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            byFieldID,
            byGaugeFieldID);
    }
    break;
    case GAMMA5:
        _LAUNCH_KERNEL(_kernelKSApplyGamma5EMT TMPARG(deviceVector, deviceGauge), block, threads,
            pTarget,
            pSource,
            pGauge,
            pU1,
            bDagger,
            bShiftCenter,
            fCharge,
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


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::appApplyGammaKSEMEvenOdd(
    void* pTargetBuffer,
    const void* pGaugeBuffer,
    const void* pEMFieldBuffer,
    Real fCharge,
    EGammaMatrix eGamma,
    UBOOL bShiftCenter,
    UBOOL bEven,
    UBOOL bDagger,
    Real fGammaCoeff,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const Real* pU1 = (const Real*)pEMFieldBuffer;

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
            appGeneral(_T("not implemented yet.\n"));
            //_LAUNCH_KERNEL(_kernelKSApplyGammaEta1234EM, block, threads, 
            //    pTarget,
            //    pSource,
            //    pGauge,
            //    pU1,
            //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
            //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            //    appGetLattice()->m_pIndexCache->m_pEtaMu,
            //    bDagger,
            //    fCharge,
            //    fGammaCoeff,
            //    eOCT,
            //    fRealCoeff,
            //    cCmpCoeff,
            //    static_cast<BYTE>(iDir));
        }
        else
        {
            _LAUNCH_KERNEL(_kernelKSApplyGamma1234EMEvenOddT TMPARG(deviceSU3Vector, deviceSU3), block, threads,
                pTarget,
                pGauge,
                pU1,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEven,
                bDagger,
                fCharge,
                fGammaCoeff,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                static_cast<BYTE>(iDir));
        }
    }
    break;
    case SIGMA12:
    case SIGMA31:
    case SIGMA41:
    case SIGMA23:
    case SIGMA42:
    case SIGMA43:
        appCrucial(_T("Sigma ij should not support even odd!\n"));
        break;
        //case GAMMA51:
        //case GAMMA52:
        //case GAMMA53:
        //case GAMMA54:
        //{
            //eta shift is:
            //x->y
            //y->x
            //z->t
            //t->z

        //    const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        //    SCHAR etaShift = -1;
        //    if (bShiftCenter)
        //    {
        //        if (byMissingDir < 2)
        //        {
        //            etaShift = 0;
        //        }
        //        else
        //        {
        //            etaShift = 2;
        //        }
        //    }

        //    _LAUNCH_KERNEL(_kernelKSApplyGamma51234EM, block, threads, 
        //        pTarget,
        //        pSource,
        //        pGauge,
        //        pU1,
        //        bDagger,
        //        byMissingDir,
        //        etaShift,
        //        fCharge,
        //        fGammaCoeff,
        //        eOCT,
        //        fRealCoeff,
        //        cCmpCoeff,
        //        byFieldID,
        //        byGaugeFieldID);
        //}
        //break;
        //case GAMMA5:
        //    _LAUNCH_KERNEL(_kernelKSApplyGamma5EM, block, threads, 
        //        pTarget,
        //        pSource,
        //        pGauge,
        //        pU1,
        //        bDagger,
        //        bShiftCenter,
        //        fCharge,
        //        fGammaCoeff,
        //        eOCT,
        //        fRealCoeff,
        //        cCmpCoeff,
        //        byFieldID,
        //        byGaugeFieldID);
        //    break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelGamma<deviceVector, deviceGauge, vectorN>::GammaKSForceEM(
    void* pForce,
    const void* pGaugeBuffer,
    const void* pEMFieldBuffer,
    Real fCharge,
    const deviceVector* const* pRationalFields,
    const Real* pRationalNumerator,
    UINT uiRationalDegree,
    Real fCoeff,
    EGammaMatrix eGamma,
    SCHAR* devicePathBuffer,
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
        _LAUNCH_KERNEL(_kernelFermionKSForceGamma1234T TMPARG(deviceVector, deviceGauge), block, threads,
            (const deviceGauge*)pGaugeBuffer,
            //(const Real*)pEMFieldBuffer,
            (deviceGauge*)pForce,
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pRationalFields,
            pRationalNumerator,
            uiRationalDegree,
            byDir,
            fCoeff,
            //fCharge,
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
        SCHAR byDirs[2] = { 0, 1 };
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

            SCHAR dimordered[2];
            for (INT order = 0; order < 2; ++order)
            {
                dimordered[0] = byDirs[0] + 1;
                dimordered[1] = bPlus2 ? static_cast<SCHAR>(byDirs[1] + 1) : static_cast<SCHAR>(-byDirs[1] - 1);
                if (1 == order)
                {
                    dimordered[0] = bPlus2 ? static_cast<SCHAR>(byDirs[1] + 1) : static_cast<SCHAR>(-byDirs[1] - 1);
                    dimordered[1] = byDirs[0] + 1;
                }
                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 2, cudaMemcpyHostToDevice));

                _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_SigmaIJEMT TMPARG(deviceVector, deviceGauge), block, threads,
                    (const deviceGauge*)pGaugeBuffer,
                    (const Real*)pEMFieldBuffer,
                    (deviceGauge*)pForce,
                    pRationalFields,
                    pRationalNumerator,
                    uiRationalDegree,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    byFieldID,
                    byGaugeFieldID,
                    fCoeff * F(0.25),
                    fCharge,
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

            SCHAR dim123[3];
            BYTE byDimIndex = 0;
            for (SCHAR byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
            {
                if (byCubeDir != static_cast<SCHAR>(byMissingDir))
                {
                    dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                    byDimIndex++;
                }
            }

            SCHAR dimordered[3];
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

                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 3, cudaMemcpyHostToDevice));

                _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_Gamma51234EMT TMPARG(deviceVector, deviceGauge), block, threads,
                    (const deviceGauge*)pGaugeBuffer,
                    (const Real*)pEMFieldBuffer,
                    (deviceGauge*)pForce,
                    pRationalFields,
                    pRationalNumerator,
                    uiRationalDegree,
                    byFieldID,
                    byGaugeFieldID,
                    fCoeff * OneOver24,
                    fCharge,
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

            const SCHAR dim1234[4] =
            {
                static_cast<SCHAR>(bPlus1234[0] ? 1 : -1),
                static_cast<SCHAR>(bPlus1234[1] ? 2 : -2),
                static_cast<SCHAR>(bPlus1234[2] ? 3 : -3),
                4
            };

            SCHAR dimordered[4];
            SCHAR dim234[3];
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
                    checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(SCHAR) * 4, cudaMemcpyHostToDevice));
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLink_Gamma5EMT TMPARG(deviceVector, deviceGauge), block, threads,
                        (const deviceGauge*)pGaugeBuffer,
                        (const Real*)pEMFieldBuffer,
                        (deviceGauge*)pForce,
                        pRationalFields,
                        pRationalNumerator,
                        uiRationalDegree,
                        byFieldID,
                        byGaugeFieldID,
                        fCoeff * OneOver192,
                        fCharge,
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

#pragma endregion

#pragma endregion

template class CFieldFermionKSTKernelGamma<CLGComplex, CLGComplex, 1>;
template class CFieldFermionKSTKernelGamma<deviceSU2Vector, deviceSU2, 2>;
template class CFieldFermionKSTKernelGamma<deviceSU3Vector, deviceSU3, 3>;

template class CFieldFermionKSTKernelGamma<deviceSU4Vector, deviceSU4, 4>;
//template class CFieldFermionKSTKernelGamma<deviceSU5Vector, deviceSU5, 5>;
//template class CFieldFermionKSTKernelGamma<deviceSU6Vector, deviceSU6, 6>;
//template class CFieldFermionKSTKernelGamma<deviceSU7Vector, deviceSU7, 7>;
//template class CFieldFermionKSTKernelGamma<deviceSU8Vector, deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================