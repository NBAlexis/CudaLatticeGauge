//=============================================================================
// FILENAME : CFieldFermionKSSU3RigidAcc.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [12/27/2023 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CFieldFermionKSSU3Gamma.h"
#include "CFieldFermionKSSU3RigidAcc.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3RigidAcc)

#pragma region DOperator

#pragma region kernel

/**
* copy of _kernelDFermionKS
* but change coefficients for partial _i to (1+gz), change m to (1+gz)m
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSRigidAcc(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    Real f2am,
    Real fG,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    const Real f1pgz = F(1.0) + (sSite4.z - _DC_Centerz) * fG;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        const SSmallInt4 x_p_mu_site = __deviceSiteIndexToInt4(x_p_mu_Fermion.m_uiSiteIndex);
        const SSmallInt4 x_m_mu_site = __deviceSiteIndexToInt4(x_m_mu_Fermion.m_uiSiteIndex);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceSU3& x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (3 != idir)
        {
            u_phi_x_p_m.MulReal(F(1.0) + F(0.5) * fG * (sSite4.z + x_p_mu_site.z - 2 * _DC_Centerz));
        }
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (3 != idir)
        {
            u_dagger_phi_x_m_m.MulReal(F(1.0) + F(0.5) * fG * (sSite4.z + x_m_mu_site.z - 2 * _DC_Centerz));
        }
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        u_phi_x_p_m.MulReal(eta_mu);
        result.Add(u_phi_x_p_m);
    }

    pResultData[uiSiteIndex].MulReal(f2am * f1pgz);
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


__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceRigidAcc(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    Real fG,
    BYTE byFieldId)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SSmallInt4 x_p_mu_site = __deviceSiteIndexToInt4(x_p_mu_Fermion.m_uiSiteIndex);
        const Real f1pgz = F(1.0) + F(0.5) * fG * (sSite4.z + x_p_mu_site.z - 2 * _DC_Centerz);

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[uiR];
            const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];

            deviceSU3Vector toContract = pGauge[linkIndex].MulVector(phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
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
            if (3 != idir)
            {
                thisTerm.MulReal(f1pgz);
            }
            pForce[linkIndex].Sub(thisTerm);
        }
    }

}

#pragma endregion


void CFieldFermionKSSU3RigidAcc::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKSRigidAcc << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        CCommonData::m_fG,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    if (m_bUseImaginaryGamma3)
    {
        CFieldFermionKSSU3Gamma::appApplyGammaKS(pTarget, pSource, pGauge,
            GAMMA3, FALSE, bDagger, F(0.5) * CCommonData::m_fG,
            eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    else
    {
        CLGComplex toMul = _make_cuComplex(F(1.0), F(0.0));

        switch (eOCT)
        {
        case EOCT_None:
            toMul = _make_cuComplex(F(0.0), F(-0.5) * CCommonData::m_fG);
            break;
        case EOCT_Minus:
            toMul = _make_cuComplex(F(0.0), F(0.5) * CCommonData::m_fG);
            break;
        case EOCT_Real:
            toMul = _make_cuComplex(F(0.0), F(-0.5) * CCommonData::m_fG * fRealCoeff);
            break;
        case EOCT_Complex:
            toMul = _cuCmulf(cCmpCoeff, _make_cuComplex(F(0.0), F(-0.5) * CCommonData::m_fG));
            break;
        default:
            break;
        }
        CFieldFermionKSSU3Gamma::appApplyGammaKS(pTarget, pSource, pGauge,
            GAMMA3, FALSE, bDagger, F(1.0),
            EOCT_Complex, F(1.0), toMul, m_byFieldId, byGaugeFieldId);
    }
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
void CFieldFermionKSSU3RigidAcc::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    preparethread;
    _kernelDFermionKSForceRigidAcc << <block, threads >> > (
        (const deviceSU3*)pGaugeBuffer,
        (deviceSU3*)pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        CCommonData::m_fG,
        m_byFieldId);

    CFieldFermionKSSU3Gamma::GammaKSForce(pForce,
        pGaugeBuffer,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        F(0.5) * CCommonData::m_fG,
        GAMMA3,
        m_pDevicePathBuffer,
        m_byFieldId,
        byGaugeFieldId);
}

#pragma endregion

CFieldFermionKSSU3RigidAcc::CFieldFermionKSSU3RigidAcc() 
    : CFieldFermionKSSU3()
    , m_bUseImaginaryGamma3(TRUE)
    , m_pDevicePathBuffer(NULL)
{
    m_bDiagonalMass = TRUE;

    checkCudaErrors(cudaMalloc((void**)&m_pDevicePathBuffer, sizeof(INT) * 4));
}

CFieldFermionKSSU3RigidAcc::~CFieldFermionKSSU3RigidAcc()
{
    checkCudaErrors(cudaFree(m_pDevicePathBuffer));
}

void CFieldFermionKSSU3RigidAcc::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

void CFieldFermionKSSU3RigidAcc::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);

    INT iImGamma3 = 1;
    if (params.FetchValueINT(_T("ImaginaryGamma3"), iImGamma3))
    {
        m_bUseImaginaryGamma3 = (0 != iImGamma3);
    }
}

CCString CFieldFermionKSSU3RigidAcc::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKSSU3::GetInfos(tab);
    sRet = sRet + tab + _T("Imaginary Gamma3 : ") + appToString(m_bUseImaginaryGamma3) + _T("\n");
    sRet = sRet + tab + _T("G : ") + appToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================