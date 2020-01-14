//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

#define _MATRIX_BOUND 2

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKS)

#pragma region DOperator

#pragma region kernel

/**
* Dw phi(x) = phi(x) - kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)
* U act on su3
* gamma act on spinor
*
* If bDagger, it is gamma5, D, gamma5
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    const Real gamma5 = (1 == ((pEtaTable[uiSiteIndex] >> 4) & 1)) ? F(-1.0) : F(1.0);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        pResultData[uiSiteIndex].MulReal(gamma5);
    }

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        Real gammaMu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        x_m_mu_Gauge_element.Dagger();

        deviceSU3Vector x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        deviceSU3Vector x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];

        if (bDDagger)
        {
            x_p_mu_Fermion_element.MulReal(gamma5);
            x_m_mu_Fermion_element.MulReal(gamma5);
        }

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(x_p_mu_Fermion_element);
        u_phi_x_p_m.MulReal(gammaMu);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(u_phi_x_p_m);
        }
        else
        {
            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(u_phi_x_p_m);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(x_m_mu_Fermion_element);
        u_dagger_phi_x_m_m.MulReal(gammaMu);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);

    if (bDDagger)
    {
        pResultData[uiSiteIndex].MulReal(gamma5);
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

/**
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDKS(
    const deviceSU3Vector* __restrict__ pInverseD,
    const deviceSU3Vector* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3* pForce,
    Real fKai,
    BYTE byFieldId)
{
    intokernaldir;

    const deviceSU3Vector x_Left(pInverseDDdagger[uiSiteIndex]);
    const deviceSU3Vector x_Right(pInverseD[uiSiteIndex]);

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        Real gammaMu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        const deviceSU3Vector x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceSU3Vector x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

        deviceSU3 x_Gauge_element = pGauge[linkIndex];

        deviceSU3Vector right1(x_p_mu_Right);
        right1.Sub(right1.MulRealC(gammaMu));
        deviceSU3 mid = deviceSU3::makeSU3ContractV(x_Left, right1);

        deviceSU3Vector right2(x_Right);
        right2.Add(right2.MulRealC(gammaMu));
        mid.Add(deviceSU3::makeSU3ContractV(right2, x_p_mu_Left));

        deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
        forceOfThisLink.Ta();
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            forceOfThisLink.MulReal(-fKai);
        }
        else
        {
            forceOfThisLink.MulReal(fKai);
        }

        pForce[linkIndex].Add(forceOfThisLink);
    }
}

#pragma endregion

void CFieldFermionKS::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKS << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_fKai,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

}

void CFieldFermionKS::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceSU3Vector* pDphiBuffer = (deviceSU3Vector*)pDphi;
    const deviceSU3Vector* pDDphiBuffer = (deviceSU3Vector*)pDDphi;

    preparethread;
    _kernelDKS<< <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pForceSU3,
        m_fKai, m_byFieldId);
}

#pragma endregion

#pragma region Kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPrintFermionKS(const deviceSU3Vector* __restrict__ pData)
{
    intokernal;

    printf("%d=(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi))\n",
        uiSiteIndex,
        pData[uiSiteIndex].m_ve[0].x, pData[uiSiteIndex].m_ve[0].y,
        pData[uiSiteIndex].m_ve[1].x, pData[uiSiteIndex].m_ve[1].y,
        pData[uiSiteIndex].m_ve[2].x, pData[uiSiteIndex].m_ve[2].y
    );
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Sub(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulCompC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulRealC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionKS(const deviceSU3Vector* __restrict__ pMe, const deviceSU3Vector* __restrict__ pOther, CLGComplex* result)
{
    intokernal;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplexKS(deviceSU3Vector* pMe, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].MulComp(a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealKS(deviceSU3Vector* pMe, Real a)
{
    intokernal;
    pMe[uiSiteIndex].MulReal(a);
}

/**
* Initialize
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionKS(deviceSU3Vector* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernalInt4;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeOneSU3Vector();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeRandomGaussian(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeRandomZ4(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceKS(deviceSU3Vector* pDeviceData, UINT uiDesiredSite, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeOneSU3VectorColor(byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceKS(deviceSU3Vector* pDeviceData, UINT uiDesiredT, BYTE byColor, BYTE byFieldID, UINT uiVolumn)
{
    intokernalInt4;
    const Real fDenominator = F(1.0) / uiVolumn;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][uiBigIdx];

    if (uiDesiredT == sSite4.w && !sIdx.IsDirichlet())
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeOneSU3VectorColor(byColor);
        pDeviceData[uiSiteIndex].MulReal(fDenominator);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
}

#pragma endregion

CFieldFermionKS::CFieldFermionKS()
    : CFieldFermion()
    , m_fKai(F(0.125))
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount));
}

CFieldFermionKS::~CFieldFermionKS()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

/**
*
*/
void CFieldFermionKS::InitialField(EFieldInitialType eInitialType)
{
    preparethread;

    _kernelInitialFermionKS << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

void CFieldFermionKS::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionKS::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 6 * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    InitialWithByte(data);
    free(data);
}

void CFieldFermionKS::InitialWithByte(BYTE* byData)
{
    deviceSU3Vector* readData = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * m_uiSiteCount);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real thisSite[24];
        memcpy(thisSite, byData + i * sizeof(Real) * 24, sizeof(Real) * 24);
        for (UINT k = 0; k < 3; ++k)
        {
            readData[i].m_ve[k] = _make_cuComplex(
                thisSite[2 * k],
                thisSite[2 * k + 1]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldFermionKS::InitialOtherParameters(CParameters& params)
{
    params.FetchValueReal(_T("Hopping"), m_fKai);
    if (m_fKai < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionKS: Kai is nearly 0, such that Dphi \approx phi! This will cause problem!\n"));
    }
    CCommonData::m_fKai = m_fKai;
}

void CFieldFermionKS::DebugPrintMe() const
{
    preparethread;
    _kernelPrintFermionKS << <block, threads >> > (m_pDeviceData);
}

void CFieldFermionKS::CopyTo(CField* U) const
{
    if (NULL == U || EFT_FermionStaggered != U->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return;
    }

    CField::CopyTo(U);

    CFieldFermionKS* pField = dynamic_cast<CFieldFermionKS*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    pField->m_byFieldId = m_byFieldId;
    pField->m_fKai = m_fKai;
}

void CFieldFermionKS::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggered != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return;
    }
    const CFieldFermionKS* pField = dynamic_cast<const CFieldFermionKS*>(x);

    preparethread;
    _kernelAxpyPlusFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKS::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggered != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return;
    }
    const CFieldFermionKS* pField = dynamic_cast<const CFieldFermionKS*>(x);

    preparethread;
    _kernelAxpyMinusFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKS::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggered != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return;
    }
    const CFieldFermionKS* pField = dynamic_cast<const CFieldFermionKS*>(x);

    preparethread;
    _kernelAxpyRealFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionKS::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggered != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return;
    }
    const CFieldFermionKS* pField = dynamic_cast<const CFieldFermionKS*>(x);

    preparethread;
    _kernelAxpyComplexFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

CLGComplex CFieldFermionKS::Dot(const CField* x) const
{
    if (NULL == x || EFT_FermionStaggered != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only copy to CFieldFermionKS!"));
        return _make_cuComplex(0, 0);
    }
    const CFieldFermionKS* pField = dynamic_cast<const CFieldFermionKS*>(x);

    preparethread;
    _kernelDotFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldFermionKS::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyComplexKS << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKS::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealKS << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKS::ApplyGamma(EGammaMatrix eGamma)
{
    appCrucial(_T("Not implemented yet...\n"));
}
/**
* generate phi by gaussian random.
* phi = D phi
*/
void CFieldFermionKS::PrepareForHMC(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    preparethread;
    _kernelInitialFermionKS << <block, threads >> > (
        pPooled->m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));

    pPooled->Return();

    if (NULL != appGetFermionSolver() && !appGetFermionSolver()->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }
    //cache a inverse DDdagger field
    if (CCommonData::m_bStoreLastSolution)
    {
        CFieldCache* pCache = appGetLattice()->m_pFieldCache;
        CField* pField = pCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField);
        if (NULL == pField)
        {
            pField = GetCopy();
            pCache->CacheField(CFieldCache::CachedInverseDDdaggerField, pField);
        }
        else
        {
            CopyTo(pField);
        }
    }
}

//Kai should be part of D operator
void CFieldFermionKS::D(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

//Kai should be part of D operator
void CFieldFermionKS::Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

void CFieldFermionKS::DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKS can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

UBOOL CFieldFermionKS::InverseD(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_D);
}

UBOOL CFieldFermionKS::InverseDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_Ddagger);
}

UBOOL CFieldFermionKS::InverseDDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_DDdagger);
}

UBOOL CFieldFermionKS::CalculateForce(
    const CFieldGauge* pGauge,
    CFieldGauge* pForce,
    ESolverPhase ePhase) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);

    CField* pDDaggerPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField* pDPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField* pCachedField = CCommonData::m_bStoreLastSolution ?
        appGetLattice()->m_pFieldCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField)
        : NULL;

    if (NULL == pDDaggerPhi || EFT_FermionWilsonSquareSU3 != pDDaggerPhi->GetFieldType()
        || NULL == pDPhi || EFT_FermionWilsonSquareSU3 != pDPhi->GetFieldType())
    {
        appCrucial(_T("Pooled field not found!\n"));
        if (NULL != pDDaggerPhi)
        {
            pDDaggerPhi->Return();
        }
        if (NULL != pDPhi)
        {
            pDPhi->Return();
        }
        return FALSE;
    }
    CFieldFermionWilsonSquareSU3* pDDaggerPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pDDaggerPhi);
    CFieldFermionWilsonSquareSU3* pDPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pDPhi);
    //if (!pDDaggerPhiWilson->InverseDDdagger(pGaugeSU3))
    if (!appGetFermionSolver()->Solve(
        pDDaggerPhiWilson, this, pGaugeSU3,
        EFO_F_DDdagger, ePhase, pCachedField))
    {
        appCrucial(_T("Sparse Linear Solver failed...\n"));
        pDDaggerPhi->Return();
        pDPhi->Return();
        return FALSE;
    }
    //phi 2 = D^{-1}phi = D+ (DD+)^{-1} phi
    //It is faster to calcuate D+ phi2 then D^{-1} phi
    pDDaggerPhiWilson->CopyTo(pDPhiWilson);
    if (NULL != pCachedField)
    {
        //The gauge field is changing slowly, and D depends only on gauge, also change slowly
        //Use the last solution as start point will accelerate the solver, so we cache it
        pDDaggerPhiWilson->CopyTo(pCachedField);
    }
    pDPhiWilson->Ddagger(pGaugeSU3);
    DerivateDOperator(
        pForceSU3->m_pDeviceData,
        pDPhiWilson->m_pDeviceData,
        pDDaggerPhiWilson->m_pDeviceData,
        pGaugeSU3->m_pDeviceData);

    pDDaggerPhi->Return();
    pDPhi->Return();

    return TRUE;
}

void CFieldFermionKS::InitialAsSource(const SFermionSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        _kernelMakePointSourceKS << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_byColorIndex);
    }
    case EFS_Wall:
    {
        preparethread;
        _kernelMakeWallSourceKS << <block, threads >> > (
            m_pDeviceData,
            sourceData.m_sSourcePoint.w,
            sourceData.m_byColorIndex,
            m_byFieldId,
            appGetLattice()->m_pIndexCache->m_uiSiteXYZ);
    }
    break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

void CFieldFermionKS::SetKai(Real fKai)
{
    m_fKai = fKai;
    CCommonData::m_fKai = fKai;
}

void CFieldFermionKS::SaveToFile(const CCString& fileName) const
{
    UINT uiSize = 0;
    BYTE* saveData = CopyDataOut(uiSize);
    appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    free(saveData);
}

BYTE* CFieldFermionKS::CopyDataOut(UINT& uiSize) const
{
    deviceSU3Vector* toSave = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 6);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real oneSite[6];
        for (UINT k = 0; k < 3; ++k)
        {
            oneSite[2 * k] = static_cast<Real>(toSave[i].m_ve[k].x);
            oneSite[2 * k + 1] = static_cast<Real>(toSave[i].m_ve[k].y);
        }
        memcpy(saveData + sizeof(Real) * i * 6, oneSite, sizeof(Real) * 6);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

TArray<CFieldFermion*> CFieldFermionKS::GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const
{
    TArray<CFieldFermion*> ret;
    for (UINT j = 0; j < 12; ++j)
    {
        ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
        if (NULL == ret[j])
        {
            appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
            _FAIL_EXIT;
        }
    }

    for (BYTE s = 0; s < 4; ++s)
    {
        for (BYTE c = 0; c < 3; ++c)
        {
            SFermionSource sourceData;
            sourceData.m_eSourceType = EFS_Point;
            sourceData.m_sSourcePoint = site;
            sourceData.m_byColorIndex = c;
            sourceData.m_bySpinIndex = s;

            ret[s * 3 + c]->InitialAsSource(sourceData);

            if (NULL != appGetFermionSolver() && !appGetFermionSolver()->IsAbsoluteAccuracy())
            {
                ret[s * 3 + c]->m_fLength = ret[s * 3 + c]->Dot(ret[s * 3 + c]).x;
            }
            ret[s * 3 + c]->InverseD(pGauge);
        }
    }
    return ret;
}

CCString CFieldFermionKS::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKS\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

#pragma region Field Matrix Operation

CFieldMatrixOperationKS::CFieldMatrixOperationKS()
{
    m_pHostResBuffer = (deviceSU3Vector**)malloc(sizeof(deviceSU3Vector*) * _kFieldMatrixMaxDim);
    m_pHostLeftBuffer = (deviceSU3Vector**)malloc(sizeof(deviceSU3Vector*) * _kFieldMatrixMaxDim);
    checkCudaErrors(cudaMalloc((void**)&m_pResBuffer, sizeof(deviceSU3Vector*) * _kFieldMatrixMaxDim));
    checkCudaErrors(cudaMalloc((void**)&m_pLeftBuffer, sizeof(deviceSU3Vector*) * _kFieldMatrixMaxDim));
}

CFieldMatrixOperationKS::~CFieldMatrixOperationKS()
{
    free(m_pHostResBuffer);
    free(m_pHostLeftBuffer);

    checkCudaErrors(cudaFree(m_pResBuffer));
    checkCudaErrors(cudaFree(m_pLeftBuffer));
}

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* V(v1,v2,...,vk) = W(w1,w2,...,wm) (m11, ..., m1k)
*                                   (..., ..., ...)
*                                   (mm1, ..., mmk)
* I think this is expansive... the FLOP of Ax is about 100n, but this has m x k x n
*/
__global__ void _CLG_LAUNCH_BOUND_(_MATRIX_BOUND)
_kernelMatrixMultiply(
    deviceSU3Vector** pRes,
    deviceSU3Vector** pLeft,
    const CLGComplex* __restrict__ pMatrix,
    UINT uiDimX, UINT uiDimY) //x=m,y=k
{
    intokernalE(3);

    CLGComplex result[CFieldMatrixOperation::_kFieldMatrixMaxDim];

    for (UINT i = 0; i < uiDimY; ++i)
    {
        result[i] = _make_cuComplex(F(0.0), F(0.0));
        for (UINT j = 0; j < uiDimX; ++j)
        {
            result[i] = _cuCaddf(result[i], _cuCmulf(
                j < uiDimY ?
                pRes[j][uiSiteIndex].m_ve[elementIdx]
                : pLeft[j - uiDimY][uiSiteIndex].m_ve[elementIdx],
                pMatrix[j * uiDimY + i]
            ));
        }
    }

    for (UINT i = 0; i < uiDimY; ++i)
    {
        pRes[i][uiSiteIndex].m_ve[elementIdx] = result[i];
    }
}

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* v1 = (m11, ..., m1m)  w1
* ..   (..., ..., ...)  ..
* vk   (mk1, ..., mkm)  wk
*                       wk+1
*                       ...
*                       wm
*/

void CFieldMatrixOperationKS::VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY; ++i)
    {
        CFieldFermionKS* pF = dynamic_cast<CFieldFermionKS*>(res[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKS only work with CFieldFermionKS!\n"));
            return;
        }
        m_pHostResBuffer[i] = pF->m_pDeviceData;
    }

    for (UINT i = 0; i < uiDimX - uiDimY; ++i)
    {
        const CFieldFermionKS* pF = dynamic_cast<const CFieldFermionKS*>(left[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKS only work with CFieldFermionKS!\n"));
            return;
        }
        m_pHostLeftBuffer[i] = pF->m_pDeviceData;
    }

    checkCudaErrors(cudaMemcpy(m_pResBuffer, m_pHostResBuffer, sizeof(deviceSU3Vector*) * uiDimY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pLeftBuffer, m_pHostLeftBuffer, sizeof(deviceSU3Vector*) * (uiDimX - uiDimY), cudaMemcpyHostToDevice));

    preparethreadE(3);
    _kernelMatrixMultiply << <block, threads >> > (m_pResBuffer, m_pLeftBuffer, deviceMatrix, uiDimX, uiDimY);
}

#pragma endregion

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================