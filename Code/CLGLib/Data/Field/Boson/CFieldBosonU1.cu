//=============================================================================
// FILENAME : CFieldBosonU1NoGauge.h
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [04/23/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBosonU1)

#pragma region kernel D and force

/**
* sum_mu (U_mu(n)phi(n+mu) + U_{-mu}(n)phi(n-mu))
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDBoson(
    const CLGComplex* __restrict__ pDeviceData,
    const CLGComplex* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pBosonMove,
    CLGComplex* pResultData,
    EOperatorCoefficientType eCoeff,
    Real fRealCoeff, 
    CLGComplex cCmpCoeff,
    BYTE byFieldId)
{
    intokernaldir;

    CLGComplex result = _zeroc;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Boson = pBosonMove[2 * linkIndex];
        const SIndex& x_m_mu_Boson = pBosonMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const CLGComplex& x_Gauge_element = NULL == pGauge ? _onec : pGauge[linkIndex];
        CLGComplex x_m_mu_Gauge_element = NULL == pGauge ? _onec : pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element = _cuConjf(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        CLGComplex u_phi_x_p_m = _cuCmulf(x_Gauge_element, pDeviceData[x_p_mu_Boson.m_uiSiteIndex]);

        //U^{dagger}(x-mu) phi(x-mu)
        CLGComplex u_dagger_phi_x_m_m = _cuCmulf(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Boson.m_uiSiteIndex]);

        result = _cuCaddf(result, u_phi_x_p_m);
        result = _cuCaddf(result, u_dagger_phi_x_m_m);
    }

    pResultData[uiSiteIndex] = result;

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex] = cuCmulf_cr(pResultData[uiSiteIndex], fRealCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex] = _cuCmulf(pResultData[uiSiteIndex], cCmpCoeff);
        break;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSquareBoson(
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pResultData)
{
    intokernal;
    pResultData[uiSiteIndex] = cuCmulf_cr(pDeviceData[uiSiteIndex], __cuCabsSqf(pDeviceData[uiSiteIndex]));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDBosonForceU1(
    const CLGComplex* __restrict__ pPhi,
    const CLGComplex* __restrict__ pGauge,
    const SIndex* __restrict__ pBosonMove,
    CLGComplex* pForce,
    BYTE byFieldId)
{
    intokernaldir;

    const CLGComplex& phi_dagger = _cuConjf(pPhi[uiSiteIndex]);

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Boson = pBosonMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        const CLGComplex& phi_p_mu = pPhi[x_p_mu_Boson.m_uiSiteIndex];

        const CLGComplex& x_Gauge_element = pGauge[linkIndex];

        //U phi(n+mu)phi^+(n)
        CLGComplex forceOfThisLink = _cuCmulf(x_Gauge_element, _cuCmulf(phi_p_mu, phi_dagger));

        //TA
        //pForce[linkIndex] = _cuCsubf(pForce[linkIndex], _make_cuComplex(F(0.0), forceOfThisLink.x));
        pForce[linkIndex].y = pForce[linkIndex].y - forceOfThisLink.x;
    }
}

#pragma endregion

void CFieldBosonU1::D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);

    if (NULL != pGauge && EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = (NULL == pGauge) ? NULL : dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldBosonU1* pPooled = dynamic_cast<CFieldBosonU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    _kernelDBoson << <block, threads >> > (
        pPooled->m_pDeviceData,
        (NULL == pFieldU1) ? NULL : pFieldU1->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        m_pDeviceData,
        eCoeffType,
        fRealCoeff,
        cCompCoeff,
        m_byFieldId);

    checkCudaErrors(cudaDeviceSynchronize());
    pPooled->Return();
}

void CFieldBosonU1::ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldGauge** pGaugeForce, const CFieldBoson* const* pBoson)
{

}

#pragma region other kernels

/**
* Initialize
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialBoson(CLGComplex* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernal;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _zeroc;
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _onec;
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = _deviceRandomGaussC(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = _deviceRandomZ4(_deviceGetFatIndex(uiSiteIndex, 0));
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
_kernelAddBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSubBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCsubf(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMulBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCmulf(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelConjMulBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCmulf(_cuConjf(pMe[uiSiteIndex]), pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], _cuCmulf(pOther[uiSiteIndex], a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMulComplexBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, UBOOL bConj)
{
    intokernal;
    if (bConj)
    {
        pMe[uiSiteIndex].y = -pMe[uiSiteIndex].y;
    }
    pMe[uiSiteIndex] = _cuCmulf(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealBoson(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], cuCmulf_cr(pOther[uiSiteIndex], a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotBoson(const CLGComplex* __restrict__ pMe, const CLGComplex* __restrict__ pOther,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex* result
#endif
)
{
    intokernal;
#if !_CLG_DOUBLEFLOAT
    result[uiSiteIndex] = _cToDouble(_cuCmulf(_cuConjf(pMe[uiSiteIndex]), pOther[uiSiteIndex]));
#else
    result[uiSiteIndex] = _cuCmulf(_cuConjf(pMe[uiSiteIndex]), pOther[uiSiteIndex]);
#endif
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyBoson(CLGComplex* pMe, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCmulf(pMe[uiSiteIndex], a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealBoson(CLGComplex* pMe, Real a)
{
    intokernal;
    pMe[uiSiteIndex] = cuCmulf_cr(pMe[uiSiteIndex], a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBosonConjugate(CLGComplex* pDeviceData)
{
    intokernal;
    pDeviceData[uiSiteIndex] = _cuConjf(pDeviceData[uiSiteIndex]);
}


#pragma endregion

CFieldBosonU1::CFieldBosonU1()
    : CFieldBoson()
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount));
}

CFieldBosonU1::~CFieldBosonU1()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldBosonU1::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialBoson << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

void CFieldBosonU1::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin
     && eFieldType != EFFT_CLGBinFloat
     && eFieldType != EFFT_CLGBinDouble)
    {
        appCrucial(_T("CFieldBosonU1::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * m_uiSiteCount);
    if (eFieldType == EFFT_CLGBinFloat)
    {
        uiSize = static_cast<UINT>(sizeof(FLOAT) * 2 * m_uiSiteCount);
    }
    else if (eFieldType == EFFT_CLGBinDouble)
    {
        uiSize = static_cast<UINT>(sizeof(DOUBLE) * 2 * m_uiSiteCount);
    }
    UINT uiReadSize = uiSize;
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiReadSize);
    if (NULL == data)
    {
        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
        _FAIL_EXIT;
    }

    if (uiSize != uiReadSize)
    {
        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), uiSize, uiReadSize);
        _FAIL_EXIT;
    }

#if _CLG_DOUBLEFLOAT
    if (eFieldType == EFFT_CLGBinFloat)
    {
        FLOAT* data1 = (FLOAT*)data;
        DOUBLE* data2 = (DOUBLE*)(malloc(sizeof(DOUBLE) * 2 * m_uiSiteCount));
        for (UINT i = 0; i < m_uiSiteCount * 2; ++i)
        {
            data2[i] = static_cast<DOUBLE>(data1[i]);
        }
        InitialWithByte((BYTE*)data2);
        free(data);
        free(data2);
    }
#else
    if (eFieldType == EFFT_CLGBinDouble)
    {
        DOUBLE* data1 = (DOUBLE*)data;
        FLOAT* data2 = (FLOAT*)(malloc(sizeof(FLOAT) * 2 * m_uiSiteCount));
        for (UINT i = 0; i < m_uiSiteCount * 2; ++i)
        {
            data2[i] = static_cast<FLOAT>(data1[i]);
        }
        InitialWithByte((BYTE*)data2);
        free(data);
        free(data2);
    }
#endif
    else
    {
        InitialWithByte(data);
        free(data);
    }

    appCrucial(_T("CFieldBosonU1::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
}

void CFieldBosonU1::InitialWithByte(BYTE* byData)
{
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    Real* hostBuffer = (Real*)byData;
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        for (UINT k = 0; k < 3; ++k)
        {
            readData[i] = _make_cuComplex(
                hostBuffer[2 * k],
                hostBuffer[2 * k + 1]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldBosonU1::InitialOtherParameters(CParameters& params)
{
    CFieldBoson::InitialOtherParameters(params);
}

BYTE* CFieldBosonU1::CopyDataOut(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    Real* fsaveData = (Real*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        fsaveData[2 * i] = toSave[i].x;
        fsaveData[2 * i + 1] = toSave[i].y;
    }
    free(toSave);
    return saveData;
}

BYTE* CFieldBosonU1::CopyDataOutFloat(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    FLOAT* fsaveData = (FLOAT*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        fsaveData[2 * i] = static_cast<FLOAT>(toSave[i].x);
        fsaveData[2 * i + 1] = static_cast<FLOAT>(toSave[i].y);
    }
    free(toSave);
    return saveData;
}

BYTE* CFieldBosonU1::CopyDataOutDouble(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    DOUBLE* fsaveData = (DOUBLE*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        fsaveData[2 * i] = static_cast<DOUBLE>(toSave[i].x);
        fsaveData[2 * i + 1] = static_cast<DOUBLE>(toSave[i].y);
    }
    free(toSave);
    return saveData;
}

void CFieldBosonU1::DebugPrintMe() const
{
    CLGComplex* toprint = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
    {
        if (0 == (uiSite % _HC_Lt))
        {
            appGeneral(_T("\n"));
        }
        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" (%d,%d,%d,%d) = %f %s %f I"),
            site4.x, site4.y, site4.z, site4.w,
            toprint[uiSite].x,
            toprint[uiSite].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(toprint[uiSite].y)
        );
    }
    appPopLogDate();

    appSafeFree(toprint);
}

void CFieldBosonU1::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);

    preparethread;
    _kernelAddBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldBosonU1::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);

    preparethread;
    _kernelSubBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldBosonU1::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);

    preparethread;
    _kernelAxpyRealBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldBosonU1::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);

    preparethread;
    _kernelAxpyComplexBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldBosonU1::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || EFT_BosonU1 != other->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(other);

    preparethread;
    _kernelMulComplexBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, bDagger);
}

cuDoubleComplex CFieldBosonU1::Dot(const CField* x) const
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return make_cuDoubleComplex(0.0, 0.0);
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);
    preparethread;
    _kernelDotBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldBosonU1::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyBoson << <block, threads >> > (m_pDeviceData, a);
}

void CFieldBosonU1::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealBoson << <block, threads >> > (m_pDeviceData, a);
}

void CFieldBosonU1::FieldMultply(const CFieldBoson* x, UBOOL bConj)
{
    if (NULL == x || EFT_BosonU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonU1* pField = dynamic_cast<const CFieldBosonU1*>(x);

    preparethread;
    if (bConj)
    {
        _kernelConjMulBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
    }
    else
    {
        _kernelMulBoson << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
    }
}

void CFieldBosonU1::Dagger()
{
    preparethread;
    _kernelBosonConjugate << <block, threads >> > (m_pDeviceData);
}

void CFieldBosonU1::CopyTo(CField* U) const
{
    if (NULL == U || EFT_BosonU1 != U->GetFieldType())
    {
        appCrucial(_T("EFT_BosonU1 can only copy to EFT_BosonU1!"));
        return;
    }

    CFieldBoson::CopyTo(U);

    CFieldBosonU1* pField = dynamic_cast<CFieldBosonU1*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
}

void CFieldBosonU1::MakeRandomMomentum()
{
    preparethread;
    _kernelInitialBoson << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================