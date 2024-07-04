//=============================================================================
// FILENAME : CFieldBosonvn.h
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [07/04/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE


#pragma region kernel D and force

/**
* sum_mu (U_mu(n)phi(n+mu) + U_{-mu}(n)phi(n-mu))
*/
template<typename deviceDataBoson, typename deviceDataGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDBosonVN(
    const deviceDataBoson* __restrict__ pDeviceData,
    const deviceDataGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pBosonMove,
    deviceDataBoson* pResultData,
    EOperatorCoefficientType eCoeff,
    Real fRealCoeff, 
    CLGComplex cCmpCoeff,
    BYTE byFieldId)
{
    intokernaldir;

    deviceDataBoson result = _makeZero<deviceDataBoson>();

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
        const deviceDataGauge& x_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[linkIndex];
        deviceDataGauge x_m_mu_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceDataBoson u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Boson.m_uiSiteIndex]);

        //U^{dagger}(x-mu) phi(x-mu)
        deviceDataBoson u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Boson.m_uiSiteIndex]);

        _add(result, u_phi_x_p_m);
        _add(result, u_dagger_phi_x_m_m);
    }

    pResultData[uiSiteIndex] = result;

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fRealCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCmpCoeff);
        break;
    }
}


template<typename deviceDataBoson, typename deviceDataGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDBosonForceVN(
    const deviceDataBoson* __restrict__ pBoson,
    const deviceDataGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pBosonMove,
    deviceDataGauge* pGaugeForce,
    BYTE byGaugeFieldId,
    BYTE byFieldId)
{
    intokernaldir;

    deviceDataBoson phi_dagger = _daggerC(pBoson[uiSiteIndex]);

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Boson = pBosonMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        const deviceDataBoson& phi_p_mu = pBoson[x_p_mu_Boson.m_uiSiteIndex];

        const deviceDataGauge& x_Gauge_element = pGauge[linkIndex];

        //U phi(n+mu)phi^+(n)
        deviceDataGauge forceOfThisLink = _makeContract<deviceDataGauge, deviceDataBoson>(_mulVec<deviceDataGauge, deviceDataBoson>(x_Gauge_element, phi_p_mu), phi_dagger);

        //TA
        //pForce[linkIndex] = _cuCsubf(pForce[linkIndex], _make_cuComplex(F(0.0), forceOfThisLink.x));
        //pGaugeForce[linkIndex].y = pGaugeForce[linkIndex].y - F(1.0) * forceOfThisLink.y;
        _ta<deviceDataGauge>(forceOfThisLink);
        _sub(pGaugeForce[linkIndex], forceOfThisLink);
    }
}


#pragma endregion

#pragma region other kernels

/**
* Initialize
*/
template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialBosonVN(deviceDataBoson* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernal;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _makeZero<deviceDataBoson>();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _makeId<deviceDataBoson>();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = _makeGaussian<deviceDataBoson>(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = _makeZ4<deviceDataBoson>(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther)
{
    intokernal;
    _add(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelSubBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther)
{
    intokernal;
    _sub(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther)
{
    intokernal;
    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulComplexBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther, UBOOL bConj)
{
    intokernal;
    if (bConj)
    {
        _dagger(pMe[uiSiteIndex]);
    }
    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealBosonVN(deviceDataBoson* pMe, const deviceDataBoson* __restrict__ pOther, Real a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotBosonVN(const deviceDataBoson* __restrict__ pMe, const deviceDataBoson* __restrict__ pOther, cuDoubleComplex* result)
{
    intokernal;
    result[uiSiteIndex] = _cToDouble(_dot(pMe[uiSiteIndex], pOther[uiSiteIndex]));
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyBosonVN(deviceDataBoson* pMe, CLGComplex a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealBosonVN(deviceDataBoson* pMe, Real a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonConjugateVN(deviceDataBoson* pDeviceData)
{
    intokernal;
    _dagger(pDeviceData[uiSiteIndex]);
}


#pragma endregion

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);

    if (NULL != pGauge && GetFieldType() != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldBosonVU can only play with gauge UN!"));
        return;
    }

    //try external gauge field
    if (NULL == pGauge && m_byGaugeFieldIds.Num() > 0)
    {
        const CField* externelfield = appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]);
        if (NULL != externelfield)
        {
            pGauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]));
            if (pGauge->IsDynamic())
            {
                appCrucial(_T("CFieldBosonUN: A dynamic field is configured for this UN, but not for the action!\n"));
                pGauge = NULL;
            }
        }
    }

    CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pPooled = dynamic_cast<CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    _kernelDBosonVN << <block, threads >> > (
        pPooled->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
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

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const
{
    if (m_byGaugeFieldIds.Num() < 1)
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }
    INT gaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pGauge, m_byGaugeFieldIds[0]);
    if (gaugeidx < 0 || gaugeidx >= m_byGaugeFieldIds.Num())
    {
        appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
        return;
    }

    const CFieldGauge* gauge = pGauge[gaugeidx];
    CFieldGauge* gaugeforce = pGaugeForce[gaugeidx];

    if (NULL == gauge || GetFieldType() != gauge->GetFieldType())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    if (NULL == gaugeforce || GetFieldType() != gaugeforce->GetFieldType())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    preparethread;
    _kernelDBosonForceVN << <block, threads >> > (
        m_pDeviceData,
        (const deviceDataGauge*)gauge->GetData(),
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        (deviceDataGauge*)gaugeforce->GetData(),
        gauge->m_byFieldId,
        m_byFieldId
        );

    //gaugeu1->DebugPrintMe();
    checkCudaErrors(cudaDeviceSynchronize());
}

template<typename deviceDataBoson, typename deviceDataGauge>
CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CFieldBosonVN()
    : CFieldBoson()
    //, m_fCharge(F(0.1))
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount));
}

template<typename deviceDataBoson, typename deviceDataGauge>
CFieldBosonVN<deviceDataBoson, deviceDataGauge>::~CFieldBosonVN()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialBosonVN << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin
     && eFieldType != EFFT_CLGBinFloat
     && eFieldType != EFFT_CLGBinDouble)
    {
        appCrucial(_T("CFieldBosonU1::InitialFieldWithFile: Not support %s File\n"), __ENUM_TO_STRING(EFieldFileType, eFieldType).c_str());
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * FloatN() * m_uiSiteCount);
    if (eFieldType == EFFT_CLGBinFloat)
    {
        uiSize = static_cast<UINT>(sizeof(FLOAT) * FloatN() * m_uiSiteCount);
    }
    else if (eFieldType == EFFT_CLGBinDouble)
    {
        uiSize = static_cast<UINT>(sizeof(DOUBLE) * FloatN() * m_uiSiteCount);
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
        DOUBLE* data2 = (DOUBLE*)(malloc(sizeof(DOUBLE) * FloatN() * m_uiSiteCount));
        for (UINT i = 0; i < m_uiSiteCount * FloatN(); ++i)
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
        FLOAT* data2 = (FLOAT*)(malloc(sizeof(FLOAT) * FloatN() * m_uiSiteCount));
        for (UINT i = 0; i < m_uiSiteCount * FloatN(); ++i)
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

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialWithByte(BYTE* byData)
{
    deviceDataBoson* readData = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * m_uiSiteCount);
    Real* hostBuffer = (Real*)byData;
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        for (UINT j = 0; j < FloatN(); ++j)
        {
            _setelement(readData[i], j, hostBuffer[FloatN() * i + j]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

template<typename deviceDataBoson, typename deviceDataGauge>
BYTE* CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyDataOut(UINT& uiSize) const
{
    deviceDataBoson* toSave = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * FloatN());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    Real* fsaveData = (Real*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        for (UINT j = 0; j < FloatN(); ++j)
        {
            fsaveData[FloatN() * i + j] = _element(toSave[i], j);
        }
    }
    free(toSave);
    return saveData;
}

template<typename deviceDataBoson, typename deviceDataGauge>
BYTE* CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyDataOutFloat(UINT& uiSize) const
{
    deviceDataBoson* toSave = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * FloatN());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    FLOAT* fsaveData = (FLOAT*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        for (UINT j = 0; j < FloatN(); ++j)
        {
            fsaveData[FloatN() * i + j] = static_cast<FLOAT>(_element(toSave[i], j));
        }
    }
    free(toSave);
    return saveData;
}

template<typename deviceDataBoson, typename deviceDataGauge>
BYTE* CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyDataOutDouble(UINT& uiSize) const
{
    deviceDataBoson* toSave = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * FloatN());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    DOUBLE* fsaveData = (DOUBLE*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        for (UINT j = 0; j < FloatN(); ++j)
        {
            fsaveData[FloatN() * i + j] = static_cast<DOUBLE>(_element(toSave[i], j));
        }
    }
    free(toSave);
    return saveData;
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DebugPrintMe() const
{
    deviceDataBoson* toprint = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * m_uiSiteCount);
    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
    {
        if (0 == (uiSite % _HC_Lt))
        {
            appGeneral(_T("\n"));
        }
        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" (%d,%d,%d,%d) = %s, "),
            site4.x, site4.y, site4.z, site4.w,
            appToString(toprint[uiSite]).c_str());
    }
    appPopLogDate();

    appSafeFree(toprint);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::AxpyPlus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);

    preparethread;
    _kernelAddBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::AxpyMinus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);

    preparethread;
    _kernelSubBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Axpy(Real a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);

    preparethread;
    _kernelAxpyRealBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);

    preparethread;
    _kernelAxpyComplexBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(other);

    preparethread;
    _kernelMulComplexBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, bDagger);
}

template<typename deviceDataBoson, typename deviceDataGauge>
cuDoubleComplex CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Dot(const CField* x) const
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return make_cuDoubleComplex(0.0, 0.0);
    }
    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(x);
    preparethread;
    _kernelDotBosonVN << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyBosonVN << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealBosonVN << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Dagger()
{
    preparethread;
    _kernelBosonConjugateVN << <block, threads >> > (m_pDeviceData);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyTo(CField* U) const
{
    if (NULL == U || GetFieldType() != U->GetFieldType())
    {
        appCrucial(_T("EFT_BosonU1 can only copy to EFT_BosonU1!"));
        return;
    }

    CFieldBoson::CopyTo(U);

    CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pField = dynamic_cast<CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceDataBoson) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::MakeRandomMomentum()
{
    if (m_bConstant)
    {
        Zero();
        return;
    }

    preparethread;
    _kernelInitialBosonVN << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);
}

__CLGIMPLEMENT_CLASS(CFieldBosonU1)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================