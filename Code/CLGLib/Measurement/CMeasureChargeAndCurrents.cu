//=============================================================================
// FILENAME : CMeasureChargeAndCurrents.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.h"
#include "CMeasureChargeAndCurrents.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChargeAndCurrents)

#pragma region device functions

__device__ void _deviceMeasure0(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    //do nothing
}

__device__ void _deviceMeasure1(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA1].MulWilsonC(element);
}

__device__ void _deviceMeasure2(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA2].MulWilsonC(element);
}

__device__ void _deviceMeasure3(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA3].MulWilsonC(element);
}

__device__ void _deviceMeasure4(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA4].MulWilsonC(element);
}

/**
* -i kappa gamma_4 sigma_12
*/
__device__ void _deviceMeasure5(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[SIGMA12E].MulWilsonC(element);
    element = __chiralGamma[GAMMA4].MulWilsonC(element);
    element.MulComp(_make_cuComplex(F(0.0), -fKappa));
}

/**
* gamma_1 + y Omega gamma_4
*/
__device__ void _deviceMeasure6(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    const Real fYOmega = static_cast<Real>(site.y - _DC_Centery) * fOmega;
    deviceWilsonVectorSU3 toAdd = __chiralGamma[GAMMA4].MulWilsonC(element);
    toAdd.MulReal(fYOmega);
    element = __chiralGamma[GAMMA1].MulWilsonC(element);
    element.Add(toAdd);
}

/**
* gamma_2 - x Omega gamma_4
*/
__device__ void _deviceMeasure7(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    const Real fXOmega = static_cast<Real>(site.x - _DC_Centerx) * fOmega;
    deviceWilsonVectorSU3 toAdd = __chiralGamma[GAMMA4].MulWilsonC(element);
    toAdd.MulReal(fXOmega);
    element = __chiralGamma[GAMMA2].MulWilsonC(element);
    element.Sub(toAdd);
}

/**
* gamma_4 gamma_5
*/
__device__ void _deviceMeasure8(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA4].MulWilsonC(element);
    element = __chiralGamma[GAMMA5].MulWilsonC(element);
}

__constant__ _deviceMeasureFunc _cMeasureFuncs[CMeasureChargeAndCurrents::_kGammaInInterests] =
{ 
    _deviceMeasure0,
    _deviceMeasure1,
    _deviceMeasure2,
    _deviceMeasure3,
    _deviceMeasure4,
    _deviceMeasure5,
    _deviceMeasure6,
    _deviceMeasure7,
    _deviceMeasure8
};

#pragma endregion

#pragma region kernles 

/**
* blockIdx.x = 0 -> _kGammaInInterests - 1
* threadIdx.x = 0 -> 11
* byArrayIdx = 0 -> Lx - 1
*/
__global__ void
_CLG_LAUNCH_(12, 8)
_kernel_Gammas(
    deviceWilsonVectorSU3** pSources,
    SSmallInt4 sSite4,
    Real fKappa,
    Real fOmega,
    BYTE byArrayIdx,
    BYTE byFieldId,
    deviceWilsonVectorSU3* res)
{
    UINT uiMeasureFunction = blockIdx.x;
    UINT uiArrayIndex = byArrayIdx * CMeasureChargeAndCurrents::_kGammaInInterests + uiMeasureFunction;
    //s * 3 + c
    UINT uiC = threadIdx.x;
    UINT uiS = threadIdx.y;
    const UINT uiCS = uiS * 3 + uiC;

    const UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        res[uiArrayIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    deviceWilsonVectorSU3 right_element(pSources[uiCS][uiSiteIndex]);
    (*_cMeasureFuncs[uiMeasureFunction])(sSite4, fKappa, fOmega, right_element);

    //Note that, it is not res[byArrayIdx] = result
    //It is res[c,s] = delta_{cc}delta_ss result[c,s]
    res[uiArrayIndex].m_d[uiS].m_ve[uiC] = right_element.m_d[uiS].m_ve[uiC];
}

/**
* blockIdx.x = 0 -> _kGammaInInterests - 1
* threadIdx.x = 0 -> Lx - 1
*/
__global__ void
_CLG_LAUNCH_(32, CMeasureChargeAndCurrents::_kGammaInInterests)
_kernel_Trace_Gammas(
    const deviceWilsonVectorSU3* __restrict__ pOperator,
    CLGComplex* pResLine)
{
    UINT uiIdx = threadIdx.x * CMeasureChargeAndCurrents::_kGammaInInterests + blockIdx.x;
    pResLine[uiIdx] = pOperator[uiIdx].Sum();
}

#pragma endregion

CMeasureChargeAndCurrents::CMeasureChargeAndCurrents()
    : CMeasure()
    , m_pHostDataBuffer(NULL)
    , m_pDeviceDataBuffer(NULL)
    , m_pOperatorData(NULL)
{
    ////see: https://stackoverflow.com/questions/9000388/device-function-pointers
    //_deviceMeasureFunc* ph_ptr = (_deviceMeasureFunc*)malloc(_kGammaInInterests * sizeof(_deviceMeasureFunc));

    //checkCudaErrors(cudaMalloc((void**)&m_pMeasureFunctions, _kGammaInInterests * sizeof(_deviceMeasureFunc)));

    //checkCudaErrors(cudaMemcpyFromSymbol(&ph_ptr, pfunc1, sizeof(_deviceMeasureFunc) * 8));

    //checkCudaErrors(cudaMemcpy(m_pMeasureFunctions, ph_ptr, _kGammaInInterests * sizeof(_deviceMeasureFunc), cudaMemcpyHostToDevice));
    //free(ph_ptr);
}

CMeasureChargeAndCurrents::~CMeasureChargeAndCurrents()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBuffer));
    }
    if (NULL != m_pOperatorData)
    {
        checkCudaErrors(cudaFree(m_pOperatorData));
    }
    //if (NULL != m_pMeasureFunctions)
    //{
    //    checkCudaErrors(cudaFree(m_pMeasureFunctions));
    //}
}

void CMeasureChargeAndCurrents::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _kGammaInInterests);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBuffer, sizeof(CLGComplex) * _HC_Lx * _kGammaInInterests));
    checkCudaErrors(cudaMalloc((void**)&m_pOperatorData, sizeof(deviceWilsonVectorSU3) * _HC_Lx * _kGammaInInterests));

    Reset();
}

void CMeasureChargeAndCurrents::SourceSanningSingleField(const class CFieldGauge* pGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& sourceSite)
{
    if (12 != sources.Num())
    {
        appCrucial(_T("Wilson Dirac SU3 need 12 sources!\n"));
        return;
    }
    deviceWilsonVectorSU3* pDevicePtr[12];
    for (INT i = 0; i < 12; ++i)
    {
        CFieldFermionWilsonSquareSU3* fermionfield = dynamic_cast<CFieldFermionWilsonSquareSU3*>(sources[i]);
        pDevicePtr[i] = fermionfield->m_pDeviceData;
    }

    dim3 _blocks(_kGammaInInterests, 1, 1);
    dim3 _thread1(12, 1, 1);

    deviceWilsonVectorSU3** ppDevicePtr;
    checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
    checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));

    Real fOmega = F(0.0);
    const CFieldFermionWilsonSquareSU3DR* pRA = dynamic_cast<const CFieldFermionWilsonSquareSU3DR*>(sources[0]);
    if (NULL != pRA)
    {
        fOmega = static_cast<Real>(pRA->GetFermionOmega());
    }
    
    //sourceSite.x = 1 to lx - 1
    _kernel_Gammas << <_blocks, _thread1 >> > (
        ppDevicePtr,
        sourceSite,
        CCommonData::m_fKai,
        fOmega,
        static_cast<BYTE>(sourceSite.x), //array idx
        GetFermionFieldId(),
        //m_pMeasureFunctions,
        m_pOperatorData
        );

    checkCudaErrors(cudaFree(ppDevicePtr));

    if (sourceSite.x == static_cast<SBYTE>(_HC_Lx) - 1)
    {
        //all sites calculated
        ++m_uiConfigurationCount;
        dim3 _thread2(_HC_Lx, 1, 1);
        _kernel_Trace_Gammas << <_blocks, _thread2 >> > (m_pOperatorData, m_pDeviceDataBuffer);

        checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(CLGComplex) * _HC_Lx * _kGammaInInterests, cudaMemcpyDeviceToHost));

        for (UINT j = 0; j < _kGammaInInterests; ++j)
        {
            for (UINT i = 1; i < _HC_Lx; ++i)
            {
                m_lstAllRes.AddItem(m_pHostDataBuffer[i * _kGammaInInterests + j]);
            }
        }

        if (m_bShowResult)
        {
            appDetailed(_T("\n\n ==================== Densities (%d con)============================ \n\n"), m_uiConfigurationCount);

            static CCString msgs[_kGammaInInterests] =
            {
                _T(" ----------- Chiral condensate ------------- \n"),
                _T(" ----------- J1 ------------- \n"),
                _T(" ----------- J2 ------------- \n"),
                _T(" ----------- J3 ------------- \n"),
                _T(" ----------- J4 ------------- \n"),

                _T(" ----------- sigma12 ------------- \n"),

                _T(" ----------- Jx ------------- \n"),
                _T(" ----------- Jy ------------- \n"),

                _T(" ----------- n5 ------------- \n"),
            };
            

            for (UINT j = 0; j < _kGammaInInterests; ++j)
            {
                appDetailed(msgs[j]);
                for (UINT i = 1; i < _HC_Lx; ++i)
                {
                    appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, 
                        m_pHostDataBuffer[i * _kGammaInInterests + j].x,
                        m_pHostDataBuffer[i * _kGammaInInterests + j].y);
                }

                appDetailed(_T("\n"));
            }

            appDetailed(_T("\n\n ================================================ \n\n"));
        }
    }
}

void CMeasureChargeAndCurrents::Report()
{
    appPushLogDate(FALSE);

    assert(m_uiConfigurationCount * (_HC_Lx - 1) * _kGammaInInterests == static_cast<UINT>(m_lstAllRes.Num()));
    TArray<Real> tmpSum;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Charge Current Densities (%d con)============================\n"), m_uiConfigurationCount);
    appGeneral(_T("==========================================================================\n"));

    static CCString msgs[_kGammaInInterests] =
    {
        _T(" ----------- Chiral condensate ------------- \n"),
        _T(" ----------- J1 ------------- \n"),
        _T(" ----------- J2 ------------- \n"),
        _T(" ----------- J3 ------------- \n"),
        _T(" ----------- J4 ------------- \n"),

        _T(" ----------- sigma12 ------------- \n"),

        _T(" ----------- Jx ------------- \n"),
        _T(" ----------- Jy ------------- \n"),

        _T(" ----------- n5 ------------- \n"),
    };

    for (UINT n = 0; n < _kGammaInInterests; ++n)
    {
        appGeneral(msgs[n]);

        appGeneral(_T("\n{\n"));

        TArray<Real> average;
        for (UINT i = 0; i < (_HC_Lx - 1); ++i)
        {
            average.AddItem(F(0.0));
        }

        for (UINT j = 0; j < m_uiConfigurationCount; ++j)
        {
            appGeneral(_T("{"));
            for (UINT i = 0; i < (_HC_Lx - 1); ++i)
            {
                const UINT uiThisSiteThisFunc = n * (_HC_Lx - 1) + i;
                const UINT uiOneConf = _kGammaInInterests * (_HC_Lx - 1);
                const Real fV = m_lstAllRes[uiOneConf * j + uiThisSiteThisFunc].x;
                average[i] = average[i] + fV;
                appGeneral(_T("%2.12f, "), fV);
            }
            appGeneral(_T("},\n"));
        }

        appGeneral(_T("}\n\n"));

        appDetailed(msgs[n]);
        appGeneral(_T("----------- Average -----------\n"));
        appGeneral(_T("{"));
        for (UINT i = 0; i < (_HC_Lx - 1); ++i)
        {
            appGeneral(_T("%2.12f, "), average[i] / m_uiConfigurationCount);
        }
        appGeneral(_T("}\n\n"));
    }

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));

    appPopLogDate();
}

void CMeasureChargeAndCurrents::Reset()
{
    CMeasure::Reset();
    m_lstAllRes.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================