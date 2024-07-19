//=============================================================================
// FILENAME : CFieldBosonVNKernel.h
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [07/20/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldBosonVNKernel.h"

__BEGIN_NAMESPACE

#if 0


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
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = _makeZero<deviceDataBoson>();
        return;
    }

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
        if (!x_p_mu_Boson.IsDirichlet())
        {
            const deviceDataGauge& x_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[linkIndex];
            //U(x,mu) phi(x+ mu)
            const deviceDataBoson u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Boson.m_uiSiteIndex]);
            _add(result, u_phi_x_p_m);
        }
        
        if (!x_m_mu_Boson.IsDirichlet())
        {
            deviceDataGauge x_m_mu_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
            if (x_m_mu_Gauge.NeedToDagger())
            {
                _dagger(x_m_mu_Gauge_element);
            }

            //U^{dagger}(x-mu) phi(x-mu)
            const deviceDataBoson u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Boson.m_uiSiteIndex]);
            _add(result, u_dagger_phi_x_m_m);
        }
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
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const deviceDataBoson& phi_dagger = pBoson[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Boson = pBosonMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));
        if (x_p_mu_Boson.IsDirichlet())
        {
            continue;
        }

        const deviceDataBoson& phi_p_mu = pBoson[x_p_mu_Boson.m_uiSiteIndex];

        const deviceDataGauge& x_Gauge_element = pGauge[linkIndex];

        //U phi(n+mu)phi^+(n)
        deviceDataGauge forceOfThisLink = _makeContract<deviceDataGauge, deviceDataBoson>(phi_dagger, _mulVec<deviceDataGauge, deviceDataBoson>(x_Gauge_element, phi_p_mu));

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
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    //if (sIdx.IsDirichlet())
    //{
    //    pDevicePtr[uiSiteIndex] = _makeZero<deviceDataBoson>();
    //    return;
    //}

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
    case EFIT_RandomGenerator:
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
        printf("Wilson Fermion Field cannot be initialized with this type! %d\n", eInitialType);
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
_kernelElementVN(const deviceDataBoson* __restrict__ pMe, UINT idx, DOUBLE* result)
{
    intokernal;
    result[uiSiteIndex] = static_cast<DOUBLE>(_element(pMe[uiSiteIndex], idx));
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

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonFixBoundary(deviceDataBoson* pDeviceData, BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pDeviceData[uiSiteIndex] = _makeZero<deviceDataBoson>();
    }
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceBoson(deviceDataBoson* pDeviceData, UINT uiDesiredSite, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _makeColorVector<deviceDataBoson>(byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = _makeZero<deviceDataBoson>();
    }
}

template<typename deviceDataBoson, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonOneLink(
    const deviceDataBoson* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    deviceDataBoson* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    _deviceCoeffFunctionPointerTwoSites pfcoeff,
    const INT* __restrict__ path,
    BYTE pathLength,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = _makeZero<deviceDataBoson>();
        return;
    }

    INT pathBuffer[_kLinkMaxLength];
    deviceDataBoson result = _makeZero<deviceDataBoson>();
    SSmallInt4 siten = _deviceSmallInt4OffsetC(sSite4, path, pathLength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    deviceGauge vn = _makeId<deviceGauge>();
    if (!sn1.IsDirichlet())
    {
        const Real fCoefficient1 = (*pfcoeff)(byFieldId, sSite4, siten, sIdx, sn1) * fCoefficient;
        if (NULL != pGauge)
        {
            vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, path);
        }
        _add(result, _mulC(_mulVec(vn, pDeviceData[sn1.m_uiSiteIndex]), fCoefficient1));
    }

    _devicePathDagger(path, pathBuffer, pathLength);
    siten = _deviceSmallInt4OffsetC(sSite4, pathBuffer, pathLength);
    const UINT uiBigIndexSiten = __bi(siten);
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndexSiten];
    if (!sn2.IsDirichlet())
    {
        const Real fCoefficient2 = (*pfcoeff)(byFieldId, siten, sSite4, sn2, sIdx) * fCoefficient;
        if (NULL != pGauge)
        {
            vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, pathBuffer);
        }
        _add(result, _mulC(_mulVec(vn, pDeviceData[sn2.m_uiSiteIndex]), fCoefficient2));
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceDataBoson, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonForce_WithLink(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceDataBoson* __restrict__ pBoson,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    _deviceCoeffFunctionPointerTwoSites pfcoeff,
    const INT* __restrict__ path,
    BYTE pathLength)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    INT pathLeft[_kLinkMaxLength];
    INT pathRight[_kLinkMaxLength];
    

    for (BYTE iSeperation = 0; iSeperation <= pathLength; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, pathLength, pathLeft, pathRight, LLength, RLength);

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
            if (!sn1.IsDirichlet() && !sn2.IsDirichlet())
            {
                const Real fCoeff = (*pfcoeff)(byFieldId, siten1, siten2, sn1, sn2) * fCoefficient;

                //=================================
                // 2. Find V(n,n1), V(n,n2)
                const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
                const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

                deviceDataBoson phi1 = _mulVec(vnn1, pBoson[sn1.m_uiSiteIndex]);
                deviceDataBoson phi2 = _mulVec(vnn2, pBoson[sn2.m_uiSiteIndex]);

                deviceGauge res = _makeContract<deviceGauge, deviceDataBoson>(phi1, phi2);
                _ta(res);
                _mul(res, fCoeff);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    _add(pForce[linkIndex], res);
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    _sub(pForce[linkIndex], res);
                }
            }
        }
    }
}

/**
* f(x) phi^*(x)phi(x)
*/
template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonDiagnal(
    const deviceDataBoson* __restrict__ pDeviceData,
    deviceDataBoson* pResultData,
    BYTE byFieldId,
    Real fCoefficient,
    _deviceCoeffFunctionPointer pfcoeff,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = _makeZero<deviceDataBoson>();
        return;
    }

    const Real fCoefficient1 = (*pfcoeff)(byFieldId, sSite4, sIdx) * fCoefficient;
    deviceDataBoson result = _mulC(pDeviceData[uiSiteIndex], fCoefficient1);
    
    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}


/**
* f(n) U_mu(n)phi(n+mu) + f(n-mu) U_{-mu}(n)phi(n-mu)
*/
template<typename deviceDataBoson, typename deviceDataGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDPartialBosonVN(
    const deviceDataBoson* __restrict__ pDeviceData,
    const deviceDataGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pBosonMove,
    deviceDataBoson* pResultData,
    BYTE idir,
    Real fCoefficient,
    _deviceCoeffFunctionPointer pfcoeff,
    EOperatorCoefficientType eCoeff,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = _makeZero<deviceDataBoson>();
        return;
    }

    deviceDataBoson result = _makeZero<deviceDataBoson>();

    //idir = mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Boson = pBosonMove[2 * linkIndex];
    const SIndex& x_m_mu_Boson = pBosonMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    if (!x_p_mu_Boson.IsDirichlet())
    {
        const deviceDataGauge& x_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[linkIndex];
        //U(x,mu) phi(x+ mu)
        const deviceDataBoson u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Boson.m_uiSiteIndex]);
        const Real fCoefficient1 = (*pfcoeff)(byFieldId, sSite4, sIdx) * fCoefficient;
        _add(result, _mulC(u_phi_x_p_m, fCoefficient1));
    }

    if (!x_m_mu_Boson.IsDirichlet())
    {
        deviceDataGauge x_m_mu_Gauge_element = NULL == pGauge ? _makeId<deviceDataGauge>() : pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        const deviceDataBoson u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Boson.m_uiSiteIndex]);
        sSite4.m_byData4[idir] = sSite4.m_byData4[idir] - 1;
        const Real fCoefficient2 = (*pfcoeff)(byFieldId, sSite4, x_m_mu_Boson) * fCoefficient;
        _add(result, _mulC(u_dagger_phi_x_m_m, fCoefficient2));
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fRealCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCmpCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceDataBoson, typename deviceDataGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDPartialBosonForceVN(
    const deviceDataBoson* __restrict__ pBoson,
    const deviceDataGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pBosonMove,
    deviceDataGauge* pGaugeForce,
    BYTE idir,
    Real fCoefficient,
    _deviceCoeffFunctionPointer pfcoeff,
    BYTE byGaugeFieldId,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const deviceDataBoson& phi_dagger = pBoson[uiSiteIndex];

    //idir = mu
    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    const SIndex& x_p_mu_Boson = pBosonMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));
    if (x_p_mu_Boson.IsDirichlet())
    {
        return;
    }

    const deviceDataBoson& phi_p_mu = pBoson[x_p_mu_Boson.m_uiSiteIndex];

    const deviceDataGauge& x_Gauge_element = pGauge[linkIndex];

    //U phi(n+mu)phi^+(n)
    deviceDataGauge forceOfThisLink = _makeContract<deviceDataGauge, deviceDataBoson>(phi_dagger, _mulVec<deviceDataGauge, deviceDataBoson>(x_Gauge_element, phi_p_mu));
    const Real fCoefficient1 = (*pfcoeff)(byFieldId, sSite4, sIdx) * fCoefficient;
    _mul(forceOfThisLink, fCoefficient1);
    //TA
    //pForce[linkIndex] = _cuCsubf(pForce[linkIndex], _make_cuComplex(F(0.0), forceOfThisLink.x));
    //pGaugeForce[linkIndex].y = pGaugeForce[linkIndex].y - F(1.0) * forceOfThisLink.y;
    _ta<deviceDataGauge>(forceOfThisLink);
    _sub(pGaugeForce[linkIndex], forceOfThisLink);
}

#pragma endregion

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);

    if (NULL != pGauge && VectorN() != pGauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonVU can only play with gauge UN!"));
        return;
    }

    if (NULL == pSource || GetFieldType() != pSource->GetFieldType())
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
        return;
    }

    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pSourceVN = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(pSource);
    if (NULL == pSourceVN)
    {
        appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
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

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    _kernelDBosonVN << <block, threads >> > (
        pSourceVN->m_pDeviceData,
        NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        m_pDeviceData,
        eCoeffType,
        fRealCoeff,
        cCompCoeff,
        m_byFieldId);

    checkCudaErrors(cudaDeviceSynchronize());
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

    if (NULL == gauge || VectorN() != gauge->MatrixN())
    {
        appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
        return;
    }

    if (NULL == gaugeforce || VectorN() != gaugeforce->MatrixN())
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
TArray<DOUBLE> CFieldBosonVN<deviceDataBoson, deviceDataGauge>::Sum() const
{
    preparethread;
    TArray<DOUBLE> ret;
    for (UINT i = 0; i < FloatN(); ++i)
    {
        _kernelElementVN << <block, threads >> > (m_pDeviceData, i, _D_RealThreadBuffer);
        ret.AddItem(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));
    }
    return ret;
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
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::FixBoundary()
{
    //appGeneral(_T("CFieldBosonVN<deviceDataBoson, deviceDataGauge>::FixBoundary()\n"));
    preparethread;
    _kernelBosonFixBoundary << <block, threads >> > (m_pDeviceData, m_byFieldId);
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

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DiagnalTerm(
    const deviceDataBoson* pSource, Real fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff,
    EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelBosonDiagnal << <block, threads >> > (
        pSource,
        m_pDeviceData,
        m_byFieldId,
        fCoeffiecient,
        fpCoeff,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::OneLink(
    const deviceDataBoson* pSource,
    const deviceDataGauge* pGuage,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    _deviceCoeffFunctionPointerTwoSites fpCoeff,
    const INT* pDevicePath,
    BYTE pathLength,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    const CLGComplex& cCmpCoeff)
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelBosonOneLink << <block, threads >> > (
        pSource,
        (const deviceDataGauge*)pGuage,
        m_pDeviceData,
        m_byFieldId,
        byGaugeFieldId,
        fCoefficient,
        fpCoeff,
        pDevicePath,
        pathLength,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::OneLinkForceGauge(
    const deviceDataGauge* pGuage, 
    BYTE byGaugeFieldId, 
    deviceDataGauge* pForce, 
    Real fCoefficient,
    _deviceCoeffFunctionPointerTwoSites fpCoeff,
    const INT* pDevicePath, 
    BYTE pathLength) const
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelBosonForce_WithLink << <block, threads >> > (
        (const deviceDataGauge*)pGuage,
        (deviceDataGauge*)pForce,
        m_pDeviceData,
        m_byFieldId,
        byGaugeFieldId,
        fCoefficient,
        fpCoeff,
        pDevicePath,
        pathLength
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::PartialSq(
    const deviceDataBoson* pSource,
    const deviceDataGauge* pGuage,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    _deviceCoeffFunctionPointer fpCoeff,
    BYTE idir,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDPartialBosonVN << <block, threads >> > (
        pSource,
        pGuage,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        m_pDeviceData,
        idir,
        fCoefficient,
        fpCoeff,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        m_byFieldId
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::PartialSqForceGauge(
    const deviceDataGauge* pGuage,
    BYTE byGaugeFieldId,
    deviceDataGauge* pForce,
    Real fCoefficient,
    _deviceCoeffFunctionPointer fpCoeff,
    BYTE idir) const
{
    preparethread;
    _kernelDPartialBosonForceVN << <block, threads >> > (
        m_pDeviceData,
        pGuage,
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        pForce,
        idir,
        fCoefficient,
        fpCoeff,
        byGaugeFieldId,
        m_byFieldId
        );
}

//template<typename deviceDataBoson, typename deviceDataGauge>
//void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::OneLinkForceBoson(const deviceDataGauge* pGuage, BYTE byGaugeFieldId, deviceDataBoson* pForce, _deviceCoeffFunctionPointer fpCoeff, const INT* pDevicePath, BYTE pathLength) const
//{
//
//}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialAsSource(const SFermionBosonSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
        {
            preparethread;
            _kernelMakePointSourceBoson << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_byColorIndex);
        }
        break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionBosonSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

template<typename deviceDataBoson, typename deviceDataGauge>
UINT CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CheckHermitian(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) const
{
    const UINT uiVolume = m_uiSiteCount;
    const UINT uiRealVolume = VectorN() * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceDataBoson* hostData = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * uiVolume);
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>* v = dynamic_cast<CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < VectorN(); ++j)
        {
            SFermionBosonSource source;
            source.m_byColorIndex = static_cast<BYTE>(j);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->D(gaugeNum, bosonNum, gaugeFields, pBoson);

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceDataBoson) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * VectorN() + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                for (UINT l = 0; l < VectorN(); ++l)
                {
                    matrixElement[(VectorN() * k + l) * uiRealVolume + x] = _vn(hostData[k], l);
                }
            }
            appGeneral(_T("%d / %d have been done\n"), x, uiRealVolume);
        }
    }

    UINT uiE = 0;
    UINT uiWrong = 0;
    //List all results
    for (UINT i = 0; i < uiRealVolume * uiRealVolume; ++i)
    {
        const UINT x = i / uiRealVolume;
        const UINT y = i % uiRealVolume;
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / VectorN());
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / VectorN());
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = static_cast<BYTE>(x % VectorN());
        const BYTE cy = static_cast<BYTE>(y % VectorN());

        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
        {
            ++uiE;
            if (appAbs(matrixElement[i].x - matrixElement[daggerIdx].x) > F(0.0000001)
                || appAbs(matrixElement[i].y + matrixElement[daggerIdx].y) > F(0.0000001))
            {
                ++uiWrong;
                appGeneral(_T("[(%d, %d, %d, %d)_(%d)-(%d, %d, %d, %d)_(%d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
                    xSite.x, xSite.y, xSite.z, xSite.w, cx,
                    ySite.x, ySite.y, ySite.z, ySite.w, cy,
                    matrixElement[i].x, matrixElement[i].y,
                    matrixElement[daggerIdx].x, matrixElement[daggerIdx].y);
            }
        }
    }
    v->Return();
    appSafeFree(matrixElement);
    appSafeFree(hostData);
    appGeneral(_T("%d none zero element checked, %d wrong found...\n"), uiE, uiWrong);
    return uiWrong;
}

__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, U1, CLGComplex, CLGComplex)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU2, deviceSU2Vector, deviceSU2)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU3, deviceSU3Vector, deviceSU3)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU4, deviceSU4Vector, deviceSU4)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU5, deviceSU5Vector, deviceSU5)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU6, deviceSU6Vector, deviceSU6)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU7, deviceSU7Vector, deviceSU7)
__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldBosonVN, SU8, deviceSU8Vector, deviceSU8)

__CLGIMPLEMENT_CLASS(CFieldBosonU1)
__CLGIMPLEMENT_CLASS(CFieldBosonSU2)
__CLGIMPLEMENT_CLASS(CFieldBosonSU3)
__CLGIMPLEMENT_CLASS(CFieldBosonSU4)
__CLGIMPLEMENT_CLASS(CFieldBosonSU5)
__CLGIMPLEMENT_CLASS(CFieldBosonSU6)
__CLGIMPLEMENT_CLASS(CFieldBosonSU7)
__CLGIMPLEMENT_CLASS(CFieldBosonSU8)

#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================