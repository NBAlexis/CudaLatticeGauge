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
#include "CFieldBosonVN.h"

__BEGIN_NAMESPACE

#pragma region kernel

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
UINT CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CheckHermitian(const CFieldBoson* data, UINT uiSiteCount, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    const UINT uiVolume = uiSiteCount;
    const UINT uiRealVolume = data->VectorN() * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceDataBoson* hostData = (deviceDataBoson*)malloc(sizeof(deviceDataBoson) * uiVolume);
    CFieldBosonVN<deviceDataBoson, deviceDataGauge>* v = dynamic_cast<CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(appGetLattice()->GetPooledFieldById(data->m_byFieldId));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < data->VectorN(); ++j)
        {
            SFermionBosonSource source;
            source.m_byColorIndex = static_cast<BYTE>(j);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->D(gaugeNum, bosonNum, gaugeFields, pBoson);

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceDataBoson) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * data->VectorN() + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                for (UINT l = 0; l < data->VectorN(); ++l)
                {
                    matrixElement[(data->VectorN() * k + l) * uiRealVolume + x] = _vn(hostData[k], l);
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
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / data->VectorN());
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / data->VectorN());
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = static_cast<BYTE>(x % data->VectorN());
        const BYTE cy = static_cast<BYTE>(y % data->VectorN());

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

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::ForceOnGauge(const deviceDataBoson* data, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, deviceDataGauge* force)
{
    preparethread;
    _kernelDBosonForceVN << <block, threads >> > (
        data,
        gaugedata,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        force,
        byGaugeFieldId,
        byFieldId
        );

    //gaugeu1->DebugPrintMe();
    //checkCudaErrors(cudaDeviceSynchronize());
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::DFromSource(const deviceDataBoson* source, deviceDataBoson* target, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, 
    EOperatorCoefficientType eCoeffType, Real fRealCoeff, const CLGComplex& cCompCoeff)
{
    preparethread;
    _kernelDBosonVN << <block, threads >> > (
        source,
        gaugedata,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        target,
        eCoeffType,
        fRealCoeff,
        cCompCoeff,
        byFieldId);
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::DiagnalTerm(
    deviceDataBoson* pTarget,
    BYTE byFieldId,
    const deviceDataBoson* pSource, Real fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff,
    EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelBosonDiagnal << <block, threads >> > (
        pSource,
        pTarget,
        byFieldId,
        fCoeffiecient,
        fpCoeff,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLink(
    deviceDataBoson* pTarget,
    BYTE byFieldId,
    const deviceDataBoson* pSource,
    const deviceDataGauge* pGauge,
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
        pGauge,
        pTarget,
        byFieldId,
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
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::OneLinkForceGauge(
    const deviceDataBoson* pBoson,
    BYTE byFieldId,
    const deviceDataGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceDataGauge* pForce,
    Real fCoefficient,
    _deviceCoeffFunctionPointerTwoSites fpCoeff,
    const INT* pDevicePath,
    BYTE pathLength)
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelBosonForce_WithLink << <block, threads >> > (
        pGauge,
        pForce,
        pBoson,
        byFieldId,
        byGaugeFieldId,
        fCoefficient,
        fpCoeff,
        pDevicePath,
        pathLength
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSq(
    deviceDataBoson* pTarget,
    BYTE byFieldId,
    const deviceDataBoson* pSource,
    const deviceDataGauge* pGauge,
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
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        pTarget,
        idir,
        fCoefficient,
        fpCoeff,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        byFieldId
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::PartialSqForceGauge(
    const deviceDataBoson* pBoson,
    BYTE byFieldId,
    const deviceDataGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceDataGauge* pForce,
    Real fCoefficient,
    _deviceCoeffFunctionPointer fpCoeff,
    BYTE idir)
{
    preparethread;
    _kernelDPartialBosonForceVN << <block, threads >> > (
        pBoson,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        pForce,
        idir,
        fCoefficient,
        fpCoeff,
        byGaugeFieldId,
        byFieldId
        );
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::AllocatePathBuffer(INT** pathbuffer)
{
    checkCudaErrors(cudaMalloc((void**)pathbuffer, sizeof(INT) * _kLinkMaxLength));
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::FreePathBuffer(INT* pathbuffer)
{
    checkCudaErrors(cudaFree(pathbuffer));
}

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(INT* devicepathbuffer, const INT* hostpathbuffer, BYTE length)
{
    checkCudaErrors(cudaMemcpy(devicepathbuffer, hostpathbuffer, sizeof(INT) * length, cudaMemcpyHostToDevice));
}

#pragma region rotation

__device__ Real _deviceBosonRotationCx(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx;
    }

    return F(0.5) * (fX1 + fX2);
}

__device__ Real _deviceBosonRotationCy(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery;
    }

    return F(0.5) * (fY1 + fY2);
}

__device__ Real _deviceBosonRotationCxy(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx;
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery;
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery;
    }

    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxShift(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    return F(0.5) * (fX1 + fX2);
}

__device__ Real _deviceBosonRotationCyShift(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    return F(0.5) * (fY1 + fY2);
}

/**
* first x then y
*/
__device__ Real _deviceBosonRotationCxyShiftXYPP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lxi == site2.x || -1 == site2.y || _DC_Lxi == site1.x || -1 == site1.y)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftYXPP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lyi == site2.y || -1 == site2.x || _DC_Lyi == site1.y || -1 == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftYXPM(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (_DC_Lyi == site2.y || _DC_Lxi == site2.x || _DC_Lyi == site1.y || _DC_Lxi == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxyShiftXYMP(BYTE byFieldId, SSmallInt4 site1, SSmallInt4 site2, const SIndex& uiSiteBI1, const SIndex& uiSiteBI2)
{
    UBOOL bOppsite = FALSE;
    if (-1 == site2.y || -1 == site2.x || -1 == site1.y || -1 == site1.x)
    {
        bOppsite = TRUE;
    }

    Real fX1, fX2;
    if (uiSiteBI1.IsDirichlet())
    {
        fX1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fX1 = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fX2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fX2 = static_cast<Real>(site2.x) - _DC_Centerx + F(0.5);
    }

    Real fY1, fY2;
    if (uiSiteBI1.IsDirichlet())
    {
        fY1 = F(0.0);
    }
    else
    {
        if (site1.Out())
        {
            site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
        }
        fY1 = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    }

    if (uiSiteBI2.IsDirichlet())
    {
        fY2 = F(0.0);
    }
    else
    {
        if (site2.Out())
        {
            site2 = __deviceSiteIndexToInt4(uiSiteBI2.m_uiSiteIndex);
        }
        fY2 = static_cast<Real>(site2.y) - _DC_Centery + F(0.5);
    }

    if (bOppsite)
    {
        return F(-0.25) * ((fX1 + fX2) * (fY1 + fY2));
    }
    return F(0.25) * ((fX1 + fX2) * (fY1 + fY2));
}

__device__ Real _deviceBosonRotationCxCySq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx;
    const Real y = static_cast<Real>(site1.y) - _DC_Centery;
    return x * x + y * y;
}

__device__ Real _deviceBosonRotationCxCySqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    const Real y = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    return x * x + y * y;
}


__device__ Real _deviceBosonRotationCxSq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx;
    return x * x;
}

__device__ Real _deviceBosonRotationCySq(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real y = static_cast<Real>(site1.y) - _DC_Centery;
    return y * y;
}

__device__ Real _deviceBosonRotationCxSqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real x = static_cast<Real>(site1.x) - _DC_Centerx + F(0.5);
    return x * x;
}

__device__ Real _deviceBosonRotationCySqShift(BYTE byFieldId, SSmallInt4 site1, const SIndex& uiSiteBI1)
{
    if (site1.Out())
    {
        site1 = __deviceSiteIndexToInt4(uiSiteBI1.m_uiSiteIndex);
    }
    const Real y = static_cast<Real>(site1.y) - _DC_Centery + F(0.5);
    return y * y;
}

__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxPf = _deviceBosonRotationCx;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCyPf = _deviceBosonRotationCy;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPf = _deviceBosonRotationCxy;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxPfShift = _deviceBosonRotationCxShift;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCyPfShift = _deviceBosonRotationCyShift;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftXYPP = _deviceBosonRotationCxyShiftXYPP;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftYXPP = _deviceBosonRotationCxyShiftYXPP;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftYXPM = _deviceBosonRotationCxyShiftYXPM;
__device__ _deviceCoeffFunctionPointerTwoSites _deviceBosonRotationCxyPfShiftXYMP = _deviceBosonRotationCxyShiftXYMP;

__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxCySqPf = _deviceBosonRotationCxCySq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxCySqShiftPf = _deviceBosonRotationCxCySqShift;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxSqPf = _deviceBosonRotationCxSq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCySqPf = _deviceBosonRotationCySq;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCxSqShiftPf = _deviceBosonRotationCxSqShift;
__device__ _deviceCoeffFunctionPointer _deviceBosonRotationCySqShiftPf = _deviceBosonRotationCySqShift;

template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCx(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCy(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCyPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxy(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxyPf, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxShift(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxPfShift, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCyShift(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCyPfShift, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftXYPP(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxyPfShiftXYPP, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftYXPP(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxyPfShiftYXPP, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftYXPM(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxyPfShiftYXPM, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftXYMP(_deviceCoeffFunctionPointerTwoSites* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxyPfShiftXYMP, sizeof(_deviceCoeffFunctionPointerTwoSites)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxCySq(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxCySqPf, sizeof(_deviceCoeffFunctionPointer)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxCySqShift(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxCySqShiftPf, sizeof(_deviceCoeffFunctionPointer)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxSq(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxSqPf, sizeof(_deviceCoeffFunctionPointer)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCySq(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCySqPf, sizeof(_deviceCoeffFunctionPointer)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxSqShift(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCxSqShiftPf, sizeof(_deviceCoeffFunctionPointer)));
}
template<typename deviceDataBoson, typename deviceDataGauge>
void CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCySqShift(_deviceCoeffFunctionPointer* target)
{
    checkCudaErrors(cudaMemcpyFromSymbol(target, _deviceBosonRotationCySqShiftPf, sizeof(_deviceCoeffFunctionPointer)));
}

#pragma endregion

template class CFieldBosonVNKernel<CLGComplex, CLGComplex>;
template class CFieldBosonVNKernel<deviceSU2Vector, deviceSU2>;
template class CFieldBosonVNKernel<deviceSU3Vector, deviceSU3>;
template class CFieldBosonVNKernel<deviceSU4Vector, deviceSU4>;
template class CFieldBosonVNKernel<deviceSU5Vector, deviceSU5>;
template class CFieldBosonVNKernel<deviceSU6Vector, deviceSU6>;
template class CFieldBosonVNKernel<deviceSU7Vector, deviceSU7>;
template class CFieldBosonVNKernel<deviceSU8Vector, deviceSU8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================