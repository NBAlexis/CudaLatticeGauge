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
//  [07/21/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKST.h"

//save registor strategy uses global memory frequently, so is slower
//however, it reduces the number of registers, enabling more threads at the one time
#define _CLG_SAVE_REGISTOR_STRATEGY 0

__BEGIN_NAMESPACE

#define _AddEachDir \
if (0 == dir) \
{ \
    if (bDDagger) \
    { \
        _oppo(result); \
    } \
    pDeviceData[uiSiteIndex] = result; \
} \
__syncthreads(); \
if (1 == dir) \
{ \
    if (bDDagger) \
    { \
        _sub(pDeviceData[uiSiteIndex], result); \
    } \
    else \
    { \
        _add(pDeviceData[uiSiteIndex], result); \
    } \
} \
__syncthreads(); \
if (2 == dir) \
{ \
    if (bDDagger) \
    { \
        _sub(pDeviceData[uiSiteIndex], result); \
    } \
    else \
    { \
        _add(pDeviceData[uiSiteIndex], result); \
    } \
} \
__syncthreads(); \
if (3 == dir) \
{ \
    if (bDDagger) \
    { \
        _sub(pDeviceData[uiSiteIndex], result); \
    } \
    else \
    { \
        _add(pDeviceData[uiSiteIndex], result); \
    } \
} \
__syncthreads();


template<typename deviceVector, typename deviceGauge, INT vectorN>
UINT CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::TestAntiHermitianS(BYTE byFieldId, const CFieldGauge* pGauge)
{
    const UINT uiVolume = _HC_Volume;
    const UINT uiRealVolume = vectorN * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceVector* hostData = (deviceVector*)malloc(sizeof(deviceVector) * uiVolume);
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* v = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < vectorN; ++j)
        {
            SFermionBosonSource source;
            source.m_byColorIndex = static_cast<BYTE>(j);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->D0S(pGauge);

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceVector) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * vectorN + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                for (UINT kcolor = 0; kcolor < vectorN; ++kcolor)
                {
                    matrixElement[(vectorN * k + kcolor) * uiRealVolume + x] = _make_cuComplex(_element(hostData[k], 2 * kcolor), _element(hostData[k], 2 * kcolor + 1));
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
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / vectorN);
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / vectorN);
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = x % vectorN;
        const BYTE cy = y % vectorN;

        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
        {
            ++uiE;
            if (appAbs(matrixElement[i].x + matrixElement[daggerIdx].x) > F(0.0000001)
                || appAbs(matrixElement[i].y - matrixElement[daggerIdx].y) > F(0.0000001))
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

#pragma region DOperator

#pragma region kernel

/**
* Dks = 2am + \sum _{\mu} \eta_{\mu} (n) (U_{\mu}(n) \delta _{n,n+\mu} -U^+_{\mu}(n-\mu) \delta _{n,n-\mu})
* U act on su3
* gamma act on spinor
*
* If bDagger, it is just \eta 5(n) Dks \eta _5(n)
* A easier way is D_{ks}(m=0) is anti-Hermitian, so D^+ = - D_{ks,m=0} + 2am
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

//template<typename deviceVector, typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelDFermionKST_DOnEvenOrOdd(
//    deviceVector* pDeviceData,
//    const deviceGauge* __restrict__ pGauge,
//    const SIndex* __restrict__ pGaugeMove,
//    const SIndex* __restrict__ pFermionMove,
//    const BYTE* __restrict__ pEtaTable,
//    BYTE byFieldId,
//    Real f2am,
//    UBOOL bEven,
//    UBOOL bDDagger,
//    EOperatorCoefficientType eCoeff,
//    Real fCoeff,
//    CLGComplex cCoeff)
//{
//    intokernaldirEO;
//
//    deviceVector result = _makeZero<deviceVector>();
//
//    //idir = mu
//    for (UINT idir = 0; idir < uiDir; ++idir)
//    {
//        //x, mu
//        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
//
//        const SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];
//
//        const SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
//        const SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];
//
//        //Assuming periodic
//        //get U(x,mu), U^{dagger}(x-mu), 
//        const deviceGauge x_Gauge_element = pGauge[linkIndex];
//        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
//        if (x_m_mu_Gauge.NeedToDagger())
//        {
//            _dagger(x_m_mu_Gauge_element);
//        }
//
//        //U(x,mu) phi(x+ mu)
//        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
//        if (x_p_mu_Fermion.NeedToOpposite())
//        {
//            _oppo(u_phi_x_p_m);
//        }
//
//        //U^{dagger}(x-mu) phi(x-mu)
//        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
//        if (x_m_mu_Fermion.NeedToOpposite())
//        {
//            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
//        }
//        else
//        {
//            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
//        }
//
//        if ((pEtaTable[uiSiteIndex] >> idir) & 1)
//        {
//            _sub(result, u_phi_x_p_m);
//        }
//        else
//        {
//            _add(result, u_phi_x_p_m);
//        }
//    }
//
//    if (abs(f2am) > _CLG_FLT_MIN_)
//    {
//        _mul(pDeviceData[uiSiteIndex], f2am);
//        if (bDDagger)
//        {
//            _sub(pDeviceData[uiSiteIndex], result);
//        }
//        else
//        {
//            _add(pDeviceData[uiSiteIndex], result);
//        }
//    }
//    else
//    {
//        if (bDDagger)
//        {
//            _oppo(result);
//            pDeviceData[uiSiteIndex] = result;
//        }
//        else
//        {
//            pDeviceData[uiSiteIndex] = result;
//        }
//    }
//
//    switch (eCoeff)
//    {
//    case EOCT_Real:
//        _mul(pDeviceData[uiSiteIndex], fCoeff);
//        break;
//    case EOCT_Complex:
//        _mul(pDeviceData[uiSiteIndex], cCoeff);
//        break;
//    }
//}


template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    deviceVector result = _makeZero<deviceVector>();
    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        result = pDeviceData[uiSiteIndex];
        _mul(result, f2am);
    }

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
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
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOddDir(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalDirEO;

    const SIndex x_m_mu_Gauge = pGaugeMove[uiLinkIndex];
    const UINT doubleLinkIndex = (uiLinkIndex << 1);
    const SIndex x_p_mu_Fermion = pFermionMove[doubleLinkIndex];
    const SIndex x_m_mu_Fermion = pFermionMove[doubleLinkIndex | 1];

    //const SIndex x_p_mu_Fermion = pFermionMove[2 * uiLinkIndex];
    //const SIndex x_m_mu_Fermion = pFermionMove[2 * uiLinkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //const deviceGauge x_Gauge_element = pGauge[uiLinkIndex];
    //deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, dir)];
    //not consider projective plane here
    //if (x_m_mu_Gauge.NeedToDagger())
    //{
       // _dagger(x_m_mu_Gauge_element);
    //}

    //U(x,mu) phi(x+ mu)
    deviceVector result = _mulVec(pGauge[uiLinkIndex], pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
    const BYTE etamu = (pEtaTable[uiSiteIndex] >> dir) & 1U;
    if (_UBOOLXOR(x_p_mu_Fermion.NeedToOpposite(), etamu))
    {
        _oppo(result);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceVector u_dagger_phi_x_m_m = _dagmulVec(pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, dir)], pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_m_mu_Fermion.NeedToOpposite(), etamu))
    {
        _add(result, u_dagger_phi_x_m_m);
    }
    else
    {
        _sub(result, u_dagger_phi_x_m_m);
    }

    //__syncthreads();
    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        if (0 == dir)
        {
            _mul(pDeviceData[uiSiteIndex], f2am);
            if (bDDagger)
            {
                _sub(pDeviceData[uiSiteIndex], result);
            }
            else
            {
                _add(pDeviceData[uiSiteIndex], result);
            }
        }
        __syncthreads();
        if (1 == dir)
        {
            if (bDDagger)
            {
                _sub(pDeviceData[uiSiteIndex], result);
            }
            else
            {
                _add(pDeviceData[uiSiteIndex], result);
            }
        }
        __syncthreads();
        if (2 == dir)
        {
            if (bDDagger)
            {
                _sub(pDeviceData[uiSiteIndex], result);
            }
            else
            {
                _add(pDeviceData[uiSiteIndex], result);
            }
        }
        __syncthreads();
        if (3 == dir)
        {
            if (bDDagger)
            {
                _sub(pDeviceData[uiSiteIndex], result);
            }
            else
            {
                _add(pDeviceData[uiSiteIndex], result);
            }
        }
        __syncthreads();
    }
    else
    {
        _AddEachDir;
    }

    if (0 == dir)
    {
        switch (eCoeff)
        {
        case EOCT_Real:
            _mul(pDeviceData[uiSiteIndex], fCoeff);
            break;
        case EOCT_Complex:
            _mul(pDeviceData[uiSiteIndex], cCoeff);
            break;
        default:
            break;
        }
    }
}



template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOddDirPM(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    const UINT uiLinkIndex = (threadIdx.x + blockIdx.x * blockDim.x) >> 1U;
    const UINT uiSiteIndex = uiLinkIndex >> 2U;
    if (uiSiteIndex >= _DC_Volume)
    {
        return;
    }
    if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven))
    {
        return;
    }
    const BYTE pm = threadIdx.x & 1U;
    const BYTE term = threadIdx.x & 7U;
    const BYTE dir = uiLinkIndex & 3U;
    const SIndex x_move_Fermion = pFermionMove[(uiLinkIndex << 1U) | pm];

    deviceVector result = pm ? _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, dir)], pDeviceData[x_move_Fermion.m_uiSiteIndex])
        : _mulVec(pGauge[uiLinkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR4(x_move_Fermion.NeedToOpposite(), (pEtaTable[uiSiteIndex] >> dir) & 1U, bDDagger, pm))
    {
        _oppo(result);
    }

    if (0 == term)
    {
        if (abs(f2am) > _CLG_FLT_MIN_)
        {
            _mul(pDeviceData[uiSiteIndex], f2am);
            _add(pDeviceData[uiSiteIndex], result);
        }
        else
        {
            pDeviceData[uiSiteIndex] = result;
        }
    }
    __syncthreads();
    #pragma unroll
    for (BYTE i = 1U; i < 8U; ++i)
    {
        if (i == term)
        {
            _add(pDeviceData[uiSiteIndex], result);
        }
        __syncthreads();
    }

    if (0 == term)
    {
        switch (eCoeff)
        {
        case EOCT_Real:
            _mul(pDeviceData[uiSiteIndex], fCoeff);
            break;
        case EOCT_Complex:
            _mul(pDeviceData[uiSiteIndex], cCoeff);
            break;
        default:
            break;
        }
    }
}

//template<typename deviceVector, typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelDFermionKST_DOnEvenOrOddDirPM2(
//    deviceVector* pDeviceData,
//    const deviceGauge* __restrict__ pGauge,
//    const SIndex* __restrict__ pFermionMove,
//    const BYTE* __restrict__ pEtaTable,
//    BYTE byFieldId,
//    Real f2am,
//    UBOOL bEven, //make sure bEven is 1 or 0 before send in
//    UBOOL bDDagger,
//    EOperatorCoefficientType eCoeff,
//    Real fCoeff,
//    CLGComplex cCoeff)
//{
//    extern __shared__ BYTE _share_KST_DOnEvenOrOddDirPM2[];
//    //bEven means skip even, so, we keep odd sites
//    //which is, translate the half site to 2*s + 1
//    UINT uiSiteIndex = (((threadIdx.x + blockIdx.x * blockDim.x) >> 3U) << 1U) | bEven;
//    if (uiSiteIndex >= _DC_Volume)
//    {
//        return;
//    }
//    // x o
//    // o x
//    // If bEven, we need uiSiteIndex - 1
//    // If not even, we need uiSiteIndex + 1
//    // So uiSiteIndex + 1 - 2 * bEven
//    BYTE eta = pEtaTable[uiSiteIndex];
//    if (_UBOOLXOR((eta >> 4U) & 1U, bEven))
//    {
//        uiSiteIndex = uiSiteIndex + 1U - (bEven << 1U);
//        eta = pEtaTable[uiSiteIndex];
//    }
//    const BYTE term = threadIdx.x & 7U;
//    const BYTE pm = term & 1U;
//    const BYTE dir = term >> 1U;
//    const UINT uiLinkIndex = (uiSiteIndex << 2U) | dir;
//    const SIndex x_move_Fermion = pFermionMove[(uiLinkIndex << 1U) | pm];
//
//    deviceVector result = pm ? _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, dir)], pDeviceData[x_move_Fermion.m_uiSiteIndex])
//        : _mulVec(pGauge[uiLinkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
//    if (_UBOOLXOR4(x_move_Fermion.NeedToOpposite(), (eta >> dir) & 1U, bDDagger, pm))
//    {
//        _oppo(result);
//    }
//
//    deviceVector* share_KST_DOnEvenOrOddDirPM2 = reinterpret_cast<deviceVector*>(_share_KST_DOnEvenOrOddDirPM2);
//    const UINT sharedidx = threadIdx.x >> 3U;
//    if (0 == term)
//    {
//        if (abs(f2am) > _CLG_FLT_MIN_)
//        {
//            share_KST_DOnEvenOrOddDirPM2[sharedidx] = pDeviceData[uiSiteIndex];
//            _mul(share_KST_DOnEvenOrOddDirPM2[sharedidx], f2am);
//            _add(share_KST_DOnEvenOrOddDirPM2[sharedidx], result);
//        }
//        else
//        {
//            share_KST_DOnEvenOrOddDirPM2[sharedidx] = result;
//        }
//    }
//    __syncthreads();
//    #pragma unroll
//    for (BYTE i = 1U; i < 8U; ++i)
//    {
//        if (i == term)
//        {
//            _add(share_KST_DOnEvenOrOddDirPM2[sharedidx], result);
//        }
//        __syncthreads();
//    }
//
//    if (0 == term)
//    {
//        switch (eCoeff)
//        {
//        case EOCT_Real:
//            _mul(share_KST_DOnEvenOrOddDirPM2[sharedidx], fCoeff);
//            break;
//        case EOCT_Complex:
//            _mul(share_KST_DOnEvenOrOddDirPM2[sharedidx], cCoeff);
//            break;
//        default:
//            break;
//        }
//
//        pDeviceData[uiSiteIndex] = share_KST_DOnEvenOrOddDirPM2[sharedidx];
//    }
//}

#pragma region e-o D

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_MCR(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoeff)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }
    _mul(result, fCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_MCC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }
    _mul(result, cCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_M(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

    #pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }

    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_CR(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoeff)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }

    _mul(result, fCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_CC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }

    _mul(result, cCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_NOMC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _adddagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
        else
        {
#if _CLG_SAVE_REGISTOR_STRATEGY
            _subdagmulP(&result, pGauge + _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
            _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
        }
    }

    pDeviceData[uiSiteIndex] = result;
}

#pragma endregion
/**
 * For some strange boundary condition
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0) * eta_mu);
        }
        else
        {
            _mul(u_phi_x_p_m, eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        _mul(u_dagger_phi_x_m_m, eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        //_mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

/**
 * Calculate Force
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceT(
    //const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId)
{
    intokernaldir;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
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
                _mul(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR]);
            }

            _sub(pForce[linkIndex], thisTerm);
        }
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    if (bEachSiteEta)
    {
        _LAUNCH_KERNEL(_kernelDFermionKSPlusEtaT TMPARG(deviceVector, deviceGauge), block, threads,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelDFermionKST TMPARG(deviceVector, deviceGauge), block, threads,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd(deviceVector* pTargetBuffer,
    const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEvenOdd, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    if (bEvenOdd && abs(f2am) > _CLG_FLT_MIN_)
    {
        appGeneral(_T("Usually D operator on even usually does not have mass, but call with mass=%f\n"), f2am);
    }

    //preparethreadDir;
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOddDir, block, threads, 
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd,
    //    bDagger,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);

    //very slow
    //preparethreadE(8);
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOddDirPM, block, threads, 
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd ? 1 : 0,
    //    bDagger ? 1 : 0,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);


    preparethreadHalf;
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2, block, threads, 
    //        pTargetBuffer,
    //        pGaugeBuffer,
    //        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //        appGetLattice()->m_pIndexCache->m_pEtaMu,
    //        byFieldId,
    //        f2am,
    //        bEvenOdd ? 1 : 0,
    //        bDagger ? 1 : 0,
    //        eOCT,
    //        fRealCoeff,
    //        cCmpCoeff);

    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        switch (eOCT)
        {
        case EOCT_None:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_M TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0);
            break;
        case EOCT_Real:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_MCR TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fRealCoeff);
            break;
        case EOCT_Complex:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_MCC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                cCmpCoeff);
            break;
        default:
            appCrucial(_T("should not call this with eoc:%d\n"), eOCT);
            break;
        }
    }
    else
    {
        switch (eOCT)
        {
        case EOCT_None:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_NOMC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0);
            break;
        case EOCT_Real:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_CR TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fRealCoeff);
            break;
        case EOCT_Complex:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_CC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                cCmpCoeff);
            break;
        default:
            appCrucial(_T("should not call this with eoc:%d\n"), eOCT);
            break;
        }
    }

    //UINT block, threads; 
    //appBlockThreadsE(_HC_VolumeHalf, 8, block, threads);
    //_kernelDFermionKST_DOnEvenOrOddDirPM2 << <block, threads, sizeof(deviceVector)*(threads/8) >> > (
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd ? 1 : 0,
    //    bDagger ? 1 : 0,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);

    //if (abs(f2am) > _CLG_FLT_MIN_)
    //{
    //    if (eOCT == EOCT_None)
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_M << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            f2am,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0);
    //    }
    //    else
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_MC << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            f2am,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0,
    //            eOCT,
    //            fRealCoeff,
    //            cCmpCoeff);
    //    }
    //}
    //else
    //{
    //    if (eOCT == EOCT_None)
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_NOMC << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0);
    //    }
    //    else
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_C << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0,
    //            eOCT,
    //            fRealCoeff,
    //            cCmpCoeff);
    //    }
    //}

}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0(
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId, 
    const deviceVector* const* pRationalFields, 
    const Real* pNumerator, 
    UINT uiRationApproxOrder)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSForceT TMPARG(deviceVector, deviceGauge), block, threads,
        //pGaugeBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId);
}

#pragma endregion

#pragma region DOperator Dirichlet

#pragma region kernel

/**
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_DT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIndex = __bi(sSite4);
    deviceVector result = _makeZero<deviceVector>();

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        pResultData[uiSiteIndex] = result;
        return;
    }
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];
        deviceVector res = _makeZero<deviceVector>();

        //UBOOL btestpf = FALSE;
        if (!x_p_mu_Fermion.IsDirichlet())
        {
            const deviceGauge& x_Gauge_element = _deviceGetGaugeBCDirT(byGaugeFieldId, pGauge, uiBigIndex, idir);
            //U(x,mu) phi(x+ mu)
            res = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(res, F(-1.0));
            }
            //btestpf = TRUE;
        }
        //else
        //{
        //    printf("do we have dirichlet?\n");
        //}

        //UBOOL btestmf = FALSE;
        if (!x_m_mu_Fermion.IsDirichlet())
        {
            deviceGauge x_m_mu_Gauge_element = _deviceGetGaugeBCT(byGaugeFieldId, pGauge, x_m_mu_Gauge);
            if (x_m_mu_Gauge.NeedToDagger())
            {
                _dagger(x_m_mu_Gauge_element);
            }
            const deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
            if (x_m_mu_Fermion.NeedToOpposite())
            {
                _add(res, u_dagger_phi_x_m_m);
            }
            else
            {
                _sub(res, u_dagger_phi_x_m_m);
            }
            //btestmf = TRUE;
        }
        //else
        //{
        //    printf("do we have dirichlet?\n");
        //}
        //if (!btestpf && !btestmf)
        //{
        //    printf("both p fermion and m fermion are dirichlet?\n");
        //}

        _mul(res, eta_mu);
        _add(result, res);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    //if (__cuCabsSqf(result.m_ve[0]) + __cuCabsSqf(result.m_ve[1]) + __cuCabsSqf(result.m_ve[2]) < 1.0e-10)
    //{
    //    printf("zero result!\n");
    //}

    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

/**
 * Calculate Force
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_DT(
    //const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;
    const UINT uiBigIndex = __bi(sSite4);
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        return;
    }

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIndex, byGaugeFieldId, static_cast<BYTE>(idir)))
        {
            continue;
        }

        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        if (x_p_mu_Fermion.IsDirichlet())
        {
            continue;
        }

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceVector* phi_i = pFermionPointers[uiR];
            const deviceVector* phi_id = pFermionPointers[uiR + uiRational];
            //const deviceGauge gaugelink = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIndex, idir);
            //deviceVector toContract = _mulVec(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_i[x_p_mu_Fermion.m_uiSiteIndex], phi_id[uiSiteIndex]);

            //toContract = _mulVec(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            _sub(thisTerm, _makeContract<deviceGauge, deviceVector>(phi_id[x_p_mu_Fermion.m_uiSiteIndex], phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR]);
            }

            _sub(pForce[linkIndex], thisTerm);
        }
    }

}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_DT(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(__deviceSiteIndexToInt4(uiSiteIndex))].IsDirichlet())
    {
        //I'm Dirichlet
        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
        return;
    }

    deviceVector result = _makeZero<deviceVector>();
    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        result = pDeviceData[uiSiteIndex];
        _mul(result, f2am);
    }

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (!x_move_Fermion.IsDirichlet())
    {
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
        {
            _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (!x_move_Fermion.IsDirichlet())
    {
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
        {
            _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
    }


    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (!x_move_Fermion.IsDirichlet())
        {
            if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
            {
                _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
            }
            else
            {
                _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
            }
        }

        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (!x_move_Fermion.IsDirichlet())
        {
            if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
            {
                _add(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
            }
            else
            {
                _sub(result, _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
            }
        }
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
    pDeviceData[uiSiteIndex] = result;
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_D(UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;

    _LAUNCH_KERNEL(_kernelDFermionKS_DT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        byFieldId,
        byGaugeFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0_D(
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSForce_DT TMPARG(deviceVector, deviceGauge), block, threads,
        //pGaugeBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId,
        byGaugeFieldId);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd_D(deviceVector* pTargetBuffer,
    const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEvenOdd, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    if (bEvenOdd && abs(f2am) > _CLG_FLT_MIN_)
    {
        appGeneral(_T("Usually D operator on even usually does not have mass, but call with mass=%f\n"), f2am);
    }

    preparethreadHalf;
    _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_DT TMPARG(deviceVector, deviceGauge), block, threads,
        pTargetBuffer,
        pGaugeBuffer,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        f2am,
        bEvenOdd ? 1 : 0,
        bDagger ? 1 : 0,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

#pragma endregion


#pragma region Helper functions to implement higher orders

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLinkT(
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    BYTE byEtaIndex,
    const SCHAR* __restrict__ path,
    BYTE pathLength)
{
    intokernalInt4;
    SCHAR pathLeft[_kLinkMaxLength];
    SCHAR pathRight[_kLinkMaxLength];
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
            INT iEtaMu1 = (pEtaTable[sn1.m_uiSiteIndex] >> byEtaIndex);
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
                    _mul(res, fCoefficient * pNumerators[rfieldId]);

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
                    _mul(res, fCoefficient * pNumerators[rfieldId]);

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

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OneLinkT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    const SCHAR* __restrict__ path,
    BYTE pathLength,
    BYTE byEtaIdx,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    SCHAR pathBuffer[_kLinkMaxLength];
    deviceVector result = _makeZero<deviceVector>();

    SSmallInt4 siten = _deviceSmallInt4OffsetC(sSite4, path, pathLength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    deviceGauge vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, path);
    INT etamu = (pEtaTable[uiSiteIndex] >> byEtaIdx) & 1;
    if (sn1.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        _sub(result, _mulVec(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }

    _devicePathDagger(path, pathBuffer, pathLength);
    siten = _deviceSmallInt4OffsetC(sSite4, pathBuffer, pathLength);
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, pathBuffer);
    etamu = (pEtaTable[sn2.m_uiSiteIndex] >> byEtaIdx) & 1;
    if (sn2.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        _add(result, _mulVec(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _mulVec(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);

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

    _add(pResultData[uiSiteIndex], result);
}


template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OnlyMassT(
    const deviceVector* __restrict__ pDeviceData,
    deviceVector* pResultData,
    Real f2am,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernal;
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    _mul(pResultData[uiSiteIndex], f2am);

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OnlyMass(const deviceVector* pSource, deviceVector* pTarget, Real fm, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_OnlyMassT<deviceVector>, block, threads,
        pSource,
        pTarget,
        fm,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkS(
    const deviceVector* pSource, 
    BYTE byFieldId,
    const deviceGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceVector* pTarget,
    Real fCoefficient,
    const SCHAR* pDevicePath,
    BYTE pathLength,
    BYTE byEtaIdx,
    UBOOL bDagger,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    const CLGComplex& cCmpCoeff)
{
    appAssert(pathLength <= _kLinkMaxLength);
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_OneLinkT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fCoefficient,
        pDevicePath,
        pathLength,
        byEtaIdx,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkForceS(
    const deviceVector* pFermion, 
    BYTE byFieldId,
    const deviceGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceGauge* pForce,
    Real fCoefficient,
    const SCHAR* pDevicePath,
    BYTE pathLength,
    BYTE byEtaIdx, 
    const deviceVector* const* pRationalFields, 
    const Real* pNumerator, 
    UINT uiRationApproxOrder)
{
    appAssert(pathLength <= _kLinkMaxLength);
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSForce_WithLinkT TMPARG(deviceVector, deviceGauge), block, threads,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pForce,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId,
        byGaugeFieldId,
        fCoefficient,
        byEtaIdx,
        pDevicePath,
        pathLength
        );
}

#pragma endregion

#pragma region Matrix Operator

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
template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixMultiplyKST(
    INT vectorN,
    deviceVector** pRes,
    deviceVector** pLeft,
    const CLGComplex* __restrict__ pMatrix,
    UINT uiDimX, UINT uiDimY) //x=m,y=k
{
    intokernalE(vectorN);

    CLGComplex result[CFieldMatrixOperation::_kFieldMatrixMaxDim];

    for (UINT i = 0; i < uiDimY; ++i)
    {
        result[i] = _make_cuComplex(F(0.0), F(0.0));
        for (UINT j = 0; j < uiDimX; ++j)
        {
            const CLGComplex a = _make_cuComplex(_element(pRes[j][uiSiteIndex], 2 * elementIdx), _element(pRes[j][uiSiteIndex], 2 * elementIdx + 1));
            const CLGComplex b = _make_cuComplex(_element(pLeft[j - uiDimY][uiSiteIndex], 2 * elementIdx), _element(pLeft[j - uiDimY][uiSiteIndex], 2 * elementIdx + 1));
            result[i] = _cuCaddf(result[i], _cuCmulf(j < uiDimY ? a : b, pMatrix[j * uiDimY + i]));
        }
    }

    for (UINT i = 0; i < uiDimY; ++i)
    {
        _setelement(pRes[i][uiSiteIndex], 2 * elementIdx, result[i].x);
        _setelement(pRes[i][uiSiteIndex], 2 * elementIdx + 1, result[i].y);
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
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::VectorMultiplyMatrix(deviceVector** hostResBuffer, deviceVector** hostLeftBuffer, deviceVector** resBuffer, deviceVector** leftBuffer,
    TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY; ++i)
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(res[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        hostResBuffer[i] = pF->m_pDeviceData;
    }

    for (UINT i = 0; i < uiDimX - uiDimY; ++i)
    {
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(left[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        hostLeftBuffer[i] = pF->m_pDeviceData;
    }

    checkCudaErrors(cudaMemcpy(resBuffer, hostResBuffer, sizeof(deviceVector*) * uiDimY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(leftBuffer, hostLeftBuffer, sizeof(deviceVector*) * (uiDimX - uiDimY), cudaMemcpyHostToDevice));

    preparethreadE(vectorN);
    _LAUNCH_KERNEL(_kernelMatrixMultiplyKST<deviceVector>, block, threads, vectorN, resBuffer, leftBuffer, deviceMatrix, uiDimX, uiDimY);
}

#pragma endregion

#pragma region EM

#pragma region kernel

/**
* Dks = 2am + \sum _{\mu} \eta_{\mu} (n) (U_{\mu}(n) \delta _{n,n+\mu} -U^+_{\mu}(n-\mu) \delta _{n,n-\mu})
* U act on su3
* gamma act on spinor
*
* If bDagger, it is just \eta 5(n) Dks \eta _5(n)
* A easier way is D_{ks}(m=0) is anti-Hermitian, so D^+ = - D_{ks,m=0} + 2am
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSEMFieldT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    Real fCharge,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceGauge x_Gauge_element = pGauge[linkIndex];
        const Real fPhaseForward = pU1[linkIndex] * fCharge;
        _mul(x_Gauge_element, _make_cuComplex(_cos(fPhaseForward), _sin(fPhaseForward)));
        const UINT linkBackward = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir);
        deviceGauge x_m_mu_Gauge_element = pGauge[linkBackward];
        const Real fPhaseBackward = pU1[linkBackward] * fCharge;
        _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(fPhaseBackward), _sin(fPhaseBackward)));

        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaEMFieldEvenOddT(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    Real fCharge,
    BYTE byFieldId,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    //Even-odd version of _kernelDFermionKSPlusEtaEMFieldT
    //(the shift-center eta handling: eta_mu on the forward term, eta_mu2 on the
    //backward term, exactly as in the non even-odd PlusEta kernel)
    intokernaldirEO;

    deviceVector result = _makeZero<deviceVector>();

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceGauge x_Gauge_element = pGauge[linkIndex];
        const Real fPhaseForward = pU1[linkIndex] * fCharge;
        _mul(x_Gauge_element, _make_cuComplex(_cos(fPhaseForward), _sin(fPhaseForward)));
        const UINT linkBackward = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir);
        deviceGauge x_m_mu_Gauge_element = pGauge[linkBackward];
        const Real fPhaseBackward = pU1[linkBackward] * fCharge;
        _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(fPhaseBackward), _sin(fPhaseBackward)));

        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0) * eta_mu);
        }
        else
        {
            _mul(u_phi_x_p_m, eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        _mul(u_dagger_phi_x_m_m, eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _add(result, u_phi_x_p_m);
    }

    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        _mul(pDeviceData[uiSiteIndex], f2am);
        if (bDDagger)
        {
            _sub(pDeviceData[uiSiteIndex], result);
        }
        else
        {
            _add(pDeviceData[uiSiteIndex], result);
        }
    }
    else
    {
        if (bDDagger)
        {
            _oppo(result);
            pDeviceData[uiSiteIndex] = result;
        }
        else
        {
            pDeviceData[uiSiteIndex] = result;
        }
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pDeviceData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pDeviceData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSEMFieldEvenOddT(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    Real fCharge,
    BYTE byFieldId,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldirEO;

    deviceVector result = _makeZero<deviceVector>();

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceGauge x_Gauge_element = pGauge[linkIndex];
        const Real fPhaseForward = pU1[linkIndex] * fCharge;
        _mul(x_Gauge_element, _make_cuComplex(_cos(fPhaseForward), _sin(fPhaseForward)));
        const UINT linkBackward = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir);
        deviceGauge x_m_mu_Gauge_element = pGauge[linkBackward];
        const Real fPhaseBackward = pU1[linkBackward] * fCharge;
        _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(fPhaseBackward), _sin(fPhaseBackward)));

        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _oppo(u_phi_x_p_m);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        _mul(pDeviceData[uiSiteIndex], f2am);
        if (bDDagger)
        {
            _sub(pDeviceData[uiSiteIndex], result);
        }
        else
        {
            _add(pDeviceData[uiSiteIndex], result);
        }
    }
    else
    {
        if (bDDagger)
        {
            _oppo(result);
            pDeviceData[uiSiteIndex] = result;
        }
        else
        {
            pDeviceData[uiSiteIndex] = result;
        }
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pDeviceData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pDeviceData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}


/**
 * For some strange boundary condition
 */
 template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaEMFieldT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    Real fCharge,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceGauge x_Gauge_element = pGauge[linkIndex];
        const Real fPhaseForward = pU1[linkIndex] * fCharge;
        _mul(x_Gauge_element, _make_cuComplex(_cos(fPhaseForward), _sin(fPhaseForward)));
        const UINT linkBackward = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir);
        deviceGauge x_m_mu_Gauge_element = pGauge[linkBackward];
        const Real fPhaseBackward = pU1[linkBackward] * fCharge;
        _mul(x_m_mu_Gauge_element, _make_cuComplex(_cos(fPhaseBackward), _sin(fPhaseBackward)));

        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0) * eta_mu);
        }
        else
        {
            _mul(u_phi_x_p_m, eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        _mul(u_dagger_phi_x_m_m, eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        //u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceEMFieldT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pU1,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    Real fCharge,
    BYTE byFieldId)
{
    intokernaldir;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceVector* phi_i = pFermionPointers[uiR];
            const deviceVector* phi_id = pFermionPointers[uiR + uiRational];
            deviceGauge gauge = pGauge[linkIndex];
            Real fPhase = pU1[linkIndex] * fCharge;
            _mul(gauge, _make_cuComplex(_cos(fPhase), _sin(fPhase)));

            //deviceVector toContract = _mulVec(gauge, phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_i[x_p_mu_Fermion.m_uiSiteIndex], phi_id[uiSiteIndex]);

            //toContract = _mulVec(gauge, phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            _sub(thisTerm, _makeContract<deviceGauge, deviceVector>(phi_id[x_p_mu_Fermion.m_uiSiteIndex], phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR]);
            }

            _sub(pForce[linkIndex], thisTerm);
        }
    }

}


#pragma endregion

#pragma region Static_interface

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorEM(
    deviceVector* pTargetBuffer,
    const deviceVector* pBuffer,
    const deviceGauge* pGaugeBuffer,
    const Real* pEMFieldBuffer,
    Real f2am,
    Real fCharge,
    UBOOL bShiftCenter,
    UBOOL bDagger,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    preparethread;

    if (bShiftCenter)
    {
        _LAUNCH_KERNEL(_kernelDFermionKSPlusEtaEMFieldT TMPARG(deviceVector, deviceGauge), block, threads,
            pBuffer,
            pGaugeBuffer,
            pEMFieldBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTargetBuffer,
            f2am,
            fCharge,
            byFieldID,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelDFermionKSEMFieldT TMPARG(deviceVector, deviceGauge), block, threads,
            pBuffer,
            pGaugeBuffer,
            pEMFieldBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTargetBuffer,
            f2am,
            fCharge,
            byFieldID,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorEMEvenOdd(
    deviceVector* pTargetBuffer,
    const deviceGauge* pGaugeBuffer,
    const Real* pEMFieldBuffer,
    Real f2am,
    Real fCharge,
    UBOOL bShiftCenter,
    UBOOL bEven,
    UBOOL bDagger,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    preparethread;

    if (bShiftCenter)
    {
        //The shift-center (eta per site, i.e. rotation axis at the center) even-odd
        //kernel: eta_mu on the forward term, eta_mu2 on the backward term.
        _LAUNCH_KERNEL(_kernelDFermionKSPlusEtaEMFieldEvenOddT TMPARG(deviceVector, deviceGauge), block, threads,
            pTargetBuffer,
            pGaugeBuffer,
            pEMFieldBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            f2am,
            fCharge,
            byFieldID,
            bEven,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelDFermionKSEMFieldEvenOddT TMPARG(deviceVector, deviceGauge), block, threads,
            pTargetBuffer,
            pGaugeBuffer,
            pEMFieldBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            f2am,
            fCharge,
            byFieldID,
            bEven,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
}





template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::KSForceEM(
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    const Real* pEMFieldBuffer,
    Real fCharge,
    const deviceVector* const* pRationalFields,
    const Real* pRationalNumerator,
    UINT uiRationalDegree,
    BYTE byFieldID)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSForceEMFieldT TMPARG(deviceVector, deviceGauge), block, threads,
        pGaugeBuffer,
        pEMFieldBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pRationalNumerator,
        uiRationalDegree,
        fCharge,
        byFieldID);
}

#pragma endregion

#pragma endregion


#pragma region HISQ

#pragma region kernels

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCached(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalDir;

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * uiLinkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * uiLinkIndex + 1];

    deviceVector result = _mulVec(pNaikLink[uiLinkIndex], pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
    const BYTE etamu = ((pEtaTable[uiSiteIndex] >> dir) & 1);
    if (_UBOOLXOR(x_p_mu_Fermion.NeedToOpposite(), etamu))
    {
        _oppo(result);
    }
    if (_UBOOLXOR(x_m_mu_Fermion.NeedToOpposite(), etamu))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_m_mu_Fermion.m_uiSiteIndex, x_m_mu_Fermion.m_byDir)], pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_m_mu_Fermion.m_uiSiteIndex, x_m_mu_Fermion.m_byDir)], pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);

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

    if (0 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (1 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (2 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (3 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _add(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(fCoeff, fCoefficient * (1 - (bDDagger << 1)));
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(cCoeff, fCoefficient * (1 - (bDDagger << 1)));
        _mul(result, cCoeff);
        break;
    default:
        _mul(result, fCoefficient * (1 - (bDDagger << 1)));
        break;
    }

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2_R(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    Real fCoeff)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _add(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
    }

    _mul(result, fCoeff);
    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2_C(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    CLGComplex cCoeff)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _add(result, _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
        else
        {
            _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
        }
    }

    _mul(result, cCoeff);

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionNaik(
    Real fCoeff,
    const deviceVector* __restrict__ n,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pFermionMove,
    deviceGauge* res,
    BYTE byFieldId)
{
    intokernalDir;

    SIndex sn1 = pFermionMove[uiLinkIndex << 1];
    deviceGauge resv = _makeContract<deviceGauge, deviceVector>(n[sn1.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(resv, fCoeff);
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;

    //This is not eta, but even-odd
    if (pEtaTable[uiSiteIndex] >> 4)
    {
        eta = eta + 1;
    }
    if (sn1.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _sub(res[uiLinkIndex], resv);
    }
    else
    {
        _add(res[uiLinkIndex], resv);
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionNaik2(
    Real fCoeff,
    const deviceVector* __restrict__ n,
    const deviceVector* __restrict__ n_p_mu,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pFermionMove,
    deviceGauge* res,
    BYTE byFieldId)
{
    intokernalDir;

    SIndex sn1 = pFermionMove[uiLinkIndex << 1];
    deviceGauge resv = _makeContract<deviceGauge, deviceVector>(n_p_mu[sn1.m_uiSiteIndex], n[uiSiteIndex]);
    _sub(resv, _makeContract<deviceGauge, deviceVector>(n[sn1.m_uiSiteIndex], n_p_mu[uiSiteIndex]));
    _mul(resv, fCoeff);
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;

    if (sn1.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _sub(res[uiLinkIndex], resv);
    }
    else
    {
        _add(res[uiLinkIndex], resv);
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorNaik(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer, 
    BYTE byFieldId, BYTE byGaugeFieldId, Real fNaik, Real fEpsilon,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    const deviceGauge* pNaikLink = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetNaikLink()->GetData());
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelDFermionKSNaikCached TMPARG(deviceVector, deviceGauge), block, threads,
        (const deviceVector*)pBuffer,
        pNaikLink,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (deviceVector*)pTargetBuffer,
        byFieldId,
        byGaugeFieldId,
        bDagger,
        fNaik,
        eOCT,
        fRealCoeff,
        cCmpCoeff
    );
    if (abs(fEpsilon) > _CLG_FLT_EPSILON)
    {
        const deviceGauge* pEffectiveGaugeLevel1 = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCached TMPARG(deviceVector, deviceGauge), block, threads,
            (const deviceVector*)pBuffer,
            pEffectiveGaugeLevel1,
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            (deviceVector*)pTargetBuffer,
            byFieldId,
            byGaugeFieldId,
            bDagger,
            fEpsilon,
            eOCT,
            fRealCoeff,
            cCmpCoeff
        );
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorNaikOnEvenOrOdd(deviceVector* pTargetBuffer,
    BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real fNaik, Real fEpsilon,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    const deviceGauge* pNaikLink = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetNaikLink()->GetData());
    //preparethreadDir;
    //_LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd, block, threads, 
    //    (deviceSU3Vector*)pTargetBuffer,
    //    pNaikLink,
    //    appGetLattice()->m_pIndexCache->m_pNaikCache[m_byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    m_byFieldId,
    //    byGaugeFieldId,
    //    bEven,
    //    bDagger,
    //    m_fNaik,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff
    //    );

    preparethreadHalf;
    if (EOCT_None == eOCT || EOCT_Real == eOCT)
    {

        Real fCoef = bDagger ? (-fNaik) : (fNaik);
        if (EOCT_Real == eOCT)
        {
            fCoef *= fRealCoeff;
        }
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_R TMPARG(deviceVector, deviceGauge), block, threads,
            (deviceVector*)pTargetBuffer,
            pNaikLink,
            appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            byFieldId,
            byGaugeFieldId,
            bEven,
            fCoef
        );
    }
    else if (EOCT_Complex == eOCT)
    {
        CLGComplex cCoeff = cuCmulf_cr(cCmpCoeff, bDagger ? -fNaik : fNaik);
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_C TMPARG(deviceVector, deviceGauge), block, threads,
            (deviceVector*)pTargetBuffer,
            pNaikLink,
            appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            byFieldId,
            byGaugeFieldId,
            bEven,
            cCoeff
        );
    }
    if (abs(fEpsilon) > _CLG_FLT_EPSILON)
    {
        const deviceGauge* pEffectiveGaugeLevel1 = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
        if (EOCT_None == eOCT || EOCT_Real == eOCT)
        {
            Real fCoef = bDagger ? (-fEpsilon) : (fEpsilon);
            if (EOCT_Real == eOCT)
            {
                fCoef *= fRealCoeff;
            }
            _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_R TMPARG(deviceVector, deviceGauge), block, threads,
                (deviceVector*)pTargetBuffer,
                pEffectiveGaugeLevel1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                byFieldId,
                byGaugeFieldId,
                bEven,
                fCoef
            );
        }
        else if (EOCT_Complex == eOCT)
        {
            CLGComplex cCoeff = cuCmulf_cr(cCmpCoeff, bDagger ? -fEpsilon : fEpsilon);
            _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_C TMPARG(deviceVector, deviceGauge), block, threads,
                (deviceVector*)pTargetBuffer,
                pEffectiveGaugeLevel1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                byFieldId,
                byGaugeFieldId,
                bEven,
                cCoeff
            );
        }
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::NaikConnection(Real fNumerator, const deviceVector* rfield, deviceGauge* naikforce, BYTE byFieldId)
{
    _RECORD(CFieldFermionKSTKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        naikforce,
        byFieldId
    );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::NaikConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId)
{
    _RECORD(CFieldFermionKSTKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik2 TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        rfield_p_mu,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        naikforce,
        byFieldId
    );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::AddConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId)
{
    _RECORD(CFieldFermionKSTKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik2 TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        rfield_p_mu,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        naikforce,
        byFieldId
    );
}


#pragma endregion

template class CFieldFermionKSTKernel<CLGComplex, CLGComplex, 1>;
template class CFieldFermionKSTKernel<deviceSU2Vector, deviceSU2, 2>;
template class CFieldFermionKSTKernel<deviceSU3Vector, deviceSU3, 3>;

template class CFieldFermionKSTKernel<deviceSU4Vector, deviceSU4, 4>;
//template class CFieldFermionKSTKernel<deviceSU5Vector, deviceSU5, 5>;
//template class CFieldFermionKSTKernel<deviceSU6Vector, deviceSU6, 6>;
//template class CFieldFermionKSTKernel<deviceSU7Vector, deviceSU7, 7>;
//template class CFieldFermionKSTKernel<deviceSU8Vector, deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================