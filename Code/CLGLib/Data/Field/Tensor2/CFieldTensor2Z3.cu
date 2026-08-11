//=============================================================================
// FILENAME : CFieldTensor2Z3.cu
//
// DESCRIPTION:
// Implementation of the discrete Z3 tensor2 (plaquette) field.
//
// The generic CFieldTensor2Kernel<deviceZN<3>> initializer would silently
// write the group identity for EFIT_RandomGaussian / EFIT_RandomGenerator,
// which is unacceptable for a discrete field; therefore a dedicated kernel is
// provided here. File loading canonicalizes every element to the nearest Z3
// root and rejects values beyond the tolerance.
//
// REVISION:
//  [08/09/26]
//=============================================================================
#include "CLGLib_Private.h"
#include "CFieldTensor2Z3.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldTensor2Z3)

#pragma region kernels

/**
* Initialize a Z3 tensor2 field.
*
* The layout is component-major, data[plaqutteIndex * siteCount + siteIndex].
* The random seed table has only Dir streams per site, so one thread owns one
* site and draws its six plaquette values from the same stream sequentially
* (identical to the generic tensor2 initializer).
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialTensor2Z3(
    deviceZN<3>* pDevicePtr,
    UINT uiSiteCount,
    UINT uiPlaqutteCount,
    EFieldInitialType eInitialType)
{
    intokernal;

    const UINT uiSeedIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    for (UINT i = 0; i < uiPlaqutteCount; ++i)
    {
        switch (eInitialType)
        {
        case EFIT_Zero:
        case EFIT_Identity:
            {
                // both map to the group identity b=0, never to complex 0
                pDevicePtr[i * uiSiteCount + uiSiteIndex] = deviceZN<3>::makeAsK(0);
            }
            break;
        case EFIT_Random:
            {
                // draw one of the three roots; flatness (dB=0) is NOT enforced
                // here, the action checks it when AllowMonopole=0
                const UINT k = static_cast<UINT>(_deviceRandomF(uiSeedIndex) * F(3.0)) % 3;
                pDevicePtr[i * uiSiteCount + uiSiteIndex] = deviceZN<3>::makeAsK(k);
            }
            break;
        case EFIT_RandomGaussian:
        case EFIT_RandomGenerator:
            {
                // should never be reached: InitialField rejects these types on
                // the host side before launching. Make any leak visible.
                printf("Tensor2Z3: forbidden initial type %d reached device kernel!\n", static_cast<INT>(eInitialType));
                pDevicePtr[i * uiSiteCount + uiSiteIndex] = deviceZN<3>::makeAsK(0);
            }
            break;
        default:
            {
                printf("Tensor2Z3 field cannot be initialized with this type! %d\n", static_cast<INT>(eInitialType));
            }
            break;
        }
    }
}


/**
* Static device helpers for the monopole diagnostic (local to this TU).
*/
static __device__ __inline__ UINT _Z3FieldToIndex(const deviceZN<3>& z)
{
    Real fBest = -F(2.0);
    UINT kBest = 0;
    for (UINT k = 0; k < 3; ++k)
    {
        const Real fTheta = F(2.0) * PI * static_cast<Real>(k) / F(3.0);
        const CLGComplex zetaConj = _make_cuComplex(_cos(fTheta), -_sin(fTheta));
        const CLGComplex prod = _cuCmulf(z.m_me, zetaConj);
        if (prod.x > fBest)
        {
            fBest = prod.x;
            kBest = k;
        }
    }
    return kBest;
}

static __device__ __inline__ BYTE _Z3FieldPlaqIndex(BYTE mu, BYTE nu)
{
    return static_cast<BYTE>(mu * (2 * _DC_Dir - mu - 1) / 2 + (nu - mu - 1));
}

static __device__ __inline__ UINT _Z3FieldSiteMove(UINT uiSiteIndex, BYTE byDir, INT iDelta)
{
    const SSmallInt4 s = __deviceSiteIndexToInt4(uiSiteIndex);
    const INT iLen = static_cast<INT>(_constIntegers[ECI_Lx + byDir]);
    SSmallInt4 t = s;
    INT iNew = 0;
    if (0 == byDir) iNew = static_cast<INT>(s.x) + iDelta;
    else if (1 == byDir) iNew = static_cast<INT>(s.y) + iDelta;
    else if (2 == byDir) iNew = static_cast<INT>(s.z) + iDelta;
    else iNew = static_cast<INT>(s.w) + iDelta;
    iNew = (iNew + iLen) % iLen;
    if (iNew < 0) iNew += iLen;
    if (0 == byDir) t.x = static_cast<SCHAR>(iNew);
    else if (1 == byDir) t.y = static_cast<SCHAR>(iNew);
    else if (2 == byDir) t.z = static_cast<SCHAR>(iNew);
    else t.w = static_cast<SCHAR>(iNew);
    return _deviceGetSiteIndex(t);
}

/**
* Monopole diagnostic: for every cube (x, mu<nu<rho) compute the Z3 charge
*   m = b_nu_rho(x+mu) - b_nu_rho(x) - b_mu_rho(x+nu) + b_mu_rho(x)
*     + b_mu_nu(x+rho) - b_mu_nu(x)   (mod 3)
* and count the non-zero cubes per site (4 orientations: xyz, xyt, xzt, yzt).
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCountMonopoleZ3Field(
    const deviceZN<3>* __restrict__ pBoundary,
    UINT uiSiteCount,
    UINT* pCount)
{
    intokernal;

    UINT uiLocal = 0;
    const UINT aiOrient[4][3] = {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 },
    };

    for (UINT o = 0; o < 4; ++o)
    {
        const BYTE mu = static_cast<BYTE>(aiOrient[o][0]);
        const BYTE nu = static_cast<BYTE>(aiOrient[o][1]);
        const BYTE rho = static_cast<BYTE>(aiOrient[o][2]);

        const UINT bNr_xm = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(nu, rho) * uiSiteCount + _Z3FieldSiteMove(uiSiteIndex, mu, 1)]);
        const UINT bNr_x  = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(nu, rho) * uiSiteCount + uiSiteIndex]);
        const UINT bMr_xn = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(mu, rho) * uiSiteCount + _Z3FieldSiteMove(uiSiteIndex, nu, 1)]);
        const UINT bMr_x  = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(mu, rho) * uiSiteCount + uiSiteIndex]);
        const UINT bMn_xr = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(mu, nu) * uiSiteCount + _Z3FieldSiteMove(uiSiteIndex, rho, 1)]);
        const UINT bMn_x  = _Z3FieldToIndex(pBoundary[_Z3FieldPlaqIndex(mu, nu) * uiSiteCount + uiSiteIndex]);

        INT m = static_cast<INT>(bNr_xm) - static_cast<INT>(bNr_x)
              - static_cast<INT>(bMr_xn) + static_cast<INT>(bMr_x)
              + static_cast<INT>(bMn_xr) - static_cast<INT>(bMn_x);
        m %= 3;
        if (m != 0)
        {
            ++uiLocal;
        }
    }

    pCount[uiSiteIndex] = uiLocal;
}

#pragma endregion


void CFieldTensor2Z3::InitialField(EFieldInitialType eInitialType)
{
    if (EFIT_RandomGaussian == eInitialType || EFIT_RandomGenerator == eInitialType)
    {
        appCrucial(_T("CFieldTensor2Z3: EFIT_RandomGaussian / EFIT_RandomGenerator are not supported on a discrete Z3 field.\n"));
        _FAIL_EXIT;
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelInitialTensor2Z3, block, threads,
        m_pDeviceData, m_uiSiteCount, PlaqutteCountPerSite(), eInitialType);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CFieldTensor2Z3::InitialWithByte(BYTE* byData)
{
    const UINT uiElementCount = GetElementCount();
    const UINT uiFloatCount = FloatN() * m_uiSiteCount;
    if (NULL == byData || uiFloatCount != uiElementCount * 2)
    {
        appCrucial(_T("CFieldTensor2Z3::InitialWithByte: invalid data (expected %d reals, %d provided).\n"), uiFloatCount, byData == NULL ? 0 : uiFloatCount);
        _FAIL_EXIT;
        return;
    }

    const Real* pRealData = (const Real*)byData;

    // tolerance: 64 * epsilon of the source scalar type
    const DOUBLE dTolerance = 64.0 * static_cast<DOUBLE>(_CLG_FLT_EPSILON);

    // deviceZN<3> is a POD holding one CLGComplex; build the host side as a
    // plain complex array (deviceZN<3>'s member functions are __device__ only)
    CLGComplex* pHostData = new CLGComplex[uiElementCount];

    for (UINT e = 0; e < uiElementCount; ++e)
    {
        const DOUBLE dRe = static_cast<DOUBLE>(pRealData[2 * e + 0]);
        const DOUBLE dIm = static_cast<DOUBLE>(pRealData[2 * e + 1]);

        UINT kBest = 0;
        DOUBLE dBest = -1.0;
        for (UINT k = 0; k < 3; ++k)
        {
            const DOUBLE dTheta = 2.0 * PI * static_cast<DOUBLE>(k) / 3.0;
            const DOUBLE dDr = dRe - cos(dTheta);
            const DOUBLE dDi = dIm - sin(dTheta);
            const DOUBLE dDist2 = dDr * dDr + dDi * dDi;
            if (dBest < 0.0 || dDist2 < dBest)
            {
                dBest = dDist2;
                kBest = k;
            }
        }

        if (sqrt(dBest) > dTolerance)
        {
            appCrucial(_T("CFieldTensor2Z3::InitialWithByte: element %d = (%f, %f) is not a Z3 root (distance %e > tol %e).\n"),
                e, dRe, dIm, sqrt(dBest), dTolerance);
            delete[] pHostData;
            _FAIL_EXIT;
            return;
        }

        const DOUBLE dTheta = 2.0 * PI * static_cast<DOUBLE>(kBest) / 3.0;
        pHostData[e] = make_cuComplex(static_cast<Real>(cos(dTheta)), static_cast<Real>(sin(dTheta)));
    }

    checkCudaErrors(cudaMemcpy(m_pDeviceData, pHostData, sizeof(CLGComplex) * uiElementCount, cudaMemcpyHostToDevice));
    delete[] pHostData;
}


UINT CFieldTensor2Z3::CountMonopoles() const
{
    preparethread;
    UINT* pDeviceCount = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceCount, sizeof(UINT) * _HC_Volume));

    _LAUNCH_KERNEL(_kernelCountMonopoleZ3Field, block, threads,
        m_pDeviceData, m_uiSiteCount, pDeviceCount);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    UINT* pHostCount = new UINT[_HC_Volume];
    checkCudaErrors(cudaMemcpy(pHostCount, pDeviceCount, sizeof(UINT) * _HC_Volume, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(pDeviceCount));

    UINT uiTotal = 0;
    for (UINT i = 0; i < _HC_Volume; ++i)
    {
        uiTotal += pHostCount[i];
    }
    delete[] pHostCount;
    return uiTotal;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
