//=============================================================================
// FILENAME : CFieldGaugeSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of gauge SU3
//
// The SU3 Matrix is
// 0 1 2
// 3 4 5
// 6 7 8
//
// Number of threads: < 1024
// Number of blocks: V / 1024
//
// threadIdx.xyz = xyz, and we loop for t and dir
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__
void _kernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSU3 id = deviceSU3::makeSU3Id();
    deviceSU3 zero = deviceSU3::makeSU3Zero();

    gaugeSU3KernelFuncionStart

        switch (eInitialType)
        {
        case EFIT_Zero:
        {
            pDevicePtr[_deviceGetLinkIndex(coord, idir)] = zero;
        }
        break;
        case EFIT_Identity:
        {
            pDevicePtr[_deviceGetLinkIndex(coord, idir)] = id;
        }
        break;
        case EFIT_Random:
        {
            pDevicePtr[_deviceGetLinkIndex(coord, idir)] = deviceSU3::makeSU3Random(_deviceGetFatIndex(coord, idir + 1));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[_deviceGetLinkIndex(coord, idir)] = deviceSU3::makeSU3RandomGenerator(_deviceGetFatIndex(coord, idir + 1));
        }
        break;
        default:
        {
            printf("SU3 Field cannot be initialized with this type!");
        }
        break;
        }

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, _Complex a)
{
    gaugeSU3KernelFuncionStart

    UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].Mulc(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpySU3Real(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

        UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].Mulc(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpyPlusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpyMinusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelScalarMultiplySU3Complex(deviceSU3 *pDevicePtr, _Complex a)
{
    gaugeSU3KernelFuncionStart

    UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Mul(a);

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelScalarMultiplySU3Real(deviceSU3 *pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
    pDevicePtr[uiLinkIndex].Mul(a);

    gaugeSU3KernelFuncionEnd
}

/**
* debug kernel
*/
__global__ void _kernelPrintSU3(const deviceSU3 * __restrict__ pDeviceData)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            printf("link at %d: %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i\n",
                linkIndex,
                pDeviceData[linkIndex].m_me[0].x, pDeviceData[linkIndex].m_me[0].y,
                pDeviceData[linkIndex].m_me[1].x, pDeviceData[linkIndex].m_me[1].y,
                pDeviceData[linkIndex].m_me[2].x, pDeviceData[linkIndex].m_me[2].y,
                pDeviceData[linkIndex].m_me[3].x, pDeviceData[linkIndex].m_me[3].y,
                pDeviceData[linkIndex].m_me[4].x, pDeviceData[linkIndex].m_me[4].y,
                pDeviceData[linkIndex].m_me[5].x, pDeviceData[linkIndex].m_me[5].y,
                pDeviceData[linkIndex].m_me[6].x, pDeviceData[linkIndex].m_me[6].y,
                pDeviceData[linkIndex].m_me[7].x, pDeviceData[linkIndex].m_me[7].y,
                pDeviceData[linkIndex].m_me[8].x, pDeviceData[linkIndex].m_me[8].y
            );
        }
    }
}

/**
* calculate Staple and Force At Site
*/
__global__
void _kernelStapleAtSiteSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
    _Complex minusBetaOverN)
{
    intokernal;

    minusBetaOverN.x = minusBetaOverN.x * (Real)0.5;
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            UINT uiPlaqutteCount = 0;
            UINT uiPlaqutteLength = 0;

            //int2.x is linkIndex
            //int2.y is fieldIndex (may on bounday)
            //sign of int2.y is whether inverse
            __idx->_deviceGetPlaquttesAtLink(plaquttes, uiPlaqutteCount, uiPlaqutteLength, linkIndex);

            deviceSU3 res = deviceSU3::makeSU3Zero();

            for (int i = 0; i < uiPlaqutteCount; ++i)
            {
                SIndex first = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1)];
                deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
                if (first.NeedToDagger())
                {
                    toAdd.Dagger();
                }

                for (int j = 1; j < uiPlaqutteLength - 1; ++j)
                {
                    SIndex nextlink = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1) + j];
                    deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                    if (nextlink.NeedToDagger())
                    {
                        toMul.Dagger();
                    }
                    toAdd.Mul(toMul);
                }
                res.Add(toAdd);
            }
            if (NULL != pStapleData)
            {
                pStapleData[linkIndex] = res;
            }

            res.Dagger();
            //staple calculated
            deviceSU3 force(pDeviceData[linkIndex]);
            force.Mul(res);
            force.TrTa();
            force.Mul(minusBetaOverN);

            //force is additive
            pForceData[linkIndex].Add(force);
        }
    }
}

/**
* calculate Staple and eneregy At Site
*/
__global__
void _kernelPlaqutteEnergySU3(
    const deviceSU3 * __restrict__ pDeviceData,
    Real minusBetaOverN,
    Real* results)
{
    intokernal;

    Real resThisThread = 0;
    SIndex plaquttes[kMaxPlaqutteCache];
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            UINT uiPlaqutteCount = 0;
            UINT uiPlaqutteLength = 0;

            //int2.x is linkIndex
            //int2.y is fieldIndex (may on bounday)
            //sign of int2.y is whether inverse
            __idx->_deviceGetPlaquttesAtLink(plaquttes, uiPlaqutteCount, uiPlaqutteLength, linkIndex);

            deviceSU3 res = deviceSU3::makeSU3Zero();

            for (int i = 0; i < uiPlaqutteCount; ++i)
            {
                SIndex first = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1)];

                deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
                if (first.NeedToDagger())
                {
                    toAdd.Dagger();
                }

                for (int j = 1; j < uiPlaqutteLength - 1; ++j)
                {
                    SIndex nextlink = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1) + j];
                    deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                    if (nextlink.NeedToDagger())
                    {
                        toMul.Dagger();
                    }
                    toAdd.Mul(toMul);
                }
                res.Add(toAdd);
            }

            res.Dagger();
            deviceSU3 gaugeOnLink(pDeviceData[linkIndex]);
            gaugeOnLink.Mul(res);
            resThisThread += (3 - gaugeOnLink.ReTr()); //Re[Tr(1-U)] = 3 - Re[Tr(U)]
        }
    }

    results[threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z] = resThisThread * minusBetaOverN;
}

/**
*
*/
__global__
void _kernelExpMultSU3(
    const deviceSU3 * __restrict__ pMyDeviceData,
    _Complex a,
    deviceSU3 *pU)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);

            deviceSU3 expP = pMyDeviceData[linkIndex].Exp(a, _DC_ExpPrecision);
            expP.Mul(pU[linkIndex]);
            expP.Norm();
            pU[linkIndex] = expP;
        }
    }
}


/**
* Trace (P^2)
*/
__global__ 
void _kernelCalculateKinematicEnergySU3(const deviceSU3 * __restrict__ pMyDeviceData, Real* results)
{
    intokernal;

    Real resThisThread = 0;
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            resThisThread += pMyDeviceData[linkIndex].Daggerc().Mulc(pMyDeviceData[linkIndex]).ReTr();
        }
    }

    results[threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z] = resThisThread;
}

#pragma endregion

#pragma region CFunctions

extern "C" {

    void _callKernelAxpySU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, UBOOL bAdd)
    {
        preparethread;
        if (bAdd)
        {
            _kernelAxpyPlusSU3 << <block, threads >> > (pDevicePtr, x);
        }
        else 
        {
            _kernelAxpyMinusSU3 << <block, threads >> > (pDevicePtr, x);
        }
    }

    void _callKernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, const _Complex& a)
    {
        preparethread;
        _kernelAxpySU3A << <block, threads >> > (pDevicePtr, x, a);
    }

    void _callKernelAxpySU3Real(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, Real a)
    {
        preparethread;
        _kernelAxpySU3Real << <block, threads >> > (pDevicePtr, x, a);
    }

    void _callKernelScalarMultiplyComplex(deviceSU3 *pDevicePtr, const _Complex& a)
    {
        preparethread;
        _kernelScalarMultiplySU3Complex << <block, threads >> > (pDevicePtr, a);
    }

    void _callKernelScalarMultiplyReal(deviceSU3 *pDevicePtr, Real a)
    {
        preparethread;
        _kernelScalarMultiplySU3Real << <block, threads >> > (pDevicePtr, a);
    }

    void _callKernelPrint(const deviceSU3 * __restrict__ pDeviceData)
    {
        preparethread;
        _kernelPrintSU3 << < block, threads >> > (pDeviceData);
    }

    void _callKernelExpMultSU3(
        const deviceSU3 * __restrict__ pMyDeviceData,
        const _Complex& a,
        deviceSU3 *pU)
    {
        preparethread;
        _kernelExpMultSU3 << < block, threads >> > (pMyDeviceData, a, pU);
    }

    void _callKernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
    {
        preparethread;
        _kernelInitialSU3Feield << <block, threads >> > (pDevicePtr, eInitialType);
    }

    void _callKernelStapleAtSiteSU3(
        const deviceSU3 * __restrict__ pDeviceData,
        deviceSU3 *pStapleData, //can be NULL
        deviceSU3 *pForceData,
        const _Complex& minusBetaOverN)
    {
        preparethread;
        _kernelStapleAtSiteSU3 << <block, threads >> > (pDeviceData, pStapleData, pForceData, minusBetaOverN);
    }

    Real _callKernelPlaqutteEneregySU3(
        const deviceSU3 * __restrict__ pDeviceData,
        const _Complex& minusBetaOverN,
        UINT uiDataSize,
        Real* output
    )
    {
        preparethread;
        _kernelPlaqutteEnergySU3 << <block, threads >> > (pDeviceData, -minusBetaOverN.x, output);

        //Now the results are in output.
        thrust::device_ptr<Real> dp(output);
        thrust::device_vector<Real> d_x(dp, dp + uiDataSize);
        return thrust::reduce(d_x.begin(), d_x.end(), (Real)0, thrust::plus<Real>());
    }

    Real _callKernelKinematicEnergySU3(const deviceSU3 * __restrict__ pDeviceData, Real* output, UINT uiDataSize)
    {
        preparethread;
        _kernelCalculateKinematicEnergySU3 << <block, threads >> > (pDeviceData, output);
        thrust::device_ptr<Real> dp(output);
        thrust::device_vector<Real> d_x(dp, dp + uiDataSize);
        return thrust::reduce(d_x.begin(), d_x.end(), (Real)0, thrust::plus<Real>());
    }
}

#pragma endregion

void CFieldGaugeSU3::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3(m_pDeviceData, pSU3x->m_pDeviceData, TRUE);

}

void CFieldGaugeSU3::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3(m_pDeviceData, pSU3x->m_pDeviceData, FALSE);

}

void CFieldGaugeSU3::ScalarMultply(const _Complex& a)
{
    _callKernelScalarMultiplyComplex(m_pDeviceData, a);
}

void CFieldGaugeSU3::ScalarMultply(Real a)
{
    _callKernelScalarMultiplyReal(m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3Real(m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(const _Complex& a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3A(m_pDeviceData, pSU3x->m_pDeviceData, a);
}


void CFieldGaugeSU3::Zero()
{
    _callKernelInitialSU3Feield(m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeSU3::Indentity()
{
    _callKernelInitialSU3Feield(m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeSU3::MakeRandomGenerator()
{
    _callKernelInitialSU3Feield(m_pDeviceData, EFIT_RandomGenerator);
}

/**
*
*/
void CFieldGaugeSU3::InitialField(EFieldInitialType eInitialType)
{
    _callKernelInitialSU3Feield(m_pDeviceData, eInitialType);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeSU3::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, const _Complex& minusBetaOverN) const
{
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    if (NULL != pStable && EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return;
    }

    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    CFieldGaugeSU3* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSU3*>(pStable);

    _callKernelStapleAtSiteSU3(m_pDeviceData, NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData, pForceSU3->m_pDeviceData, minusBetaOverN);
}

Real CFieldGaugeSU3::CalculatePlaqutteEnergy(const _Complex& minusBetaOverN) const
{
    return _callKernelPlaqutteEneregySU3(m_pDeviceData, minusBetaOverN, _HC_ThreadCount, _D_RealThreadBuffer) / m_uiLinkeCount;
}

Real CFieldGaugeSU3::CalculateKinematicEnergy() const
{
    return _callKernelKinematicEnergySU3(m_pDeviceData, _D_RealThreadBuffer, _HC_ThreadCount);
}

CFieldGaugeSU3::CFieldGaugeSU3()
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount));
}

CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

void CFieldGaugeSU3::ExpMult(const _Complex& a, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    _callKernelExpMultSU3(m_pDeviceData, a, pUField->m_pDeviceData);
}


void CFieldGaugeSU3::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU3 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }
    CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU3::DebugPrintMe() const
{
    _callKernelPrint(m_pDeviceData);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================