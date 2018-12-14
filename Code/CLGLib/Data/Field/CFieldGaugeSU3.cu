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

#pragma region BLAS

#pragma region BLAS device



#pragma endregion BLAS device

#pragma region BLAS kernel

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
void _kernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* x, const cuComplex& a)
{
    gaugeSU3KernelFuncionStart

        UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].Scalec(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpySU3(deviceSU3 *pDevicePtr, const deviceSU3* x)
{
    gaugeSU3KernelFuncionStart

        UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir);
        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

extern "C" {
    void _callKernelAxpySU3(deviceSU3 *pDevicePtr, const deviceSU3* x)
    {
        preparethread;
        _kernelAxpySU3 << <block, threads >> > (pDevicePtr, x);
    }

    void _callKernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* x, const cuComplex& a)
    {
        preparethread;
        _kernelAxpySU3A << <block, threads >> > (pDevicePtr, x, a);
    }
}

#pragma endregion BLAS kernel

#pragma region BLAS member function

void CFieldGaugeSU3::Axpy(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3(m_pDeviceData, pSU3x->m_pDeviceData);

}

void CFieldGaugeSU3::Axpy(FLOAT a, const CField* x)
{
    Axpy(make_cuComplex(a, 0.0f), x);
}

void CFieldGaugeSU3::Axpy(const cuComplex& a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    _callKernelAxpySU3A(m_pDeviceData, pSU3x->m_pDeviceData, a);
}

#pragma endregion BLAS member function

#pragma endregion BLAS


#pragma region kernels

/**
* calculate Staple At Site
*/
__global__ 
void _kernelStapleAtSiteSU3(
    const CIndex* pIndex,
    int2* plaquttes,
    const deviceSU3 *pDeviceData, 
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData, 
    const cuComplex& minusBetaOverN)
{
    intokernal;

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
            pIndex->_deviceGetPlaquttesAtLink(plaquttes, uiPlaqutteCount, uiPlaqutteLength, linkIndex);

            deviceSU3 res = deviceSU3::makeSU3Zero();

            for (int i = 0; i < uiPlaqutteCount; ++i)
            {
                int2 first = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1)];
                UBOOL bReverse = first.y < 0;
                UINT fieldId = bReverse ? (-first.y - 1) : (first.y - 1);

                deviceSU3 toAdd(pDeviceData[first.x]);
                if (bReverse)
                {
                    toAdd.Dagger();
                }

                for (int j = 1; j < uiPlaqutteLength - 1; ++j)
                {
                    int2 nextlink = plaquttes[uiPlaqutteCount * (uiPlaqutteLength - 1) + j];
                    UBOOL bNextReverse = nextlink.y < 0;
                    UINT uiNextfieldId = bNextReverse ? (-nextlink.y - 1) : (nextlink.y - 1);

                    deviceSU3 toMul(pDeviceData[first.x]);
                    if (bNextReverse)
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
*
*/
__global__
void _kernelExpMultSU3(
    const deviceSU3 *pMyDeviceData,
    const cuComplex& a,
    UINT uiPrecision,
    deviceSU3 *pU)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);

            deviceSU3 expP = pMyDeviceData[linkIndex].Exp(a, uiPrecision);
            expP.Mul(pU[linkIndex]);
            pU[linkIndex] = expP;
        }
    }
}

#pragma endregion

#pragma region Call Kernel

extern "C" {

    void _callKernelExpMultSU3(
        const deviceSU3 *pMyDeviceData,
        const cuComplex& a,
        UINT uiPrecision,
        deviceSU3 *pU)
    {
        preparethread;
        _kernelExpMultSU3 <<< block, threads >> > (pMyDeviceData, a, uiPrecision, pU);
    }

    void _callKernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
    {
        preparethread;
        _kernelInitialSU3Feield << <block, threads >> > (pDevicePtr, eInitialType);
    }

    void _callKernelStapleAtSiteSU3(
        const CIndex* pIndex,
        int2* buffer,
        const deviceSU3 *pDeviceData,
        deviceSU3 *pStapleData, //can be NULL
        deviceSU3 *pForceData,
        const cuComplex& minusBetaOverN)
    {
        preparethread;
        _kernelStapleAtSiteSU3 << <block, threads >> > (pIndex, buffer, pDeviceData, pStapleData, pForceData, minusBetaOverN);
    }

}
#pragma endregion

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
void CFieldGaugeSU3::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, const cuComplex& minusBetaOverN) const
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

    _callKernelStapleAtSiteSU3(appGetLattice()->m_pDeviceIndex, m_pDeviceStapleIndex, m_pDeviceData, NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData, pForceSU3->m_pDeviceData, minusBetaOverN);
}

CFieldGaugeSU3::CFieldGaugeSU3()
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * _HC_Volumn * _HC_Dir));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceStapleIndex, sizeof(int2) * _constMaxPlaqIndexLength));
}

CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceStapleIndex));
    checkCudaErrors(cudaFree(m_pDeviceData));
}

void CFieldGaugeSU3::ExpMult(const cuComplex& a, UINT uiPrecision, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    _callKernelExpMultSU3(m_pDeviceData, a, uiPrecision, pUField->m_pDeviceData);
}


void CFieldGaugeSU3::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU3 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }
    CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);

    checkCudaErrors(cudaMemcpy(m_pDeviceData, pTargetField->m_pDeviceData, sizeof(deviceSU3) * _HC_Volumn * _HC_Dir, cudaMemcpyDeviceToDevice));
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================