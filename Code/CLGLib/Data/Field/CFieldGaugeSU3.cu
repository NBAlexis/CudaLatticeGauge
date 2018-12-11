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

#pragma region BLAS

#pragma region BLAS device



#pragma endregion BLAS device

#pragma region BLAS kernel

/**
* Initial SU3 Field with a value
*/
__global__
void _kernelInitialSU3Feield(CDeviceLattice* pLattice, deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
{
    const deviceSU3 id = deviceSU3::makeSU3Id();
    const deviceSU3 zero = deviceSU3::makeSU3Zero();

    gaugeSU3KernelFuncionStart

        switch (eInitialType)
        {
            case EFIT_Zero:
                {
                    pDevicePtr[_deviceGetLinkIndex(pLattice, coord, idir)] = zero;
                }
                break;
            case EFIT_Identity:
                {
                    pDevicePtr[_deviceGetLinkIndex(pLattice, coord, idir)] = id;
                }
                break;
            case EFIT_Random:
                {
                    pDevicePtr[_deviceGetLinkIndex(pLattice, coord, idir)] = deviceSU3::makeSU3Random(pLattice->GetDeviceRandom(), _deviceGetFatIndex(pLattice, coord, idir + 1));
                }
                break;
            case EFIT_RandomGenerator:
                {
                    pDevicePtr[_deviceGetLinkIndex(pLattice, coord, idir)] = deviceSU3::makeSU3RandomGenerator(pLattice->GetDeviceRandom(), _deviceGetFatIndex(pLattice, coord, idir + 1));
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
void _kernelAxpy(CDeviceLattice* pLattice, deviceSU3 *pDevicePtr, const deviceSU3* x, const cuComplex& a)
{
    gaugeSU3KernelFuncionStart

        UINT uiLinkIndex = _deviceGetLinkIndex(pLattice, coord, idir);
        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].Scalec(a));

    gaugeSU3KernelFuncionEnd
}

__global__
void _kernelAxpy(CDeviceLattice* pLattice, deviceSU3 *pDevicePtr, const deviceSU3* x)
{
    gaugeSU3KernelFuncionStart

        UINT uiLinkIndex = _deviceGetLinkIndex(pLattice, coord, idir);
    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

#pragma endregion BLAS kernel

#pragma region BLAS member function

void CFieldGaugeSU3::Zero()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pLattice, m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeSU3::Indentity()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pLattice, m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeSU3::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pLattice, m_pDeviceData, EFIT_RandomGenerator);
}

void CFieldGaugeSU3::Axpy(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);

    preparethread;
    _kernelAxpy << <block, threads >> > (m_pLattice, m_pDeviceData, pSU3x->m_pDeviceData);
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
    preparethread;
    _kernelAxpy << <block, threads >> > (m_pLattice, m_pDeviceData, pSU3x->m_pDeviceData, a);
}

#pragma endregion BLAS member function

#pragma endregion BLAS



/**
* calculate Staple At Site
*/
__global__ 
void _kernelStapleAtSiteSU3(
    const CDeviceLattice* pLattice, 
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
            UINT linkIndex = _deviceGetLinkIndex(pLattice, coord, idir);

            UINT uiPlaqutteCount = 0;
            UINT uiPlaqutteLength = 0;

            //int2.x is linkIndex
            //int2.y is fieldIndex (may on bounday)
            //sign of int2.y is whether inverse
            int2* plaquttes = pLattice->m_pIndex->_deviceGetPlaquttesAtLink(uiPlaqutteCount, uiPlaqutteLength, linkIndex);

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
    const CDeviceLattice* pLattice,
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
            UINT linkIndex = _deviceGetLinkIndex(pLattice, coord, idir);

            deviceSU3 expP = pMyDeviceData[linkIndex].Exp(a, uiPrecision);
            expP.Mul(pU[linkIndex]);
            pU[linkIndex] = expP;
        }
    }
}

/**
*
*/
CFieldGaugeSU3::CFieldGaugeSU3(CLatticeData* pLattice, EFieldInitialType eInitialType)
    : CFieldGauge(pLattice)
    , m_pDeviceData(NULL)
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * pLattice->m_uiVolumn * pLattice->m_uiDir));

    preparethreadsimple;

    _kernelInitialSU3Feield << <block, threads >> > (m_pLattice, m_pDeviceData, eInitialType);
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

    preparethread;
    _kernelStapleAtSiteSU3 << <block, threads >> > (m_pLattice, m_pDeviceData, NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData, pForceSU3->m_pDeviceData, minusBetaOverN);
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
    _kernelExpMultSU3 << <block, threads >> > (m_pLattice, m_pDeviceData, a, uiPrecision, pUField->m_pDeviceData);
}


void CFieldGaugeSU3::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU3 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }
    CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);

    checkCudaErrors(cudaMemcpy(m_pDeviceData, pTargetField->m_pDeviceData, sizeof(deviceSU3) * m_pOwner->m_uiVolumn * m_pOwner->m_uiDir, cudaMemcpyDeviceToDevice));
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================