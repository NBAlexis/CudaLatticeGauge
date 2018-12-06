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

/**
* y = a x + y
*/
__device__ __inline__ void _deviceAxpy(cuComplex* y, const cuComplex& a, const cuComplex * x)
{
    for (int i = 0; i < 9; ++i)
    {
        y[i] = cuCaddf(y[i], cuCmulf(a, x[i]));
    }
}

/**
* Initial SU3 Field with identity matrix
*/
__global__ void _kernelInitialSU3Feield(CDeviceLattice* pLattice, deviceSU3 *pDevicePtr, UBOOL bHot)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            pDevicePtr[_deviceGetLinkIndex(pLattice, coord, idir)] = bHot ? 
                deviceSU3::makeSU3Random(pLattice->GetDeviceRandom(), _deviceGetFatIndex(pLattice, coord, idir + 1))
              : deviceSU3::makeSU3Id();
        }
    }
}

/**
* calculate Staple At Site
*/
__global__ void _kernelStapleAtSite(CDeviceLattice* pLattice, deviceSU3 *pDeviceData, deviceSU3 *pStapleData)
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
            int2* plaquttes = pLattice->m_pIndex->_deviceGetPlaquttesAtLink(uiPlaqutteCount, uiPlaqutteLength, linkIndex);

            cuComplex* res;

            for (int i = 0; i < uiPlaqutteCount; ++i)
            {
                
            }
        }
    }
}

/**
*
*/
CFieldGaugeSU3::CFieldGaugeSU3(CLatticeData* pLattice, UBOOL bHot)
    : CFieldGauge(pLattice)
    , m_pDeviceData(NULL)
    , m_pDeviceStaple(NULL)
    , m_pDevicePlaquetteEnergy(NULL)
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * pLattice->m_uiVolumn * pLattice->m_uiDir));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceStaple, sizeof(deviceSU3) * pLattice->m_uiVolumn * pLattice->m_uiDir));
    checkCudaErrors(cudaMalloc((void **)&m_pDevicePlaquetteEnergy, sizeof(FLOAT) * pLattice->m_uiVolumn * pLattice->m_uiDir));

    preparethreadsimple;

    _kernelInitialSU3Feield << <block, threads >> > (m_pLattice, m_pDeviceData, bHot);

    checkCudaErrors(cudaMemcpy(m_pDeviceStaple, m_pDeviceData, sizeof(deviceSU3) * pLattice->m_uiVolumn * pLattice->m_uiDir, cudaMemcpyDeviceToDevice));
}

/**
*
*/
CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

/**
*
*/
//void CFieldGaugeSU3::axpy(const cuComplex& a, const CField *x)
//{
//    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
//    {
//        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
//        return;
//    }
//
//    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
//
//    dim3 block(m_uiLatticeDecompose[0], m_uiLatticeDecompose[1], m_uiLatticeDecompose[2]);
//    dim3 threads(m_uiLatticeDecompose[3], m_uiLatticeDecompose[4], m_uiLatticeDecompose[5]);
//
//    _kernelAxpy << <block, threads >> >(m_pDeviceData, pSU3x->m_pDeviceData, a, m_uiLatticeLength[CCommonData::kMaxDim - 1], m_uiDir, m_uiLatticeMultipy);
//}

/**
*
*/
void CFieldGaugeSU3::CalculateStaple()
{
//    CIndex* pIndex = m_pOwner->m_pIndex;
//
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================