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

__device__ __inline__ UINT _deviceGetIndexStart(const UINT* coord, const UINT* mult)
{
    return coord[0] * mult[0] + coord[1] * mult[1] + coord[2] * mult[2] + coord[3] * mult[3] + coord[4];
}

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
__global__ void _kernelInitialSU3Feield(cuComplex *pDevicePtr, UINT tLength, UINT dir, const UINT* mult)
{
    UINT coord[5];
    coord[0] = threadIdx.x + blockIdx.x * blockDim.x;
    coord[1] = threadIdx.y + blockIdx.y * blockDim.y;
    coord[2] = threadIdx.z + blockIdx.z * blockDim.z;

    for (UINT it = 0; it < tLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < idir; ++idir)
        {
            coord[4] = idir;
            UINT indexStart = _deviceGetIndexStart(coord, mult);
            pDevicePtr[indexStart + 0] = make_cuComplex(1.0f, 0.0f);
            pDevicePtr[indexStart + 1] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 2] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 3] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 4] = make_cuComplex(1.0f, 0.0f);
            pDevicePtr[indexStart + 5] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 6] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 7] = make_cuComplex(0.0f, 0.0f);
            pDevicePtr[indexStart + 8] = make_cuComplex(1.0f, 0.0f);
        }
    }
}

/**
* y = a * x + y
*/
__global__ void _kernelAxpy(cuComplex *pDevicePtry, const cuComplex *pDevicePtrx, const cuComplex& a, UINT tLength, UINT dir, const UINT* mult)
{
    UINT coord[5];
    coord[0] = threadIdx.x + blockIdx.x * blockDim.x;
    coord[1] = threadIdx.y + blockIdx.y * blockDim.y;
    coord[2] = threadIdx.z + blockIdx.z * blockDim.z;

    for (UINT it = 0; it < tLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < idir; ++idir)
        {
            coord[4] = idir;
            UINT indexStart = _deviceGetIndexStart(coord, mult);
            _deviceAxpy(pDevicePtry + indexStart, a, pDevicePtrx + indexStart);
        }
    }
}

/**
* calculate Staple At Site
*/
__global__ void _kernelStapleAtSite(cuComplex *pDeviceData, cuComplex *pStapleData, CLatticeData * pLattice, UINT tLength, const UINT* length, const UINT* mult)
{
    UINT coord[5];
    coord[0] = threadIdx.x + blockIdx.x * blockDim.x;
    coord[1] = threadIdx.y + blockIdx.y * blockDim.y;
    coord[2] = threadIdx.z + blockIdx.z * blockDim.z;

    for (UINT it = 0; it < tLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < idir; ++idir)
        {
            coord[4] = idir;

            UINT linkIndex = _deviceGetIndexStart(coord, mult) / 9;

            UINT uiPlaqutteCount = 0;
            UINT uiPlaqutteLength = 0;
            int2* plaquttes = pLattice->m_pIndex->GetPlaquttesAtLink(uiPlaqutteCount, uiPlaqutteLength, pLattice->m_uiDim, linkIndex, length, mult);

            for (int i = 0; i < uiPlaqutteCount; ++i)
            {

            }
        }
    }
}

/**
*
*/
CFieldGaugeSU3::CFieldGaugeSU3(CLatticeData* pLattice)
    : CFieldGauge(pLattice)
    , m_pDeviceData(NULL)
    , m_pDeviceStaple(NULL)
    , m_pDevicePlaquetteEnergy(NULL)
{
    for (int i = 0; i < CCommonData::kMaxDim + 1; ++i)
    {
        //9 elements for every link
        m_uiLatticeMultipy[i] *= 9;
    }

    UINT volumn = m_uiLatticeLength[CCommonData::kMaxDim - 1]
        * m_uiLatticeDecompose[0]
        * m_uiLatticeDecompose[1]
        * m_uiLatticeDecompose[2]
        * m_uiLatticeDecompose[3]
        * m_uiLatticeDecompose[4]
        * m_uiLatticeDecompose[5];

    checkCudaErrors(cudaMalloc((void **)&m_pDeviceData, sizeof(cuComplex) * volumn * m_uiDir * 9));
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceStaple, sizeof(cuComplex) * volumn * m_uiDir * 9));
    checkCudaErrors(cudaMalloc((void **)&m_pDevicePlaquetteEnergy, sizeof(FLOAT) * volumn * m_uiDir));

    dim3 block(m_uiLatticeDecompose[0], m_uiLatticeDecompose[1], m_uiLatticeDecompose[2]);
    dim3 threads(m_uiLatticeDecompose[3], m_uiLatticeDecompose[4], m_uiLatticeDecompose[5]);
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, m_uiLatticeLength[CCommonData::kMaxDim - 1], m_uiDir, m_uiLatticeMultipy);

    checkCudaErrors(cudaMemcpy(m_pDeviceStaple, m_pDeviceData, sizeof(cuComplex) * volumn * m_uiDir * 9, cudaMemcpyDeviceToDevice));
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
void CFieldGaugeSU3::axpy(const cuComplex& a, const CField *x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);

    dim3 block(m_uiLatticeDecompose[0], m_uiLatticeDecompose[1], m_uiLatticeDecompose[2]);
    dim3 threads(m_uiLatticeDecompose[3], m_uiLatticeDecompose[4], m_uiLatticeDecompose[5]);

    _kernelAxpy << <block, threads >> >(m_pDeviceData, pSU3x->m_pDeviceData, a, m_uiLatticeLength[CCommonData::kMaxDim - 1], m_uiDir, m_uiLatticeMultipy);
}

/**
*
*/
void CFieldGaugeSU3::CalculateStaple()
{
    CIndex* pIndex = m_pOwner->m_pIndex;

}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================