//=============================================================================
// FILENAME : CFieldBoundaryGaugeSUN.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldBoundaryZero.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
template<typename deviceData>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialBoundaryFeield_Zero(deviceData*pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = _makeZero<deviceData>();
}


#pragma endregion


template<typename deviceData>
CFieldBoundaryZero<deviceData>::CFieldBoundaryZero() : CFieldBoundary()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceData) * 8 * _HC_Dir));
}

template<typename deviceData>
CFieldBoundaryZero<deviceData>::~CFieldBoundaryZero()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

template<typename deviceData>
void CFieldBoundaryZero<deviceData>::InitialField(CParameters& param)
{
    CFieldBoundary::InitialField(param);

    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialBoundaryFeield_Zero << <block, thread >> > (m_pDeviceData);
}

__CLGIMPLEMENT_CLASS(CFieldBoundarySU3Vector)
__CLGIMPLEMENT_CLASS(CFieldBoundaryWilsonSquareSU3)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================