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
#include "CFieldBoundaryOne.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
template<typename deviceData>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialBoundaryFeield_Identity(deviceData* pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = _makeId<deviceData>();
}


#pragma endregion


template<typename deviceData>
CFieldBoundaryOne<deviceData>::CFieldBoundaryOne() : CFieldBoundary()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceData) * 8 * _HC_Dir));
}

template<typename deviceData>
CFieldBoundaryOne<deviceData>::~CFieldBoundaryOne()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

template<typename deviceData>
void CFieldBoundaryOne<deviceData>::InitialField(CParameters& param)
{
    CFieldBoundary::InitialField(param);

    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialBoundaryFeield_Identity << <block, thread >> > (m_pDeviceData);
}

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeU1)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU2)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU3)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU4)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU5)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU6)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU7)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU8)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================