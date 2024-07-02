//=============================================================================
// FILENAME : CFieldBoundaryGaugeSU2.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSU2Feield_Identity(deviceSU2 *pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = deviceSU2::makeSU2Id();
}


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU2)

CFieldBoundaryGaugeSU2::CFieldBoundaryGaugeSU2() : CFieldBoundary()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSU2) * 8 * _HC_Dir));
}

void CFieldBoundaryGaugeSU2::InitialField(CParameters& param)
{
    CFieldBoundary::InitialField(param);

    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialSU2Feield_Identity << <block, thread >> > (m_pDeviceData);
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================