//=============================================================================
// FILENAME : CFieldBoundaryWilsonSquareSU3.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialWilsonFeield_Identity(deviceWilsonVectorSU3 *pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
}


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldBoundaryWilsonSquareSU3)

void CFieldBoundaryWilsonSquareSU3::InitialField(CParameters& param)
{
    CFieldBoundary::InitialField(param);

    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, 1, 1);
    _kernelInitialWilsonFeield_Identity << <block, thread >> > (m_pDeviceData);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================