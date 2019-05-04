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
    UINT uiSiteIndex = threadIdx.x;
    pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
}


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldBoundaryWilsonSquareSU3)

void CFieldBoundaryWilsonSquareSU3::InitialField(CParameters& )
{
    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, 1, 1);
    _kernelInitialWilsonFeield_Identity << <block, thread >> > (m_pDeviceData);
}

CCString CFieldBoundaryWilsonSquareSU3::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldBoundaryWilsonSquareSU3\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================