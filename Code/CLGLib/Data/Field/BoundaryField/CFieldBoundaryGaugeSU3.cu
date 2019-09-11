//=============================================================================
// FILENAME : CFieldBoundaryGaugeSU3.cu
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
_kernelInitialSU3Feield_Identity(deviceSU3 *pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = deviceSU3::makeSU3Id();
}


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU3)

CFieldBoundaryGaugeSU3::CFieldBoundaryGaugeSU3()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSU3) * 8 * _HC_Dir));
}

void CFieldBoundaryGaugeSU3::InitialField(CParameters& )
{
    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialSU3Feield_Identity << <block, thread >> > (m_pDeviceData);
}

CCString CFieldBoundaryGaugeSU3::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldBoundaryGaugeSU3\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================