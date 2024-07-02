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

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
template<INT N, INT NoE>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSUNFeield_Identity(deviceSUN<N, NoE> *pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = deviceSUN<N, NoE>::makeSUNId();
}


#pragma endregion


template<INT N, INT NoE>
CFieldBoundaryGaugeSUN<N, NoE>::CFieldBoundaryGaugeSUN() : CFieldBoundary()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSUN<N, NoE>) * 8 * _HC_Dir));
}

template<INT N, INT NoE>
CFieldBoundaryGaugeSUN<N, NoE>::~CFieldBoundaryGaugeSUN()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
}

template<INT N, INT NoE>
void CFieldBoundaryGaugeSUN<N, NoE>::InitialField(CParameters& param)
{
    CFieldBoundary::InitialField(param);

    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialSUNFeield_Identity << <block, thread >> > (m_pDeviceData);
}

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU4)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU5)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU6)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU7)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU8)


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================