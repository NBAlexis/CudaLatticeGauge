//=============================================================================
// FILENAME : CFieldBoundaryGaugeU1.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [10/01/2021 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialU1Feield_Identity(CLGComplex *pDevicePtr)
{
    const UINT uiSiteIndex = threadIdx.x;
    const UINT uiBoundIndex = threadIdx.y;

    pDevicePtr[uiSiteIndex * _DC_Dir + uiBoundIndex] = _onec;
}


#pragma endregion

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeU1)

CFieldBoundaryGaugeU1::CFieldBoundaryGaugeU1()
{
    //8 faces and 4 directions
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(CLGComplex) * 8 * _HC_Dir));
}

void CFieldBoundaryGaugeU1::InitialField(CParameters& )
{
    //we currently, only support identity
    dim3 block(1, 1, 1);
    dim3 thread(8, _HC_Dir, 1);
    _kernelInitialU1Feield_Identity << <block, thread >> > (m_pDeviceData);
}

CCString CFieldBoundaryGaugeU1::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldBoundaryGaugeU1\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================