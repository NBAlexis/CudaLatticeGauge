//=============================================================================
// FILENAME : CIndex.cu
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/16/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT)
{
    if (eIT == EIndexType_Square)
    {
        (*devicePtr) = (void*)(new CIndexSquare(*pBC));
        size[0] =(UINT)sizeof(CIndexSquare);
        return;
    }
    size[0] = 0;
}

extern "C" 
{ 
    void _cCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT)
    {
        _kernelCreateIndex << <1, 1 >> > (devicePtr, pBC, size, eIT);
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================