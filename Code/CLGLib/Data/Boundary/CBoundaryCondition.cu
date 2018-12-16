//=============================================================================
// FILENAME : CBoundaryCondition.cpp
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__global__ void _kernelCreateBoundaryCondition(void** devicePtr, UINT* size, EBoundaryCondition eBC)
{
    if (eBC == EBC_TorusSquare)
    {
        (*devicePtr) = (void*)(new CBoundaryConditionTorusSquare());
        size[0] = (UINT)sizeof(CBoundaryConditionTorusSquare);
        return;
    }
    size[0] = 0;
}

extern "C" { 
    void _cCreateBC(void** devicePtr, UINT* size, EBoundaryCondition eBC)
    {
        _kernelCreateBoundaryCondition << <1, 1 >> > (devicePtr, size, eBC);
    }
}
__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================