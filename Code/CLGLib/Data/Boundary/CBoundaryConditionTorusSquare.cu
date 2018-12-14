//=============================================================================
// FILENAME : CBoundaryConditionTorusSquare.cpp
// 
// DESCRIPTION:
// This is the periodic boundary condition
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CBoundaryConditionTorusSquare)

__device__ 
uint2 CBoundaryConditionTorusSquare::_devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const
{
    UINT xyzt[4];
    xyzt[0] = site.x < 0 ? (_DC_Lx - 1) : (site.x >= _DC_Lx ? 0 : site.x);
    xyzt[1] = site.y < 0 ? (_DC_Ly - 1) : (site.y >= _DC_Ly ? 0 : site.y);
    xyzt[2] = site.z < 0 ? (_DC_Lz - 1) : (site.z >= _DC_Lz ? 0 : site.z);
    xyzt[3] = site.w < 0 ? (_DC_Lt - 1) : (site.w >= _DC_Lt ? 0 : site.w);
    uint2 ret;
    ret.x = xyzt[0] * _DC_MultX
          + xyzt[1] * _DC_MultY
          + xyzt[2] * _DC_MultZ
          + xyzt[3];
    ret.y = 0;
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
