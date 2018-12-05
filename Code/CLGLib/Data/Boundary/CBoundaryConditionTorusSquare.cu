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

__device__ uint2 CBoundaryConditionTorusSquare::GetMappedIndex(const int4 &site, const int4 &fromsite, const UINT* length, const UINT* mult)
{
    UINT xyzt[4];
    xyzt[0] = site.x < 0 ? (length[0] - 1) : (site.x >= length[0] ? 0 : site.x);
    xyzt[1] = site.y < 0 ? (length[1] - 1) : (site.y >= length[1] ? 0 : site.y);
    xyzt[2] = site.z < 0 ? (length[2] - 1) : (site.z >= length[2] ? 0 : site.z);
    xyzt[3] = site.w < 0 ? (length[3] - 1) : (site.w >= length[3] ? 0 : site.w);
    uint2 ret;
    ret.x = xyzt[0] * mult[0] + xyzt[1] * mult[1] + xyzt[2] * mult[2] + xyzt[3] * mult[3];
    ret.y = 0;
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
