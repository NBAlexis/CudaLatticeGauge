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

__device__ 
SIndex CBoundaryConditionTorusSquare::_devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const
{
    UINT xyzt[4];
    xyzt[0] = site.x < 0 ? (_DC_Lx - 1) : (site.x >= _DC_Lx ? 0 : site.x);
    xyzt[1] = site.y < 0 ? (_DC_Ly - 1) : (site.y >= _DC_Ly ? 0 : site.y);
    xyzt[2] = site.z < 0 ? (_DC_Lz - 1) : (site.z >= _DC_Lz ? 0 : site.z);
    xyzt[3] = site.w < 0 ? (_DC_Lt - 1) : (site.w >= _DC_Lt ? 0 : site.w);

    return SIndex(xyzt[0] * _DC_MultX
                + xyzt[1] * _DC_MultY
                + xyzt[2] * _DC_MultZ
                + xyzt[3]);
}

/**
* Torus do not have boundary fields
* Just call _devcieGetMappedIndex
*/
__device__ 
SIndex CBoundaryConditionTorusSquare::_devcieGetFermionMappedIndex(BYTE byFieldId, const int4 &site, const int4 &fromsite) const
{
    return _devcieGetMappedIndex(site, fromsite);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
