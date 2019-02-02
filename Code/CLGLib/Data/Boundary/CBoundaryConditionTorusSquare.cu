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
    UINT x,y,z,t;
    x = site.x < 0 ? (_DC_Lx - 1) : (site.x >= _DC_Lx ? 0 : site.x);
    y = site.y < 0 ? (_DC_Ly - 1) : (site.y >= _DC_Ly ? 0 : site.y);
    z = site.z < 0 ? (_DC_Lz - 1) : (site.z >= _DC_Lz ? 0 : site.z);
    t = site.w < 0 ? (_DC_Lt - 1) : (site.w >= _DC_Lt ? 0 : site.w);

    return SIndex(x * _DC_MultX
                + y * _DC_MultY
                + z * _DC_MultZ
                + t);
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
