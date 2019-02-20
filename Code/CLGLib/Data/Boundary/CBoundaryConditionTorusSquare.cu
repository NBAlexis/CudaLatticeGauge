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

__device__ CBoundaryConditionTorusSquare::CBoundaryConditionTorusSquare() : deviceBoundaryCondition()
{
    for (UINT i = 0; i < _kMaxFieldCount; ++i)
    {
        m_FermionBC[i].x = 1;
        m_FermionBC[i].y = 1;
        m_FermionBC[i].z = 1;
        m_FermionBC[i].w = 1;
    }
}

__device__ 
SIndex CBoundaryConditionTorusSquare::_devcieGetMappedIndex(const SSmallInt4 &site, const SIndex &fromsite) const
{
    UINT x, y, z, t;
    BYTE byOldRegion = fromsite.m_byReginId;
    
    if (site.x < 0)
    {
        x = _DC_Lx - 1;
        byOldRegion = (byOldRegion & /*1110*/0x0E) | ((~byOldRegion) & 0x01);
    }
    else if (site.x >= _DC_Lx)
    {
        x = 0;
        byOldRegion = (byOldRegion & /*1110*/0x0E) | ((~byOldRegion) & 0x01);
    }
    else 
    {
        x = site.x;
    }

    if (site.y < 0)
    {
        y = _DC_Ly - 1;
        byOldRegion = (byOldRegion & /*1101*/0x0D) | ((~byOldRegion) & 0x02);
    }
    else if (site.y >= _DC_Ly)
    {
        y = 0;
        byOldRegion = (byOldRegion & /*1101*/0x0D) | ((~byOldRegion) & 0x02);
    }
    else
    {
        y = site.y;
    }

    if (site.z < 0)
    {
        z = _DC_Lz - 1;
        byOldRegion = (byOldRegion & /*1011*/0x0B) | ((~byOldRegion) & 0x04);
    }
    else if (site.z >= _DC_Lz)
    {
        z = 0;
        byOldRegion = (byOldRegion & /*1011*/0x0B) | ((~byOldRegion) & 0x04);
    }
    else
    {
        z = site.z;
    }

    if (site.w < 0)
    {
        t = _DC_Lt - 1;
        byOldRegion = (byOldRegion & /*0111*/0x07) | ((~byOldRegion) & 0x08);
    }
    else if (site.w >= _DC_Lt)
    {
        t = 0;
        byOldRegion = (byOldRegion & /*0111*/0x07) | ((~byOldRegion) & 0x08);
    }
    else
    {
        t = site.w;
    }

    return SIndex(x * _DC_MultX
                + y * _DC_MultY
                + z * _DC_MultZ
                + t, //index
        0, //dir
        0, //tag
        0, //field
        byOldRegion); //region
}

/**
* Torus do not have boundary fields
*/
__device__ 
SIndex CBoundaryConditionTorusSquare::_devcieGetFermionMappedIndex(BYTE byFieldId, const SSmallInt4 &site, const SIndex &fromsite) const
{
    SIndex ret = _devcieGetMappedIndex(site, fromsite);

    //region is (bbbb)
    //For example, the boundary condition is anti-periodic on t-direction
    //If it is from a region (0000)
    //And the result is (1000)
    //It means it go to a minus-t region, and need a opposite
    //On the other hand, if it is from a region (1000)
    //And the result is (0000)
    //It means it go back from the minus-t region to normal one
    //So whether need to opposite, is decided whether the t-bit is 1
    SBYTE retOpposite = 1;
    if (m_FermionBC[byFieldId].x < 0 && (0 != (ret.m_byReginId & 0x01)))
    {
        retOpposite = -1;
    }
    if (m_FermionBC[byFieldId].y < 0 && (0 != (ret.m_byReginId & 0x02)))
    {
        retOpposite = retOpposite * -1;
    }
    if (m_FermionBC[byFieldId].z < 0 && (0 != (ret.m_byReginId & 0x04)))
    {
        retOpposite = retOpposite * -1;
    }
    if (m_FermionBC[byFieldId].w < 0 && (0 != (ret.m_byReginId & 0x08)))
    {
        retOpposite = retOpposite * -1;
    }

    if (retOpposite < 0)
    {
        ret.m_byTag = _kOpposite;
    }
    return ret;
}

__device__ void CBoundaryConditionTorusSquare::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
{
    assert(byFieldId < _kMaxFieldCount);

    m_FermionBC[byFieldId] = bc.m_sPeriodic;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
