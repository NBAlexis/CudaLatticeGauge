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

__device__ uint2 CBoundaryConditionTorusSquare::_devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const
{
    UINT xyzt[4];
    xyzt[0] = site.x < 0 ? (m_pLattice->m_uiLatticeLength[0] - 1) : (site.x >= m_pLattice->m_uiLatticeLength[0] ? 0 : site.x);
    xyzt[1] = site.y < 0 ? (m_pLattice->m_uiLatticeLength[1] - 1) : (site.y >= m_pLattice->m_uiLatticeLength[1] ? 0 : site.y);
    xyzt[2] = site.z < 0 ? (m_pLattice->m_uiLatticeLength[2] - 1) : (site.z >= m_pLattice->m_uiLatticeLength[2] ? 0 : site.z);
    xyzt[3] = site.w < 0 ? (m_pLattice->m_uiLatticeLength[3] - 1) : (site.w >= m_pLattice->m_uiLatticeLength[3] ? 0 : site.w);
    uint2 ret;
    ret.x = xyzt[0] * m_pLattice->m_uiLatticeMultipy[0] + xyzt[1] * m_pLattice->m_uiLatticeMultipy[1] + xyzt[2] * m_pLattice->m_uiLatticeMultipy[2] + xyzt[3];
    ret.y = 0;
    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
