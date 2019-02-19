//=============================================================================
// FILENAME : CIndexSquare.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

/**
* For square lattice, we assume dimenssion = direction
* bdir, mu: bond direction, or link direction
* pdir, nu: plaqutte direction, the second bond direction
*
* Forward
*      c
*    +-<-+
*  d |   | b ^
*    *->-+
*      a
*  a = [site][b_dir]
*  b = [site+b_dir][p_dir]
*  c = [site+p_dir][b_dir]^-1
*  d = [site][p_dir]^-1
*      a
*    *-<-+
*  d |   | b ^
*    +->-+
*      c
*  a = [site][b_dir]^-1
*  d = [site-p_dir][p_dir]^-1
*  c = [site-p_dir][b_dir]
*  b = [site-p_dir+b_dir][p_dir]
*  
*  Here we use Staple
*    [site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
*  + [site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
*  The return, int2(x = linkIndex; |y| - 1 = fieldIndex, 0 if it is not boundary, sign of y is for inverse)  
*/
__device__ void CIndexSquare::_deviceGetPlaquttesAtLink(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiLinkIndex, BYTE st) const
{
    BYTE uiDim = _DC_Dim;

    //for square, dir should equal to dim
    //assert(uiDim == _DC_Dir);

    count = 2 * (uiDim - 1);
    plaqutteLength = 4; //for square

    //For square lattice, we assume dimenssion = number of direction
    UINT uiSiteIndex = uiLinkIndex / uiDim;
    UINT uiLinkDir = uiLinkIndex % uiDim;
    BYTE uiMaxDim = (0 == (st & CIndex::kTime)) ? 3 : 4;
    BYTE uiMinDim = (0 == (st & CIndex::kSpace)) ? 3 : 0;

    //Note, 2D is z, t
    //3D is y, z, t
    //4D is x, y, z, t
    //so i = 4 - dim to 3, (for example, dim = 2, it is 2, 3; dim = 4, it is 0, 1, 2, 3).
    if (uiMinDim < 4 - uiDim)
    {
        uiMinDim = 4 - uiDim;
    }

    //uiLinkIndex is bdir
    //i is pdir
    BYTE iListIndex = 0;
    for (BYTE i = uiMinDim; i < uiMaxDim; ++i)
    {
        if (i != uiLinkDir)
        {
            SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            //=============================================
            //add forward
            //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
            retV[iListIndex] = SIndex(uiSiteIndex, i);
            ++iListIndex;

            SSmallInt4 fsite = _deviceMoveSquareSite(xyzt, i + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 1]);
            retV[iListIndex].m_byDir = uiLinkDir;
            ++iListIndex;

            fsite = _deviceMoveSquareSite(xyzt, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            //=============================================
            //add backward
            //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
            fsite = _deviceMoveSquareSite(xyzt, -(i + 1));
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 3]);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(retV[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
            ++iListIndex;

            fsite = _deviceMoveSquareSite(fsite, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = i;
            ++iListIndex;
        }
    }

    //assert(count * 3 == iListIndex);
}

__device__ void CIndexSquare::_deviceGetPlaquttesAtSite(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiSiteIndex, BYTE st) const
{
    BYTE uiDim = _DC_Dim;

    //for square, dir should equal to dim
    assert(uiDim == _DC_Dir);

    count = uiDim * (uiDim - 1) / 2;
    plaqutteLength = 4; //for square

    BYTE uiMaxDim = (0 == (st & CIndex::kTime)) ? 3 : 4;
    
    BYTE iListIndex = 0;
    for (BYTE uiLink = 0; uiLink < uiMaxDim; ++uiLink)
    {
        BYTE uiMinDim = (0 == (st & CIndex::kSpace)) ? 3 : uiLink + 1;
        for (BYTE uiPlaq = uiMinDim; uiPlaq < uiMaxDim; ++uiPlaq)
        {
            SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            retV[iListIndex] = SIndex(uiSiteIndex, uiLink);
            ++iListIndex;

            SSmallInt4 movedSite = _deviceMoveSquareSite(xyzt, uiLink + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(movedSite, retV[iListIndex - 1]);
            retV[iListIndex].m_byDir = uiPlaq;
            ++iListIndex;

            movedSite = _deviceMoveSquareSite(xyzt, uiPlaq + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(movedSite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = uiLink;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;
        }
    }

    assert(count * 4 == iListIndex);
}

__device__ void CIndexSquare::_deviceGetPlaquttesAtLinkAll(SIndex* retV, UINT uiLinkIndex) const
{
    BYTE uiDim = _DC_Dim;

    UINT uiSiteIndex = uiLinkIndex / uiDim;
    UINT uiLinkDir = uiLinkIndex % uiDim;

    //uiLinkIndex is bdir
    //i is pdir
    BYTE iListIndex = 0;
    for (BYTE i = 0; i < uiDim; ++i)
    {
        if (i != uiLinkDir)
        {
            SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            //=============================================
            //add forward
            //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
            retV[iListIndex] = SIndex(uiSiteIndex, i);
            ++iListIndex;

            SSmallInt4 fsite = _deviceMoveSquareSite(xyzt, i + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 1]);
            retV[iListIndex].m_byDir = uiLinkDir;
            ++iListIndex;

            fsite = _deviceMoveSquareSite(xyzt, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            //=============================================
            //add backward
            //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
            fsite = _deviceMoveSquareSite(xyzt, -(i + 1));
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 3]);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(retV[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
            ++iListIndex;

            fsite = _deviceMoveSquareSite(fsite, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = i;
            ++iListIndex;
        }
    }
}

__device__ void CIndexSquare::_deviceGetPlaquttesAtSiteAll(SIndex* retV, UINT uiSiteIndex) const
{
    BYTE uiDim = _DC_Dim;
    BYTE iListIndex = 0;
    for (BYTE uiLink = 0; uiLink < uiDim; ++uiLink)
    {
        for (BYTE uiPlaq = uiLink + 1; uiPlaq < uiDim; ++uiPlaq)
        {
            SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            retV[iListIndex] = SIndex(uiSiteIndex, uiLink);
            ++iListIndex;

            SSmallInt4 movedSite = _deviceMoveSquareSite(xyzt, uiLink + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(movedSite, retV[iListIndex - 1]);
            retV[iListIndex].m_byDir = uiPlaq;
            ++iListIndex;

            movedSite = _deviceMoveSquareSite(xyzt, uiPlaq + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(movedSite, retV[iListIndex - 2]);
            retV[iListIndex].m_byDir = uiLink;
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            retV[iListIndex].m_byTag = _kDagger;
            ++iListIndex;
        }
    }
}

/**
* virtual, cannot inline
*/
__device__ SIndex CIndexSquare::_deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, SBYTE iWalkDir) const
{
    SSmallInt4 fromSite = __deviceSiteIndexToInt4(uiSiteIndex);
    SIndex sFromSite = SIndex(uiSiteIndex);
    SSmallInt4 siteInt4 = _deviceMoveSquareSite(fromSite, iWalkDir);
    return m_pBoundaryCondition->_devcieGetFermionMappedIndex(uiFieldId, siteInt4, sFromSite);
}

/**
* virtual, cannot inline
*/
__device__ SIndex CIndexSquare::_deviceGaugeIndexWalk(UINT uiSiteIndex, SBYTE iWalkDir) const
{
    SSmallInt4 fromSite = __deviceSiteIndexToInt4(uiSiteIndex);
    SIndex sFromSite = SIndex(uiSiteIndex);
    SSmallInt4 siteInt4 = _deviceMoveSquareSite(fromSite, iWalkDir);
    return m_pBoundaryCondition->_devcieGetMappedIndex(siteInt4, sFromSite);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================