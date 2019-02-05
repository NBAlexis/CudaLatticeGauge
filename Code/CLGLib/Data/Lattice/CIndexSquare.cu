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
__device__ void CIndexSquare::_deviceGetPlaquttesAtLink(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiLinkIndex, UINT st) const
{
    UINT uiDim = _DC_Dim;

    //for square, dir should equal to dim
    assert(uiDim == _DC_Dir);

    count = 2 * (uiDim - 1);
    plaqutteLength = 4; //for square

    //For square lattice, we assume dimenssion = number of direction
    UINT uiSiteIndex = uiLinkIndex / uiDim;
    UINT uiLinkDir = uiLinkIndex % uiDim;
    UINT uiMaxDim = (0 == (st & CIndex::kTime)) ? 3 : 4;
    UINT uiMinDim = (0 == (st & CIndex::kSpace)) ? 3 : 0;

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
    UINT iListIndex = 0;
    for (int i = uiMinDim; i < uiMaxDim; ++i)
    {
        if (i != uiLinkDir)
        {
            int4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            //=============================================
            //add forward
            //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
            retV[iListIndex] = SIndex(uiSiteIndex, i);
            ++iListIndex;

            int4 fsite3 = _deviceMoveSquareSite(xyzt, i + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite3, xyzt);
            retV[iListIndex].m_byDir = uiLinkDir;
            ++iListIndex;

            int4 fsite2 = _deviceMoveSquareSite(xyzt, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite2, xyzt);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = SIndex::kDagger;
            ++iListIndex;

            //=============================================
            //add backward
            //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
            int4 bsite2 = _deviceMoveSquareSite(xyzt, -(i + 1));
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(bsite2, xyzt);
            retV[iListIndex].m_byDir = i;
            retV[iListIndex].m_byTag = SIndex::kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(retV[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
            ++iListIndex;

            int4 bsite4 = _deviceMoveSquareSite(bsite2, uiLinkDir + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(bsite4, bsite2);
            retV[iListIndex].m_byDir = i;
            ++iListIndex;
        }
    }

    assert(count * 3 == iListIndex);
}

__device__ void CIndexSquare::_deviceGetPlaquttesAtSite(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiSiteIndex, UINT st) const
{
    UINT uiDim = _DC_Dim;

    //for square, dir should equal to dim
    assert(uiDim == _DC_Dir);

    count = uiDim * (uiDim - 1) / 2;
    plaqutteLength = 4; //for square

    UINT uiMaxDim = (0 == (st & CIndex::kTime)) ? 3 : 4;
    
    UINT iListIndex = 0;
    for (UINT uiLink = 0; uiLink < uiMaxDim; ++uiLink)
    {
        UINT uiMinDim = (0 == (st & CIndex::kSpace)) ? 3 : uiLink + 1;
        for (UINT uiPlaq = uiMinDim; uiPlaq < uiMaxDim; ++uiPlaq)
        {
            int4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

            retV[iListIndex] = SIndex(uiSiteIndex, uiLink);
            ++iListIndex;

            int4 fsite1 = _deviceMoveSquareSite(xyzt, uiLink + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite1, xyzt);
            retV[iListIndex].m_byDir = uiPlaq;
            ++iListIndex;

            int4 fsite2 = _deviceMoveSquareSite(xyzt, uiPlaq + 1);
            retV[iListIndex] = m_pBoundaryCondition->_devcieGetMappedIndex(fsite2, xyzt);
            retV[iListIndex].m_byDir = uiLink;
            retV[iListIndex].m_byTag = SIndex::kDagger;
            ++iListIndex;

            retV[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            retV[iListIndex].m_byTag = SIndex::kDagger;
            ++iListIndex;
        }
    }

    assert(count * 4 == iListIndex);
}

/**
* virtual, cannot inline
*/
__device__ SIndex CIndexSquare::_deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, INT iWalkDir) const
{
    int4 fromSite = __deviceSiteIndexToInt4(uiSiteIndex);
    int4 siteInt4 = _deviceMoveSquareSite(fromSite, iWalkDir);
    return m_pBoundaryCondition->_devcieGetFermionMappedIndex(uiFieldId, siteInt4, fromSite);
}

/**
* virtual, cannot inline
*/
__device__ SIndex CIndexSquare::_deviceGaugeIndexWalk(UINT uiSiteIndex, INT iWalkDir) const
{
    int4 fromSite = __deviceSiteIndexToInt4(uiSiteIndex);
    int4 siteInt4 = _deviceMoveSquareSite(fromSite, iWalkDir);
    return m_pBoundaryCondition->_devcieGetMappedIndex(siteInt4, fromSite);
}

#pragma TODO(optimiztion)

#if 0
/**
* Using _deviceGetPlaquttesAtLink, will get all plaquttes connecting to a link, it is mu neq nu
* However, when calculating H, it is reductent, because
* Re(F(mu,nu)) = Re(F(nu,mu))
* So, we only need to sum over mu < nu, not mu neq nu
* Note that, beta is also redefined!
*/
__device__ void CIndexSquare::_deviceGetUniquePlaquttesAtLink(int2* retV, UINT& count, UINT& plaqutteLength, UINT uiLinkIndex) const
{
    UINT uiDim = _DC_Dim;
    //UINT* length = pLattice->m_uiLatticeLength;
    //UINT* mult = pLattice->m_uiLatticeMultipy;

    //for square, dir should equal to dim
    assert(uiDim == _DC_Dir);

    //For square lattice, we assume dimenssion = number of direction
    UINT uiSiteIndex = uiLinkIndex / uiDim;
    UINT uiLinkDir = uiLinkIndex % uiDim;

    count = (3 - uiLinkDir);
    plaqutteLength = 4; //for square

    //uiLinkIndex is bdir
    //i is pdir
    UINT iListIndex = 0;
    for (int i = uiLinkDir + 1; i < 4; ++i)
    {
        int4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);

        //=============================================
        //add forward
        //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
        retV[iListIndex].x = uiSiteIndex + i;
        retV[iListIndex].y = 1;
        ++iListIndex;

        int4 fsite3 = _deviceMoveSquareSite(xyzt, i + 1);
        uint2 fsiteIndex3 = m_pBoundaryCondition->_devcieGetMappedIndex(fsite3, xyzt);
        retV[iListIndex].x = fsiteIndex3.x + uiLinkDir;
        retV[iListIndex].y = (fsiteIndex3.y + 1);
        ++iListIndex;

        int4 fsite2 = _deviceMoveSquareSite(xyzt, uiLinkDir + 1);
        uint2 fsiteIndex2 = m_pBoundaryCondition->_devcieGetMappedIndex(fsite2, xyzt);
        retV[iListIndex].x = fsiteIndex2.x + i;
        retV[iListIndex].y = -(fsiteIndex2.y + 1);
        ++iListIndex;

        //=============================================
        //add backward
        //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
        int4 bsite2 = _deviceMoveSquareSite(xyzt, -(i + 1));
        uint2 bsiteIndex2 = m_pBoundaryCondition->_devcieGetMappedIndex(bsite2, xyzt);
        retV[iListIndex].x = bsiteIndex2.x + i;
        retV[iListIndex].y = -(bsiteIndex2.y + 1);
        ++iListIndex;

        retV[iListIndex].x = bsiteIndex2.x + uiLinkDir;
        retV[iListIndex].y = bsiteIndex2.y + 1;
        ++iListIndex;

        int4 bsite4 = _deviceMoveSquareSite(bsite2, uiLinkDir + 1);
        uint2 bsiteIndex4 = m_pBoundaryCondition->_devcieGetMappedIndex(bsite4, bsite2);
        retV[iListIndex].x = bsiteIndex4.x + i;
        retV[iListIndex].y = bsiteIndex2.y + 1;
        ++iListIndex;
    }
}
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================