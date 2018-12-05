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
__device__ int2* CIndexSquare::GetPlaquttesAtLink(UINT& count, UINT& plaqutteLength, UINT uiDim, UINT uiLinkIndex, const UINT* length, const UINT* mult, UINT st)
{
    int2* retV;
    count = 2 * (uiDim - 1);
    plaqutteLength = 4; //for square
    cudaMalloc((void **)&retV, sizeof(int2) * count * (plaqutteLength - 1));

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
    UINT elementCount = mult[4]; //for example, for SU3 it is 9
    UINT iListIndex = 0;
    for (int i = uiMinDim; i < uiMaxDim; ++i)
    {
        if (i != uiLinkDir)
        {
            int4 xyzt = GetXYZTSquare(uiSiteIndex, mult);

            //=============================================
            //add forward
            //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
            retV[iListIndex].x = uiSiteIndex + i * elementCount;
            retV[iListIndex].y = 1;
            ++iListIndex;

            int4 fsite3 = MoveSquareSite(xyzt, i + 1);
            uint2 fsiteIndex3 = m_pBoundaryCondition->GetMappedIndex(fsite3, xyzt, length, mult);
            retV[iListIndex].x = fsiteIndex3.x + uiLinkDir * elementCount;
            retV[iListIndex].y = (fsiteIndex3.y + 1);
            ++iListIndex;

            int4 fsite2 = MoveSquareSite(xyzt, uiLinkDir + 1);
            uint2 fsiteIndex2 = m_pBoundaryCondition->GetMappedIndex(fsite2, xyzt, length, mult);
            retV[iListIndex].x = fsiteIndex2.x + i * elementCount;
            retV[iListIndex].y = -(fsiteIndex2.y + 1);
            ++iListIndex;

            //=============================================
            //add backward
            //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
            int4 bsite2 = MoveSquareSite(xyzt, -(i + 1));
            uint2 bsiteIndex2 = m_pBoundaryCondition->GetMappedIndex(bsite2, xyzt, length, mult);
            retV[iListIndex].x = bsiteIndex2.x + i * elementCount;
            retV[iListIndex].y = -(bsiteIndex2.y + 1);
            ++iListIndex;

            retV[iListIndex].x = bsiteIndex2.x + uiLinkDir * elementCount;
            retV[iListIndex].y = bsiteIndex2.y + 1;
            ++iListIndex;

            int4 bsite4 = MoveSquareSite(bsite2, uiLinkDir + 1);
            uint2 bsiteIndex4 = m_pBoundaryCondition->GetMappedIndex(bsite4, bsite2, length, mult);
            retV[iListIndex].x = bsiteIndex4.x + i * elementCount;
            retV[iListIndex].y = bsiteIndex2.y + 1;
            ++iListIndex;
        }
    }
    return retV;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================