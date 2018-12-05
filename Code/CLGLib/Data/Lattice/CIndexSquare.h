//=============================================================================
// FILENAME : CIndexSquare.h
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _CINDEXSQUARE_H_
#define _CINDEXSQUARE_H_

__BEGIN_NAMESPACE

/**
*
*/
__device__ static __inline__ int4 GetXYZTSquare(UINT siteIndex, const UINT* mult)
{
    int4 xyzt;
    xyzt.x = siteIndex / mult[0];
    xyzt.y = (siteIndex % mult[0]) / mult[1];
    xyzt.z = (siteIndex % mult[1]) / mult[2];
    xyzt.w = (siteIndex % mult[2]) / mult[3];
    return xyzt;
}

/**
*
*/
__device__ static __inline__ int4 GetXYZTSquare(UINT linkIndex, UINT dirCount, const UINT* mult)
{
    return GetXYZTSquare(linkIndex / dirCount, mult);
}

/**
* manipulate site
*/
__device__ static __inline__ int4 MoveSquareSite(const int4 &in, INT dir)
{
    int4 ret = in;
    UBOOL bReverse = dir < 0;
    UINT uDir = bReverse ? ((-dir) - 1) : (dir - 1);
    if (0 == uDir)
    {
        if (bReverse)
        {
            ret.x--;
        }
        else 
        {
            ret.x++;
        }
    }
    else if (1 == uDir)
    {
        if (bReverse)
        {
            ret.y--;
        }
        else
        {
            ret.y++;
        }
    }
    else if (2 == uDir)
    {
        if (bReverse)
        {
            ret.z--;
        }
        else
        {
            ret.z++;
        }
    }
    else 
    {
        if (bReverse)
        {
            ret.w--;
        }
        else
        {
            ret.w++;
        }
    }
    return ret;
}

class CLGAPI CIndexSquare : public CIndex
{

public:

    CIndexSquare(class CLatticeData * pOwner) : CIndex(pOwner) { ; }
    virtual __device__ int2* GetPlaquttesAtLink(UINT& count, UINT& plaqutteLength, UINT uiDim, UINT uiLinkIndex, const UINT* length, const UINT* mult, UINT st = kSpaceTime);

};

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================