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
* manipulate site
*/
__device__ __inline__ static
int4 _deviceMoveSquareSite(int4 ret, INT dir)
{
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
    __device__ CIndexSquare(class deviceBoundaryCondition * devicePtr) : CIndex(devicePtr) { ; }

    __device__ virtual void _deviceGetPlaquttesAtLink(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiLinkIndex, UINT st = kSpaceTime) const;
    __device__ virtual void _deviceGetPlaquttesAtSite(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiSiteIndex, UINT st = kSpaceTime) const;

    __device__ virtual SIndex _deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, INT uiWalkDir) const;
    __device__ virtual SIndex _deviceGaugeIndexWalk(UINT uiSiteIndex, INT uiWalkDir) const;
};

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================