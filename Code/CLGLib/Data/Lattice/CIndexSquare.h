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
SSmallInt4 _deviceMoveSquareSite(SSmallInt4 ret, SBYTE dir)
{
    UBOOL bReverse = dir < 0;
    BYTE uDir = static_cast<BYTE>(bReverse ? ((-dir) - 1) : (dir - 1));
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

    /**
    * To be discarded..
    * 
    */
    __device__ virtual void _deviceGetPlaquttesAtLink(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiLinkIndex, BYTE st = kSpaceTime) const;
    __device__ virtual void _deviceGetPlaquttesAtSite(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiSiteIndex, BYTE st = kSpaceTime) const;
    __device__ virtual void _deviceGetPlaquttesAtLinkAll(SIndex* retV, UINT uiLinkIndex) const;
    __device__ virtual void _deviceGetPlaquttesAtSiteAll(SIndex* retV, UINT uiSiteIndex) const;

    __device__ virtual void _deviceGetPlaqutteCountLength(BYTE& plaqLength, BYTE& countPerSite, BYTE& countPerLink)
    {
        plaqLength = 4;
        countPerSite = static_cast<BYTE>(_DC_Dim * (_DC_Dim - 1) / 2);
        countPerLink = static_cast<BYTE>(2 * (_DC_Dim - 1));
    }

    __device__ virtual SIndex _deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, SBYTE uiWalkDir) const;
    __device__ virtual SIndex _deviceGaugeIndexWalk(UINT uiSiteIndex, SBYTE uiWalkDir) const;
};

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================