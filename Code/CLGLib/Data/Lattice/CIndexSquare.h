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

__CLG_REGISTER_HELPER_HEADER(CIndexSquare)

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
    __CLGDECLARE_CLASS(CIndexSquare)

public:
    CIndexSquare() : CIndex() { }

    /**
    * To bake the index array, the volumn is \prod _i (li + 2 * depth)
    * So we need to re calculate the thread decompose.
    * For simplicity, we just decompse using threadIdx.x and blockIdx.x
    * The return value is thread per block
    */
    static UINT GetDecompose(UINT volumn);

    virtual void BakeAllIndexBuffer(class CIndexData* pData);
    virtual void BakePlaquttes(class CIndexData* pData, BYTE byFieldId);
    virtual void BakeMoveIndex(class CIndexData* pData, BYTE byFieldId);
    virtual UINT GetPlaqutteCount() const;
};

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================