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

#pragma region Gamma KS

static __device__ __inline__ SCHAR _deviceEta2(UINT uiEta, BYTE i, BYTE j)
{
    return ((uiEta >> i) + (uiEta >> j)) & 1;
}

/**
 * eta xyz, eta yzt, eta xyt, ...
 * for 1, 3 there is a minus sign
 * missingDir:
 * 3 - xyz x:1  y:(-1)^x z:(-1)^(x+y)             res: (-1)^y
 * 2 - xyt x:1  y:(-1)^x t:(-1)^(x+y+z)           res: (-1)^(y+z)
 * 0 - yzt y:(-1)^x z:(-1)^(x+y) t:(-1)^(x+y+z)   res: (-1)^(x+z)
 * 1 - xzt x:1  z:(-1)^(x+y) t:(-1)^(x+y+z)       res: (-1)^z
 * 
 * MILC t,x,y,z
 * 3 - xyz x:(-1)^t, y:(-1)^(t+x), z:(-1)^(t+x+y)       res: (-1)^(y+t)
 * 2 - xyt x:(-1)^t, y:(-1)^(t+x), t:1                  res: (-1)^x
 * 1 - xzt x:(-1)^t, z:(-1)^(t+x+y), t:1                res: (-1)^(x+y)
 * 0 - yzt y:(-1)^(t+x), z:(-1)^(t+x+y), t:1            res: (-1)^y
 */
static __device__ __inline__ SCHAR _deviceEta3(const SSmallInt4& sSite, BYTE missingDir)
{
    missingDir += 4U * static_cast<BYTE>(_DC_MILCSTAGGEREDPHASE);

    switch (missingDir)
    {
    case 3:
        return (sSite.y + 1) & 1;
    case 2:
        return (sSite.y + sSite.z) & 1;
    case 0:
        return (sSite.x + sSite.z) & 1;
    case 1:
        return (sSite.z + 1) & 1;

    //MILC convention
    case 7:
        return (sSite.y + sSite.w + 1) & 1;
    case 6:
        return sSite.x & 1;
    case 5:
        return (sSite.x + sSite.y + 1) & 1;
    case 4:
        return sSite.y & 1;
    }
    return 0;
}

/**
* just as same as gamma53 (case 2 of above)
*/
static __device__ __inline__ Real _deviceEta124(const SSmallInt4& sSite)
{
    return _deviceEta3(sSite, 2) ? (F(-1.0)) : (F(1.0));
}

#pragma endregion

/**
* manipulate site
*/
//__device__ __inline__ static
//SSmallInt4 _deviceMoveSquareSite(SSmallInt4 ret, SCHAR dir)
//{
//    const UBOOL bReverse = dir < 0;
//    const BYTE uDir = static_cast<BYTE>(bReverse ? ((-dir) - 1) : (dir - 1));
//    if (0 == uDir)
//    {
//        if (bReverse)
//        {
//            ret.x--;
//        }
//        else 
//        {
//            ret.x++;
//        }
//    }
//    else if (1 == uDir)
//    {
//        if (bReverse)
//        {
//            ret.y--;
//        }
//        else
//        {
//            ret.y++;
//        }
//    }
//    else if (2 == uDir)
//    {
//        if (bReverse)
//        {
//            ret.z--;
//        }
//        else
//        {
//            ret.z++;
//        }
//    }
//    else 
//    {
//        if (bReverse)
//        {
//            ret.w--;
//        }
//        else
//        {
//            ret.w++;
//        }
//    }
//    return ret;
//}

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

    void BakeAllIndexBuffer(class CIndexData* pData) override;
    void BakePlaquttes(class CIndexData* pData, BYTE byFieldId) override;
    void BakeMoveIndex(class CIndexData* pData, BYTE byFieldId) override;
    void BakeEtaMuTable(class CIndexData* pData) override;
    void BakeNaikTable(class CIndexData* pData, BYTE byFieldId) override;
    UINT GetPlaqutteCount(BYTE byFieldId) const override;

};

#pragma region device Functions

//static __device__ __inline__ SSmallInt4 _deviceCoordMoving(const SSmallInt4& sFrom, BYTE i)
//{
//    const SCHAR offset = i < _DC_Dir ? -1 : 1;
//    SSmallInt4 ret = sFrom;
//    ret.m_byData4[i < _DC_Dir ? (4 - _DC_Dir + i) : (4 - _DC_Dir * 2 + i)] += offset;
//    return ret;
//}

static __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& sSite, const UINT* __restrict__ pSmallData);

//static __device__ __inline__ UBOOL _deviceIsBondDirichlet(const SIndex* __restrict__ pTable,
//    UINT uiBigIdx, BYTE byDir)
//{
//    return (pTable[uiBigIdx * _DC_Dir + byDir].m_byTag & _kDirichlet) != 0;
//}



#pragma endregion

#if _CLG_ASSUME_SQUARE_LATTICE
#define plaqLength 4U
#define plaqLengthm1 3U
#define plaqCountPerLink 6U
#define plaqCountPerSite 6U
#define plaqCountAllLink 18U
#define plaqCountAllSite 24U
#endif

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================