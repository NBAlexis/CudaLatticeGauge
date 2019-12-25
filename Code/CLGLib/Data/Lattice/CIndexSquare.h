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
    const UBOOL bReverse = dir < 0;
    const BYTE uDir = static_cast<BYTE>(bReverse ? ((-dir) - 1) : (dir - 1));
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

    void BakeAllIndexBuffer(class CIndexData* pData) override;
    void BakePlaquttes(class CIndexData* pData, BYTE byFieldId) override;
    void BakeMoveIndex(class CIndexData* pData, BYTE byFieldId) override;
    void BakeEtaMuTable(class CIndexData* pData) override;
    UINT GetPlaqutteCount() const override;
};

#pragma region device Functions

static __device__ __inline__ SSmallInt4 _deviceCoordMoving(const SSmallInt4& sFrom, BYTE i)
{
    const SBYTE offset = i < _DC_Dir ? -1 : 1;
    SSmallInt4 ret = sFrom;
    ret.m_byData4[i < _DC_Dir ? (4 - _DC_Dir + i) : (4 - _DC_Dir * 2 + i)] += offset;
    return ret;
}

static __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& sSite, const UINT* __restrict__ pSmallData);

static __device__ __inline__ UBOOL _deviceIsBondDirichlet(const BYTE* __restrict__ pTable,
    UINT uiBigIdx, BYTE byDir)
{
    return (pTable[uiBigIdx * _DC_Dir + byDir] & _kDirichlet) != 0;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CINDEXSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================