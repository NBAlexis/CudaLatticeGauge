//=============================================================================
// FILENAME : CLatticeData.h
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CLATTICEDATA_H_
#define _CLATTICEDATA_H_

__BEGIN_NAMESPACE

class CLGAPI CLatticeData
{
public:

    /**
    * Initial with CommonData
    */
    CLatticeData();

    /**
    * One thread deal with only data[x / m_pLatticeDecompose[0], y / m_pLatticeDecompose[1], z / m_pLatticeDecompose[2]]
    * with m_pLatticeDecompose[3] * m_pLatticeDecompose[4] * m_pLatticeDecompose[5] blocks
    * For 1D, m_pLatticeDecompose[1] = m_pLatticeDecompose[2] = 1
    * For 2D, m_pLatticeDecompose[2] = 1
    */
    UINT m_uiVolumn;
    UINT m_uiDim;
    UINT m_uiDir;
    UINT m_uiLatticeLength[CCommonData::kMaxDim];
    UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];

    /*
    * SU3(x=(x,y,z,t))_{n=a*3+b}=
    * m_pData[(
        (x*m_uiLatticeLength[1]*m_uiLatticeLength[2]*m_uiLatticeLength[3] 
       + y*m_uiLatticeLength[2]*m_uiLatticeLength[3] 
       + z*m_uiLatticeLength[3] 
       + t)
       * m_uiDir + dir) * elementCount + n]
    *
    * m_uiLatticeMultipy[0] = m_uiLatticeLength[1]*m_uiLatticeLength[2]*m_uiLatticeLength[3] * dir * elementcount
    * m_uiLatticeMultipy[1] = m_uiLatticeLength[2]*m_uiLatticeLength[3] * dir * elementcount
    * m_uiLatticeMultipy[2] = m_uiLatticeLength[3] * dir * elementcount
    * m_uiLatticeMultipy[3] = dir * elementcount
    * m_uiLatticeMultipy[4] = elementcount
    * for field on set, dir = 1
    */
    UINT m_uiLatticeMultipy[CCommonData::kMaxDim + 1];

    //Feel free to set and get it
    class CIndex* m_pIndex;

    STRING m_sFields[CCommonData::kMaxFieldCount];
    class CField* m_pFields[CCommonData::kMaxFieldCount];
};

__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================