//=============================================================================
// FILENAME : CCommonData.h
// 
// DESCRIPTION:
// This is the class for the common data
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CCOMMONDATA_H_
#define _CCOMMONDATA_H_

__BEGIN_NAMESPACE

enum EFieldType
{
    EFT_GaugeSU3,
    EFT_Max,
    EFT_ForceDword = 0x7fffffff,
};

class CLGAPI CCommonData
{
public:
    enum { kLatticeDecompose = 3, kMaxDim = 4, kMaxFieldCount = 8, };

    static UINT m_uiDim;
    static UINT m_uiDir;
    static UINT m_uiLatticeLength[CCommonData::kMaxDim];
    static UINT m_uiMaxThread;
    static EFieldType m_eGaugeType;

    //related with a
    //related with beta
    static FLOAT m_fG;

    static void InitialWithDefault()
    {
        m_uiDim = 4;
        m_uiDir = 4;
        m_uiLatticeLength[0] = 8;
        m_uiLatticeLength[1] = 8;
        m_uiLatticeLength[2] = 8;
        m_uiLatticeLength[3] = 16;
        m_uiMaxThread = 1024;
        m_fG = 1.0f;
        m_eGaugeType = EFT_GaugeSU3;
    }
};

__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================