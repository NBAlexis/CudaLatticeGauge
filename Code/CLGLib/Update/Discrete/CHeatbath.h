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

class CLGAPI CCommonData
{
public:

    static UINT m_iDim;
    static UINT m_iDir;
    static STRING m_sGauge;
    static STRING m_sFermion;

    //related with a
    //related with beta
    static FLOAT m_fG;
};

__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================