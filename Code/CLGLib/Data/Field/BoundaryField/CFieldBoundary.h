//=============================================================================
// FILENAME : CFieldBoundary.h
// 
// DESCRIPTION:
//
// To implement Dirichlet boundary conditions
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARY_H_
#define _CFIELDBOUNDARY_H_

__BEGIN_NAMESPACE

/**
* It is more convinient NOT to inhirent from CField.
*/
class CLGAPI CFieldBoundary : public CBase
{
public:
    CFieldBoundary() : m_byFieldId(0) {}
    ~CFieldBoundary() {}

    virtual EFieldType GetFieldType() const = 0;
    virtual void InitialField(CParameters& param)
    {
        INT iValue = 1;
        param.FetchValueINT(_T("FieldId"), iValue);
        m_byFieldId = static_cast<BYTE>(iValue);
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString ret = CBase::GetInfos(tab);
        ret = ret + tab + _T("FieldId : ") + appToString(m_byFieldId) + _T("\n");
        return ret;
    }

    BYTE m_byFieldId;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================