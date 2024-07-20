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

class CLGAPI CFieldBoundaryParent : public CBase
{
public:
    virtual void InitialField(CParameters& param) = 0;
};

/**
* It is more convinient NOT to inhirent from CField.
*/
template<typename deviceData>
class __DLL_EXPORT CFieldBoundary : public CFieldBoundaryParent
{
public:
    CFieldBoundary() 
        : m_byFieldId(0)
        , m_pDeviceData(NULL)
    {
        CCommonKernel<deviceData>::AllocateBuffer(&m_pDeviceData, 8 * _HC_Dir);
    }

    ~CFieldBoundary() 
    {
        CCommonKernel<deviceData>::FreeBuffer(&m_pDeviceData);
    }

    virtual EFieldType GetFieldType() const = 0;
    void InitialField(CParameters& param) override
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
    deviceData* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================