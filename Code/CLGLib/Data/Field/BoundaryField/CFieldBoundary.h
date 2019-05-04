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
    CFieldBoundary() {}
    ~CFieldBoundary() {}

    virtual EFieldType GetFieldType() const = 0;
    virtual void InitialField(CParameters& param) = 0;
    virtual CCString GetInfos(const CCString &tab) const = 0;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================