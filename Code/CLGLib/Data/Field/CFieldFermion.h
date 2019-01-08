//=============================================================================
// FILENAME : CFieldFermion.h
// 
// DESCRIPTION:
// This is the class for all fermion fields
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDFERMION_H_
#define _CFIELDFERMION_H_

#define intokernal_fermion \
intokernal \
UN_USE(uiDir);

__BEGIN_NAMESPACE

class CLGAPI CFieldFermion : public CField
{
public:
    CFieldFermion();

    virtual void PrepareForHMC(const CFieldGauge* pGauge) = 0;
    virtual void CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce) = 0;

    virtual void ApplyOperator(EFieldOperator op, const CField* otherfield)
    {
        switch (op)
        {
        case EFO_F_D:
            D(otherfield);
            break;
        case EFO_F_Ddagger:
            Ddagger(otherfield);
            break;
        case EFO_F_DDdagger:
            DDdagger(otherfield);
            break;
        case EFO_F_InverseDDdagger:
            InverseDDdagger(otherfield);
            break;
        default:
            appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
            break;
        }
    }

    virtual void D(const CField* pGauge) = 0;
    virtual void Ddagger(const CField* pGauge) = 0;
    virtual void DDdagger(const CField* pGauge) = 0;
    virtual void InverseDDdagger(const CField* pGauge) = 0;

protected:

    UINT m_uiLinkeCount;
    UINT m_uiSiteCount;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================