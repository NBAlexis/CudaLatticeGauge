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

__BEGIN_NAMESPACE

class CLGAPI CFieldFermion : public CField
{
public:
    CFieldFermion();

    virtual void InitialFieldWithFile(const CCString& , EFieldFileType ) 
    {
        appCrucial(_T("Not supported for fermions"));
    }

    virtual void PrepareForHMC(const CFieldGauge* pGauge) = 0;
    virtual void CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce) = 0;

    virtual UBOOL ApplyOperator(EFieldOperator op, const CField* otherfield)
    {
        switch (op)
        {
        case EFO_F_D:
            D(otherfield);
            return TRUE;
        case EFO_F_Ddagger:
            Ddagger(otherfield);
            return TRUE;
        case EFO_F_DDdagger:
            DDdagger(otherfield);
            return TRUE;
        case EFO_F_InverseD:
            return InverseD(otherfield);
        case EFO_F_InverseDdagger:
            return InverseDdagger(otherfield);
        case EFO_F_InverseDDdagger:
            return InverseDDdagger(otherfield);
        default:
            appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
            return FALSE;
        }
    }

    virtual void D(const CField* pGauge) = 0;
    virtual void Ddagger(const CField* pGauge) = 0;
    virtual void DDdagger(const CField* pGauge) = 0;
    virtual UBOOL InverseD(const CField* pGauge) = 0;
    virtual UBOOL InverseDdagger(const CField* pGauge) = 0;
    virtual UBOOL InverseDDdagger(const CField* pGauge) = 0;

protected:

    UINT m_uiLinkeCount;
    UINT m_uiSiteCount;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================