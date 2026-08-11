//=============================================================================
// FILENAME : CFieldTensor2.h
//
// DESCRIPTION:
// This is the class for all tensor2 fields.
// A tensor2 field is a plaquette field, every site has 6 elements,
// ordered as xy, xz, xt, yz, yt, zt (the same order as _plaq_idx),
// stored as data[plaqutteIndex * siteCount + siteIndex] (component-major).
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#ifndef _CFIELDTENSOR2_H_
#define _CFIELDTENSOR2_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldTensor2 : public CField
{
public:
    CFieldTensor2()
        : CField()
        , m_uiSiteCount(_HC_Volume)
    {

    }

    UBOOL ApplyOperator(EFieldOperator, INT, INT, INT, const CFieldGauge* const*, const CFieldBoson* const*, const CFieldTensor2* const*, EOperatorCoefficientType, Real, Real, void*) override
    {
        appCrucial("CFieldTensor2: Do Operator implimented yet\n");
        return FALSE;
    }

    UBOOL IsTensor2Field() const override { return TRUE; }

    /**
    * float (Real) count per site, it is 6 x element dimension
    */
    virtual UINT FloatN() const = 0;

    UINT GetSiteCount() const { return m_uiSiteCount; }

    /**
    * number of plaquettes of one site, xy, xz, xt, yz, yt, zt
    */
    static UINT PlaqutteCountPerSite() { return _HC_Dir * (_HC_Dir - 1) / 2; }

    /**
    * total element count = site count x plaquette count per site
    */
    UINT GetElementCount() const { return m_uiSiteCount * PlaqutteCountPerSite(); }

    void CopyParamTo(CField* U) const override
    {
        CField::CopyParamTo(U);

        CFieldTensor2* pOther = dynamic_cast<CFieldTensor2*>(U);
        pOther->m_uiSiteCount = m_uiSiteCount;
    }

protected:

    UINT m_uiSiteCount;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDTENSOR2_H_

//=============================================================================
// END OF FILE
//=============================================================================
