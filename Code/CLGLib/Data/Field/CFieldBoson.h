//=============================================================================
// FILENAME : CFieldBoson.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOSON_H_
#define _CFIELDBOSON_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldBoson : public CField
{
public:
    CFieldBoson()
        : CField()
    {

    }

    virtual void PrepareForHMC(const CFieldGauge* pGauge) = 0;

    /**
    * Calculate force can fail due to solver
    */
    //virtual UBOOL CalculateForceOnGauge(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const = 0;

    virtual void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DD(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    
    virtual void Square(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;

#pragma region real operators

#pragma endregion

    UBOOL IsGaugeField() const override { return FALSE; }
    UBOOL IsFermionField() const override { return FALSE; }
    UBOOL IsBosonField() const override  { return TRUE; }

    /**
     * For even odd preconditioner
     */
    virtual void WriteEvenSites(const CFieldFermion*, const CFieldGauge*, UBOOL) { appCrucial(_T("Not implemented.\n")); }
    virtual void WriteBackEvenSites(CFieldFermion*, const CFieldGauge*, UBOOL) const { appCrucial(_T("Not implemented.\n")); }

    UINT GetSiteCount() const { return m_uiSiteCount; }

protected:

    UINT m_uiSiteCount;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSON_H_

//=============================================================================
// END OF FILE
//=============================================================================