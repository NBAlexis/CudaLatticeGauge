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
        , m_uiSiteCount(_HC_Volume)
    {

    }

    virtual void MakeRandomMomentum() = 0;

    UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, 
        EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* pOtherParameters = NULL) override;

    /**
    * Calculate force can fail due to solver
    */
    //virtual UBOOL CalculateForceOnGauge(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const = 0;

    virtual void D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldGauge** pGaugeForce, const CFieldBoson* const* pBoson) = 0;
    //virtual void DD(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    

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
    virtual void FieldMultply(const CFieldBoson* x, UBOOL bConj = TRUE) = 0;

    UINT GetSiteCount() const { return m_uiSiteCount; }

protected:

    UINT m_uiSiteCount;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSON_H_

//=============================================================================
// END OF FILE
//=============================================================================