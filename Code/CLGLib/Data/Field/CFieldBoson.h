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
        , m_bConstant(FALSE)
    {

    }

    virtual void MakeRandomMomentum() = 0;

    UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, 
        EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* pOtherParameters = NULL) override;

    /**
    * NOTE: This is minus term in D
    * S = - (Dphi)^2 + c1 phi^2 + c2 phi^4 + ...
    * This is not D, but the minus terms in D
    * (partial phi)^2 = phi^2 - phi(n)[U_{mu}(n)phi(n+mu) + U_{-mu}(n)phi(n-mu)]
    * This is [U_{mu}(n)phi(n+mu) + U_{-mu}(n)phi(n-mu)]
    */
    virtual void D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const = 0;

    virtual UINT CheckHermitian(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) const = 0;
    virtual void InitialAsSource(const SFermionBosonSource& sourceData) = 0;

protected:

    virtual void DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;

public:
    
    
    void InitialOtherParameters(CParameters& param) override
    {
        CField::InitialOtherParameters(param);

        INT iConst = 0;
        param.FetchValueINT(_T("Constant"), iConst);
        m_bConstant = (0 != iConst);
    }

#pragma region real operators

#pragma endregion

    UBOOL IsBosonField() const override { return TRUE; }
    virtual UINT VectorN() const = 0;
    virtual UINT FloatN() const = 0;
    
    UINT GetSiteCount() const { return m_uiSiteCount; }

    virtual TArray<DOUBLE> Sum() const = 0;

    void CopyTo(CField* U) const override
    {
        CField::CopyTo(U);

        CFieldBoson* pOther = dynamic_cast<CFieldBoson*>(U);
        pOther->m_uiSiteCount = m_uiSiteCount;
        pOther->m_bConstant = m_bConstant;
    }

protected:

    UINT m_uiSiteCount;

public:

    //Use for debug
    UBOOL m_bConstant;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSON_H_

//=============================================================================
// END OF FILE
//=============================================================================