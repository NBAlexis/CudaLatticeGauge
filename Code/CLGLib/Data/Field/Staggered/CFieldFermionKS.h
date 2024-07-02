//=============================================================================
// FILENAME : CFieldFermionKSSU3.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKS_H_
#define _CFIELDFERMIONKS_H_

__BEGIN_NAMESPACE


class CLGAPI CFieldFermionKS : public CFieldFermion
{
public:

    CFieldFermionKS()
        : CFieldFermion()
        , m_bEachSiteEta(FALSE)
        , m_f2am(F(0.01))
        , m_bDiagonalMass(FALSE)
        , m_pMDNumerator(NULL)
    {
        
    }

    enum { _kKSLinkLength = 3 };
    void InitialOtherParameters(CParameters& params) override;
    void Zero() override { InitialField(EFIT_Zero); }
    void Identity() override
    {
        appCrucial(_T("Not supported for CFieldFermionKS!"));
    }

    virtual void ApplyGammaKS(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EGammaMatrix eGamma)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            ApplyGammaKSS(pGauge, eGamma);
            return;
        }
        appCrucial(_T("ApplyGammaKS not implemented\n"));
    }

    //================= test anti-hermitian =========
    virtual UINT TestAntiHermitian(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) const
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            return TestAntiHermitianS(pGauge);
        }
        appCrucial(_T("TestAntiHermitian not implemented\n"));
        return 1;
    }

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D_MD(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
    {
        RationalApproximation(EFO_F_DDdagger, gaugeNum, bosonNum, gaugeFields, pBoson, &m_rMD);
    }

    virtual void D0(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            D0S(pGauge);
            return;
        }
        appCrucial(_T("D0 not implemented\n"));
    }

    virtual void D_MC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
    {
        RationalApproximation(EFO_F_DDdagger, gaugeNum, bosonNum, gaugeFields, pBoson, &m_rMC);
    }

    void SetMass(Real f2am)
    {
        m_f2am = f2am;
    }

    Real GetMass() const { return m_f2am; }

    void CopyTo(CField* U) const override
    {
        CField::CopyTo(U);
        CFieldFermionKS* pField = dynamic_cast<CFieldFermionKS*>(U);
        pField->m_f2am = m_f2am;
        pField->m_bDiagonalMass = m_bDiagonalMass;
        pField->m_rMC = m_rMC;
        pField->m_rMD = m_rMD;
        pField->m_bEachSiteEta = m_bEachSiteEta;
    }

public:

    #pragma region Help functions to implement higher orders

    virtual void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) = 0;

    virtual void OneLink(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            OneLinkS(pGauge->GetData(), pGauge->m_byFieldId, pTarget, fCoefficient, pDevicePath, pathLength, byEtaIdx, bDagger, eOCT, fRealCoeff, cCmpCoeff);
            return;
        }
        appCrucial(_T("OneLink not implemented\n"));
    }

    virtual void OneLinkForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, CFieldGauge* const* pGaugeForce, CFieldBoson* const* pBosonForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const
    {
        if (SingleField())
        {
            INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
            const CFieldGauge* pGauge = gaugeFields[idx];
            CFieldGauge* pForce = pGaugeForce[idx];
            OneLinkForceS(pGauge->GetData(), pGauge->m_byFieldId, pForce->GetData(), fCoefficient, pDevicePath, pathLength, byEtaIdx);
            return;
        }
        appCrucial(_T("OneLinkForce not implemented\n"));
    }

    #pragma endregion

    //For test use only!
    void TestSetEtaShift(UBOOL bShift) { m_bEachSiteEta = bShift; }
    UBOOL TestIsEtaShift() const { return m_bEachSiteEta; }
    virtual void PrepareForHMCOnlyRandomize() = 0;
    virtual void PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) = 0;

    CCString GetInfos(const CCString& tab) const override;

    //For some strange boundary condition
    //Normally, eta_{\mu}(n+\mu)=eta_{\mu}, so set this = FALSE
    UBOOL m_bEachSiteEta;

    Real m_f2am;

    //in case mass term is not a number
    UBOOL m_bDiagonalMass;

protected:

#pragma region single field case

    virtual void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
    {
        appCrucial(_T("ApplyGammaKSS not implemented\n"));
    }

    //================= test anti-hermitian =========
    virtual UINT TestAntiHermitianS(const CFieldGauge* pGauge) const
    {
        appCrucial(_T("TestAntiHermitianS not implemented\n"));
        return 1;
    }

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D0S(const CField* pGauge)
    {
        appCrucial(_T("D0S not implemented\n"));
    }

    //============================
    //Override these two functions for KS
    virtual void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
    {
        appCrucial(_T("DerivateD0 not implemented\n"));
    }

    virtual void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
    {
        appCrucial(_T("DOperatorKS not implemented\n"));
    }

    //============================

    /**
     * Do not override me
     */
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, m_f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    virtual void OneLinkS(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
    {
        appCrucial(_T("OneLinkS not implemented\n"));
    }

    virtual void OneLinkForceS(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const
    {
        appCrucial(_T("OneLinkForceS not implemented\n"));
    }

#pragma endregion

    // r(x) = x^{1/4} use to prepare for Nf=2
    // r(x) = x^{3/8} use as s quark for Nf=2+1
    // r(x) = (x+dm/x)^{-1/4} use as u,d quark for Nf=2+1
    CRatinalApproximation m_rMC;

    // r(x) = x^{-1/2} use to calculate force and action for Nf=2
    // r(x) = x^{-3/4} use to s quark for Nf=2+1 for Nf=2
    // r(x) = (x+dm/x)^{1/2} use as u,d quark for Nf=2+1
    CRatinalApproximation m_rMD;

    Real* m_pMDNumerator;
};

#pragma region Some help functions to implement higher orders

/**
 * Same as CFieldFermionKSSU3R
 * full is a list of path directions with length = iLength
 * it will be divided into two list, where l is full[0, iSep], r is (full[iSep, iLength])^dagger
 * l, r should be allocated on device
 */
static __device__ __inline__ void _deviceSeperate(const INT* __restrict__ full, INT iSep, UINT iLength, INT* l, INT* r, BYTE& LL, BYTE& RL)
{
    LL = static_cast<BYTE>(iSep);
    RL = static_cast<BYTE>(iLength - iSep);

    for (INT i = 0; i < LL; ++i)
    {
        l[i] = -full[iSep - i - 1];
    }

    for (INT i = 0; i < RL; ++i)
    {
        r[i] = full[iSep + i];
    }
}

static __device__ __inline__ void _devicePathDagger(const INT* __restrict__ path, INT* res, UINT iLength)
{
    for (UINT i = 0; i < iLength; ++i)
    {
        res[i] = -path[iLength - i - 1];
    }
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKS_H_

//=============================================================================
// END OF FILE
//=============================================================================