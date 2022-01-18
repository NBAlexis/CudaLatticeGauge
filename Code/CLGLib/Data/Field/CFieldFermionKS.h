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

    virtual void ApplyGammaKS(const CFieldGauge* pGauge, EGammaMatrix eGamma) = 0;

    //================= test anti-hermitian =========
    virtual UINT TestAntiHermitian(const CFieldGauge* pGauge) const = 0;

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D_MD(const CField* pGauge) = 0;
    virtual void D0(const CField* pGauge) = 0;
    virtual void D_MC(const CField* pGauge) = 0;

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
        pField->m_rMC = m_rMC;
        pField->m_rMD = m_rMD;
        pField->m_bEachSiteEta = m_bEachSiteEta;
    }

    //============================
    //Override these two functions for KS
    virtual void DerivateD0(void* pForce, const void* pGaugeBuffer) const = 0;
    virtual void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const = 0;
    //============================

    /**
     * Do not override me
     */
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, m_f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    /**
     * DerivateDOperator not implemented for KS
     */
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override
    {
        appCrucial(_T("Do not call KS DerivateDOperator!"));
        _FAIL_EXIT;
    }

    #pragma region Help functions to implement higher orders

    virtual void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) = 0;
    virtual void OneLink(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) = 0;
    virtual void OneLinkForce(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const = 0;

    #pragma endregion

    //For test use only!
    void TestSetEtaShift(UBOOL bShift) { m_bEachSiteEta = bShift; }
    UBOOL TestIsEtaShift() const { return m_bEachSiteEta; }

    //For some strange boundary condition
    //Normally, eta_{\mu}(n+\mu)=eta_{\mu}, so set this = FALSE
    UBOOL m_bEachSiteEta;

    Real m_f2am;

protected:

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