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
#pragma once

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
    {
        
    }

    ~CFieldFermionKS()
    {

    }

    void InitialOtherParameters(CParameters& params) override;

    //================= test anti-hermitian =========
    virtual UINT TestAntiHermitian(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields) const
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
    void SetMass(Real f2am)
    {
        m_f2am = f2am;
    }

    Real GetMass() const { return m_f2am; }

    void CopyParamTo(CField* U) const override
    {
        CFieldFermion::CopyParamTo(U);
        CFieldFermionKS* pField = dynamic_cast<CFieldFermionKS*>(U);
        pField->m_f2am = m_f2am;
        pField->m_bDiagonalMass = m_bDiagonalMass;
        pField->m_bEachSiteEta = m_bEachSiteEta;
    }

    DOUBLE EnergyS(const CFieldGauge* pGauge) const override;

public:

    #pragma region Help functions to implement higher orders

    virtual void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const  = 0;

    virtual void OneLink(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, void* pTarget, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
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
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const
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

    //================= test anti-hermitian =========
    virtual UINT TestAntiHermitianS(const CFieldGauge* pGauge) const
    {
        appCrucial(_T("TestAntiHermitianS not implemented\n"));
        return 1;
    }

    //These are truely D or InverseD etc.

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

    virtual void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
    {
        appCrucial(_T("D0OperatorKSOnEvenOrOdd not implemented\n"));
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

    virtual void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const
    {
        appCrucial(_T("CalculateForceEvenOddS not implemented\n"));
    }

    /**
    * calculate f0 by using one shift solver
    */
    virtual void CalculateForceEvenOddS_SingleTermOfRational(const CFieldGauge* pGauge, const CFieldFermionKS* phi, CFieldGauge* pForce, Real fCoef, INT iRationalTermIndex) const
    {
        appCrucial(_T("CalculateForceEvenOddS_SingleTermOfRational not implemented\n"));
    }

    /**
    * calculate f0 by using one shift solver
    */
    virtual void CalculateForceS_SingleTermOfRational(const CFieldGauge* pGauge, const CFieldFermionKS* phi, const CFieldFermionKS* phid, CFieldGauge* pForce, Real fCoef, INT iRationalTermIndex) const
    {
        appCrucial(_T("CalculateForceS_SingleTermOfRational not implemented\n"));
    }

public:

    virtual void OneLinkS(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
    {
        appCrucial(_T("OneLinkS not implemented\n"));
    }

    virtual void OneLinkForceS(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const
    {
        appCrucial(_T("OneLinkForceS not implemented\n"));
    }

#pragma endregion

public:


    //Real* m_pMDNumerator;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKS_H_

//=============================================================================
// END OF FILE
//=============================================================================