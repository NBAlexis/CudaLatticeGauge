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

__DEFINE_ENUM(ERational,

    //For WD Nf=2, do not use rational
    ER_NoRational,
    //For WD Nf=1, when prepare for MC, use MD which is (D^+D)^{1/2} instead of D^1, when update, use even - psedufermion field
    ER_WDMDRational,
    //For Staggered or WD Nf!=2, for example Nf=1, when prepare for MC, use MD which is (D^+D)^{1/2}, when update, use (D^+D)^{-1/2}
    ER_AllRational,
    )

class CLGAPI CFieldFermion : public CField
{
public:
    CFieldFermion();

    virtual void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            PrepareForHMCS(pGauge);
            return;
        }
        appCrucial(_T("CFieldFermion PrepareForHMC not implemented\n"));
    }

    /**
    * Calculate force can fail due to solver
    */
    virtual UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, CFieldGauge* const* pGaugeForce, CFieldBoson* const* pBosonForce, ESolverPhase ePhase) const
    {
        if (SingleField())
        {
            INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
            return CalculateForceS(gaugeFields[idx], pGaugeForce[idx], ePhase);
        }
        appCrucial(_T("CalculateForce not implemented\n"));
        return FALSE;
    }

    UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* pOtherParameters = NULL) override
    {
        switch (op)
        {
        case EFO_F_D:
            D(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_Ddagger:
            Ddagger(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DD:
            DD(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DDdagger:
            DDdagger(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_InverseD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseD(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields);
        case EFO_F_InverseDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDdagger(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields);
        case EFO_F_InverseDDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDDdagger(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields);
        case EFO_F_InverseDD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDD(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields);
        case EFO_F_RationalD:
            {
                //It seems this was never used
                INT* iRAIndex = (INT*)(pOtherParameters);
                if (NULL == iRAIndex)
                {
                    appCrucial(_T("ApplyOperator, the operator %s with CRatinalApproximation null.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                    return FALSE;
                }
                return RationalApproximation(EFO_F_D, gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, *iRAIndex);
            }
        case EFO_F_D_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_Ddagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DdaggerWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DdaggerWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DD_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DDdagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDdaggerWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDdaggerWithMass(gaugeNum, bosonNum, tensor2Num, pGauge, pBoson, tensor2Fields, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        default:
            appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
            return FALSE;
        }
    }

    virtual void ApplyGamma(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EGammaMatrix eGamma)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            ApplyGammaS(pGauge, eGamma);
            return;
        }
        appCrucial(_T("ApplyGamma not implemented\n"));
    }

    virtual void D(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("D not implemented\n"));
    }

    virtual void Ddagger(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DdaggerS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("Ddagger not implemented\n"));
    }

    virtual void DDdagger(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDdaggerS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDdagger not implemented\n"));
    }
    virtual void DD(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DD not implemented\n"));
    }
    virtual void DWithMass(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DWithMass not implemented\n"));
    }
    virtual void DdaggerWithMass(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DdaggerWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DdaggerWithMass not implemented\n"));
    }
    virtual void DDWithMass(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDWithMass not implemented\n"));
    }
    virtual void DDdaggerWithMass(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDdaggerWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDdaggerWithMass not implemented\n"));
    }

    virtual UBOOL InverseD(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields);
    virtual UBOOL InverseDdagger(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields);
    virtual UBOOL InverseDD(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields);
    virtual UBOOL InverseDDdagger(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields);
    virtual void InitialAsSource(const SFermionBosonSource& sourceData) = 0;

    virtual TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const = 0;
    virtual UBOOL RationalApproximation(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, INT iRationalIndex, UBOOL bAction = FALSE);
    virtual UBOOL RationalApproximationPooled(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, INT iRationalIndex, TArray<CField*>& solutions) const;

    UBOOL IsFermionField() const override { return TRUE; }
    UINT GetSiteCount() const { return m_uiSiteCount; }

    virtual TArray<CField*> CalculateRationalFields(const CFieldGauge* pGauge) const
    {
        appCrucial(_T("CalculateRationalFields not implemented\n"));
        return TArray<CField*>();
    }

    virtual void D_MC(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields)
    {
        RationalApproximation(EFO_F_DDdagger, gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields, m_iMCIndex);
    }

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D_MD(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, UBOOL bAction = TRUE)
    {
        RationalApproximation(EFO_F_DDdagger, gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields, m_iMDIndex, bAction);
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

    virtual void D0OnEvenOrOdd(UBOOL bEven, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            D0OnEvenOrOddS(pGauge, bEven);
            return;
        }
        appCrucial(_T("D0OnEven not implemented\n"));
    }

    /**
    * The naik link is so special, it uses effective gauge at level-1, so we have to preserve an interface for it
    * The epsilon term is also special due to inherent-machinism, we should use component instead of inherent at the first place, but for now, we use this ugly implementation
    */
    virtual void CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const
    {
        appCrucial(_T("CalculateF0AndNaik not implemented\n"));
    }

protected:

#pragma region single field case

    virtual void ApplyGammaS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
    {
        appCrucial(_T("ApplyGammaS not implemented\n"));
    }

    virtual void PrepareForHMCS(const CFieldGauge* pGauge)
    {
        appCrucial(_T("PrepareForHMCS not implemented\n"));
    }

    virtual UBOOL CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const
    {
        appCrucial(_T("CalculateForceS not implemented\n"));
        return FALSE;
    }

    virtual void DS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DS not implemented\n"));
    }

    virtual void DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DdaggerS not implemented\n"));
    }

    virtual void DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DDdaggerS not implemented\n"));
    }

    virtual void DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DDS not implemented\n"));
    }

    virtual void DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DWithMassS not implemented\n"));
    }

    virtual void DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DdaggerWithMassS not implemented\n"));
    }

    virtual void DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DDWithMassS not implemented\n"));
    }

    virtual void DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        appCrucial(_T("DDdaggerWithMassS not implemented\n"));
    }

    virtual void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
    {
        appCrucial(_T("DOperator not implemented\n"));
    }

    virtual void DerivateDOperator(DOUBLE fCoeff, void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
    {
        appCrucial(_T("DerivateDOperator not implemented\n"));
    }

    /**
     * Use to calculate force of rational apporixmation
     * D with only hopping terms (with only terms which has gauge links)
     */
    virtual void D0S(const CField* pGauge)
    {
        appCrucial(_T("D0S not implemented\n"));
    }

    virtual void D0OnEvenOrOddS(const CField* pGauge, UBOOL bEven)
    {
        appCrucial(_T("D0OnEvenOrOddS not implemented\n"));
    }

    virtual DOUBLE EnergyS(const CFieldGauge* pGauge) const = 0;

#pragma endregion

    //UINT m_uiLinkeCount;
    UINT m_uiSiteCount;

    //Multi-GPU (Phase 2): number of halo SITE slots appended after the local sites
    //in m_pDeviceData. 0 on single-GPU / unsplit builds. Mirrors CFieldGauge's
    //m_uiHaloLinkCount; set at allocation time by the concrete subclass.
    UINT m_uiHaloSiteCount;

public:

    //Multi-GPU (Phase 2): halo site capacity appended to the fermion device buffer.
    UINT GetHaloSiteCount() const { return m_uiHaloSiteCount; }

    //Multi-GPU (Phase 2): bytes of one per-site element (sizeof(deviceVector)), so
    //CHaloManager can pack/exchange a fermion halo without knowing the element type.
    //Default 0 means "halo not wired for this type yet" -> RefillHalo safely skips it
    //(same as the pre-Phase-2 behaviour). Concrete bases override with sizeof(element).
    virtual UINT GetSiteElementBytes() const { return 0; }

    /**
    * Make sure gauge smearing is already done, and is put into gaugeFields
    */
    virtual DOUBLE Energy(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields) const
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            return EnergyS(pGauge);
        }
        appCrucial(_T("Energy not implemented\n"));
        return 0.0;
    }

    UBOOL m_bEvenPseudofermion;

    void InitialOtherParameters(CParameters& params) override;

    CCString GetInfos(const CCString& tab) const override;

    // r(x) = x^{1/4} use to prepare for Nf=2
    // r(x) = x^{3/8} use as s quark for Nf=2+1
    // r(x) = (x+dm/x)^{-1/4} use as u,d quark for Nf=2+1
    //CRatinalApproximation m_rMC;
    INT m_iMCIndex;

    // r(x) = x^{-1/2} use to calculate force and action for Nf=2
    // r(x) = x^{-3/4} use to s quark for Nf=2+1 for Nf=2
    // r(x) = (x+dm/x)^{1/2} use as u,d quark for Nf=2+1
    //CRatinalApproximation m_rMD;
    INT m_iMDIndex;

    ERational m_eRational;

    // In simulation, the action will handle effective gauge automatically, 
    // But in measurement, if one need to use effective gauge for D operator, one need this flag to control it
    // OK, now I found out that the measurement will also handle it (manully).
    // UBOOL m_bDoperatorUseEffectiveGauge;

    void CopyParamTo(CField* U) const override
    {
        CField::CopyParamTo(U);
        CFieldFermion* pField = dynamic_cast<CFieldFermion*>(U);

        pField->m_uiSiteCount = m_uiSiteCount;
        pField->m_uiHaloSiteCount = m_uiHaloSiteCount;
        pField->m_bEvenPseudofermion = m_bEvenPseudofermion;
        pField->m_iMCIndex = m_iMCIndex;
        pField->m_iMDIndex = m_iMDIndex;
        pField->m_eRational = m_eRational;
    }
};


__DEFINE_ENUM(EMeasureDiagnal,
    EMD_D,
    EMD_InverseD,
    EMD_Gamma1,
    EMD_Gamma2,
    EMD_Gamma3,
    EMD_Gamma4,
    EMD_Gamma5,
    EMD_Sigma12,
    EMD_Sigma13,
    EMD_Sigma14,
    EMD_Sigma23,
    EMD_Sigma24,
    EMD_Sigma34,
    EMD_Gamma51,
    EMD_Gamma52,
    EMD_Gamma53,
    EMD_Gamma54,

    EMD_Oribital,
    EMD_Spin,

    EMD_Max,
    )

extern CLGAPI void ExportDiagnalWilsonSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields,
    const class CFieldFermionWilsonSquareSU3* pFermion);

extern CLGAPI void ExportDiagnalStaggeredSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields,
    const class CFieldFermionKSSU3* pFermion);


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================