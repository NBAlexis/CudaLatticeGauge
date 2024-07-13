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

__DEFINE_ENUM(EFermionBosonSource,
    EFS_Point,
    EFS_Wall,
    EFS_MomentumWall,
)

//this is a host only structure
struct SFermionBosonSource
{
    EFermionBosonSource m_eSourceType;
    SSmallInt4 m_sSourcePoint;
    //For staggered fermion, spin index is 0-7, for shift of x,y,z
    BYTE m_bySpinIndex;
    BYTE m_byColorIndex;
    CLGComplex m_cOtherParameters1;
    CLGComplex m_cOtherParameters2;
    CLGComplex m_cOtherParameters3;
    CLGComplex m_cOtherParameters4;
};

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

    UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* pOtherParameters = NULL) override
    {
        switch (op)
        {
        case EFO_F_D:
            D(gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_Ddagger:
            Ddagger(gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DD:
            DD(gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DDdagger:
            DDdagger(gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_InverseD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseD(gaugeNum, bosonNum, pGauge, pBoson);
        case EFO_F_InverseDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDdagger(gaugeNum, bosonNum, pGauge, pBoson);
        case EFO_F_InverseDDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDDdagger(gaugeNum, bosonNum, pGauge, pBoson);
        case EFO_F_InverseDD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDDdagger(gaugeNum, bosonNum, pGauge, pBoson);
        case EFO_F_RationalD:
            {
                class CRatinalApproximation* pRA = (class CRatinalApproximation*)(pOtherParameters);
                if (NULL == pRA)
                {
                    appCrucial(_T("ApplyOperator, the operator %s with CRatinalApproximation null.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                    return FALSE;
                }
                return RationalApproximation(EFO_F_D, gaugeNum, bosonNum, pGauge, pBoson, pRA);
            }
        case EFO_F_D_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DWithMass(gaugeNum, bosonNum, pGauge, pBoson, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DWithMass(gaugeNum, bosonNum, pGauge, pBoson, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_Ddagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DdaggerWithMass(gaugeNum, bosonNum, pGauge, pBoson, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DdaggerWithMass(gaugeNum, bosonNum, pGauge, pBoson, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DD_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDWithMass(gaugeNum, bosonNum, pGauge, pBoson, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDWithMass(gaugeNum, bosonNum, pGauge, pBoson, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DDdagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDdaggerWithMass(gaugeNum, bosonNum, pGauge, pBoson, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDdaggerWithMass(gaugeNum, bosonNum, pGauge, pBoson, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        default:
            appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
            return FALSE;
        }
    }

    virtual void ApplyGamma(EGammaMatrix eGamma) = 0;

    virtual void D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("D not implemented\n"));
    }

    virtual void Ddagger(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DdaggerS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("Ddagger not implemented\n"));
    }

    virtual void DDdagger(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDdaggerS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDdagger not implemented\n"));
    }
    virtual void DD(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDS(pGauge, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DD not implemented\n"));
    }
    virtual void DWithMass(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DWithMass not implemented\n"));
    }
    virtual void DdaggerWithMass(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DdaggerWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DdaggerWithMass not implemented\n"));
    }
    virtual void DDWithMass(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDWithMass not implemented\n"));
    }
    virtual void DDdaggerWithMass(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0))
    {
        if (SingleField())
        {
            const CFieldGauge* pGauge = GetDefaultGauge(gaugeNum, gaugeFields);
            DDdaggerWithMassS(pGauge, fMass, eCoeffType, fCoeffReal, fCoeffImg);
            return;
        }
        appCrucial(_T("DDdaggerWithMass not implemented\n"));
    }

    virtual UBOOL InverseD(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson);
    virtual UBOOL InverseDdagger(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson);
    virtual UBOOL InverseDD(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson);
    virtual UBOOL InverseDDdagger(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson);
    virtual void InitialAsSource(const SFermionBosonSource& sourceData) = 0;

    virtual TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const = 0;
    virtual UBOOL RationalApproximation(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const class CRatinalApproximation* pRational);

    UBOOL IsFermionField() const override { return TRUE; }
    UINT GetSiteCount() const { return m_uiSiteCount; }

protected:

#pragma region single field case

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

    virtual void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
    {
        appCrucial(_T("DerivateDOperator not implemented\n"));
    }

#pragma endregion

    UINT m_uiLinkeCount;
    UINT m_uiSiteCount;
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
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson,
    const class CFieldFermionWilsonSquareSU3* pFermion);

extern CLGAPI void ExportDiagnalStaggeredSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson,
    const class CFieldFermionKSSU3* pFermion);


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================