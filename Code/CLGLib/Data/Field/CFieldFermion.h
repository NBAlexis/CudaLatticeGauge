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

__DEFINE_ENUM(EFermionSource,
    EFS_Point,
    EFS_Wall,
    EFS_MomentumWall,
)

//this is a host only structure
struct SFermionSource
{
    EFermionSource m_eSourceType;
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

    virtual void PrepareForHMC(const CFieldGauge* pGauge) = 0;

    /**
    * Calculate force can fail due to solver
    */
    virtual UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const = 0;

    UBOOL ApplyOperator(EFieldOperator op, const CField * otherfield, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* pOtherParameters = NULL) override
    {
        switch (op)
        {
        case EFO_F_D:
            D(otherfield, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_Ddagger:
            Ddagger(otherfield, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DD:
            DD(otherfield, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_DDdagger:
            DDdagger(otherfield, eCoeffType, fCoeffReal, fCoeffImg);
            return TRUE;
        case EFO_F_InverseD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseD(otherfield);
        case EFO_F_InverseDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDdagger(otherfield);
        case EFO_F_InverseDDdagger:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDDdagger(otherfield);
        case EFO_F_InverseDD:
            if (EOCT_None != eCoeffType)
            {
                appCrucial(_T("ApplyOperator, the operator %s with coefficient is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                return FALSE;
            }
            return InverseDDdagger(otherfield);
        case EFO_F_RationalD:
            {
                class CRatinalApproximation* pRA = (class CRatinalApproximation*)(pOtherParameters);
                if (NULL == pRA)
                {
                    appCrucial(_T("ApplyOperator, the operator %s with CRatinalApproximation null.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
                    return FALSE;
                }
                return RationalApproximation(EFO_F_D, otherfield, pRA);
            }
        case EFO_F_D_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DWithMass(otherfield, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DWithMass(otherfield, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_Ddagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DdaggerWithMass(otherfield, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DdaggerWithMass(otherfield, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DD_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDWithMass(otherfield, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDWithMass(otherfield, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        case EFO_F_DDdagger_WithMass:
            {
                Real* fMass = (Real*)pOtherParameters;
                if (NULL == fMass)
                {
                    DDdaggerWithMass(otherfield, CCommonData::m_fShiftedMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
                else
                {
                    DDdaggerWithMass(otherfield, *fMass, eCoeffType, fCoeffReal, fCoeffImg);
                }
            }
            return TRUE;
        default:
            appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
            return FALSE;
        }
    }

    virtual void ApplyGamma(EGammaMatrix eGamma) = 0;

    virtual void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DD(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DdaggerWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DDWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;
    virtual void DDdaggerWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) = 0;

    virtual UBOOL InverseD(const CField* pGauge) = 0;
    virtual UBOOL InverseDdagger(const CField* pGauge) = 0;
    virtual UBOOL InverseDD(const CField* pGauge) = 0;
    virtual UBOOL InverseDDdagger(const CField* pGauge) = 0;
    virtual void InitialAsSource(const SFermionSource& sourceData) = 0;
    virtual TArray<CFieldFermion*> GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const = 0;

    virtual UBOOL RationalApproximation(EFieldOperator op, const CField* pGauge, const class CRatinalApproximation* pRational);

#pragma region real operators

    /**
    * It seems no need to create a "Foprt" class like Bridge++
    * TODO: this only support single gauge field...
    */
    virtual void DOperator(void* pTargetBuffer, const void * pBuffer, const void * pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const = 0;
    virtual void DerivateDOperator(void* pForce, const void * pDphi, const void * pDDphi, const void * pGaugeBuffer) const = 0;

#pragma endregion

    UBOOL IsGaugeField() const override { return FALSE; }
    UBOOL IsFermionField() const override { return TRUE; }

    /**
     * For even odd preconditioner
     */
    virtual void WriteEvenSites(const CFieldFermion*, const CFieldGauge*, UBOOL) { appCrucial(_T("Not implemented.\n")); }
    virtual void WriteBackEvenSites(CFieldFermion*, const CFieldGauge*, UBOOL) const { appCrucial(_T("Not implemented.\n")); }

    UINT GetSiteCount() const { return m_uiSiteCount; }

protected:

    virtual UBOOL InverseD_eo(const CField*)
    {
        appCrucial(_T("Not implemented.\n"));
        return FALSE;
    }

    virtual UBOOL InverseDdagger_eo(const CField*)
    {
        appCrucial(_T("Not implemented.\n"));
        return FALSE;
    }

    virtual UBOOL InverseDDdagger_eo(const CField*)
    {
        appCrucial(_T("Not implemented.\n"));
        return FALSE;
    }

    UINT m_uiLinkeCount;
    UINT m_uiSiteCount;
    SBYTE m_byEvenFieldId;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMION_H_

//=============================================================================
// END OF FILE
//=============================================================================