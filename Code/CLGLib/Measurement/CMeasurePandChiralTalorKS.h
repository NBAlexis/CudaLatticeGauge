//=============================================================================
// FILENAME : CMeasurePandChiralTalorKS.h
// 
// DESCRIPTION:
// This is measurement for Talor expansion of rotation frame
// NOTE : The fermion part is not fully implemented!...
//
// REVISION:
//  [16/11/2022 nbale]
//=============================================================================

#ifndef _CMEASUREPANDCHIRALTALORKS_H_
#define _CMEASUREPANDCHIRALTALORKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasurePandChiralTalorKS)

__DEFINE_ENUM(ECPCTTraceTypeKS,
    ECPCTTTKS_D,
    ECPCTTTKS_MD,
    ECPCTTTKS_DMD,
    ECPCTTTKS_MDMD,
    ECPCTTTKS_DMDMD,

    ECPCTTTKS_Max,
    );

class CLGAPI CMeasurePandChiralTalorKS : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasurePandChiralTalorKS)
public:
    CMeasurePandChiralTalorKS()
        : CMeasureStochastic()
    {
    }

    void OnConfigurationAcceptedZ4SingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }

protected:
    
    static void ApplyM(CFieldFermionKSSU3* pTarget, const CFieldFermionKSSU3* pSource, const CFieldGaugeSU3* pGauge);

public:

#if _CLG_DOUBLEFLOAT
    CLGComplex m_cTmpSum[ECPCTTT_Max];
    TArray<CLGComplex> m_lstTraceRes[ECPCTTT_Max];
    TArray<CLGComplex> m_lstPolyakov;
    TArray<Real> m_lstPolyakovSOmega;
    TArray<Real> m_lstPolyakovSOmegaSq;
#else
    cuDoubleComplex m_cTmpSum[ECPCTTTKS_Max];
    TArray<cuDoubleComplex> m_lstTraceRes[ECPCTTTKS_Max];
    TArray<cuDoubleComplex> m_lstPolyakov;
    TArray<DOUBLE> m_lstPolyakovSOmega;
    TArray<DOUBLE> m_lstPolyakovSOmegaSq;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPANDCHIRALTALORKS_H_

//=============================================================================
// END OF FILE
//=============================================================================