//=============================================================================
// FILENAME : CGaugeSmearingASQTAD.h
// 
// DESCRIPTION:
// 
// Upto fat 7, it is:
// 
// 1/8 + 6/16 + (6*4)/64 + (6*4*2)/384
//
// REVISION:
//  [mm/dd/yy]
//  [11/14/2024 nbale]
//=============================================================================
#ifndef _CGAUGESMEARINGASQTAD_H_
#define _CGAUGESMEARINGASQTAD_H_

__BEGIN_NAMESPACE

//How to calculate (U^+U)^(-1/2)
__DEFINE_ENUM(EProjectOneOverTwo, 
    EPOOT_RationalApproximation)

template<typename gaugetype, INT matrixN>
class __DLL_EXPORT CGaugeSmearingASQTAD : public CGaugeSmearing
{
public:
    CGaugeSmearingASQTAD()
        : CGaugeSmearing()
        , m_fOriginal(F(0.125))
        , m_fFat3(F(0.0625))
        , m_fFat5(F(0.015625))
        , m_fFat7(F(0.00260417))
        , m_fLepage(F(0.0))
        , m_pGaugeNotProjected(NULL)
        , m_pP3_1(NULL)
        , m_pP3_2(NULL)
        , m_pP3_3(NULL)
        , m_pP5_1_1(NULL)
        , m_pP5_2_1(NULL)
        , m_pP5_3_1(NULL)
        , m_pP5_1_2(NULL)
        , m_pP5_2_2(NULL)
        , m_pP5_3_2(NULL)
        , m_pProjForcePtr(NULL)
        , m_pQ(NULL)
        , m_pQ2(NULL)
        , m_pInverseSqrtQ(NULL)
        , m_pCaylayHamiltonCache(NULL)
        , m_bProj(TRUE)
        , m_bProjDet(FALSE)
        , m_bUseCaylayHamilton(TRUE)
        , m_uiRationalApproximationOrder(0)
        , m_pDeviceRationalApproximation(NULL)
        , m_pDeviceDet(NULL)
        , m_pDeviceDetPhase(NULL)
    {

    }

    ~CGaugeSmearingASQTAD();

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeSmearing(class CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject = TRUE) override;
    CCString GetInfos(const CCString& sTab) const override;

    /**
    * since the gauge is changed by gauge smearing, however, to calculate the force one need original gauge, we store it in "GaugeSmearing"
    * therefore we don't need to send in the pGauge
    */
    void DerivateOnU(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const override;

    virtual EFieldType GetFieldType() const = 0;
    virtual UBOOL UseCaylayHamilton() const { return FALSE; }

    Real m_fOriginal;
    Real m_fFat3;
    Real m_fFat5;
    Real m_fFat7;
    Real m_fLepage;

    CFieldGauge* m_pGaugeNotProjected;
    CFieldGauge* m_pP3_1;
    CFieldGauge* m_pP3_2;
    CFieldGauge* m_pP3_3;

    CFieldGauge* m_pP5_1_1;
    CFieldGauge* m_pP5_2_1;
    CFieldGauge* m_pP5_3_1;
    CFieldGauge* m_pP5_1_2;
    CFieldGauge* m_pP5_2_2;
    CFieldGauge* m_pP5_3_2;

    TArray<CFieldGauge*> m_pProjForce;
    gaugetype** m_pProjForcePtr;
    CFieldGauge* m_pQ;
    CFieldGauge* m_pQ2;
    CFieldGauge* m_pInverseSqrtQ;
    DOUBLE* m_pCaylayHamiltonCache;

    //UR = U(U^+U)^{-1/2}
    UBOOL m_bProj;

    //UR' = UR/det[U]^{1/3}
    UBOOL m_bProjDet;

    UBOOL m_bUseCaylayHamilton;

    //========= Used for projection =============
    //CRatinalApproximation m_ra; we directly use a device buffer
    UINT m_uiRationalApproximationOrder;
    Real* m_pDeviceRationalApproximation;
    CLGComplex* m_pDeviceDet;
    DOUBLE* m_pDeviceDetPhase;

    /**
    * When NULL != pPreserve, pPreserve is used to store pGauge, and pGaugeOrignal is no need
    * When NULL == pPreserve, pGaugeOrignal is used as orignal one (NOTE that, this is the orignal one ONLY at the first level)
    */
    static void Fat357Lepage(class CFieldGauge* pGauge, class CFieldGauge* preserve, class CFieldGauge* pp3_1, class CFieldGauge* pp3_2, class CFieldGauge* pp3_3,
        class CFieldGauge* pp5_1_1, class CFieldGauge* pp5_1_2, class CFieldGauge* pp5_2_1, class CFieldGauge* pp5_2_2, class CFieldGauge* pp5_3_1, class CFieldGauge* pp5_3_2,
        Real fOrignal, Real fFat3, Real fFat5, Real fFat7, Real fLepage, const class CFieldGauge* pGaugeOrignal);

    static void CalculateFat3Only(
        const CFieldGauge* pGauge,
        CFieldGauge* p3_1, CFieldGauge* p3_2, CFieldGauge* p3_3);

    static void CalculateFat5Only(
        const CFieldGauge* pGauge,
        const CFieldGauge* p3_1, const CFieldGauge* p3_2, const CFieldGauge* p3_3,
        CFieldGauge* p5_1_1, CFieldGauge* p5_1_2, CFieldGauge* p5_2_1,
        CFieldGauge* p5_2_2, CFieldGauge* p5_3_1, CFieldGauge* p5_3_2);

    static void SmearingForce(const class CFieldGauge* pOrignalGauge, class CFieldGauge* pf0,
        const class CFieldGauge* p3_1, const class CFieldGauge* p3_2, const class CFieldGauge* p3_3,
        const class CFieldGauge* p5_1_1, const class CFieldGauge* p5_1_2, const class CFieldGauge* p5_2_1, const class CFieldGauge* p5_2_2, const class CFieldGauge* p5_3_1, const class CFieldGauge* p5_3_2,
        Real fOrignal, Real fFat3, Real fFat5, Real fFat7, Real fLepage);

    static void ProjectCaylayHamilton(class CFieldGauge* pGauge, class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pQ, class CFieldGauge* pQ2, class CFieldGauge* pInverseSqrtQ, DOUBLE* constants, UBOOL bToSU3, DOUBLE* detphase);

    static void ProjectCaylayHamiltonForce(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pf0,
        const class CFieldGauge* pQ, const class CFieldGauge* pQ2, const class CFieldGauge* pInverseSqrtQ, const DOUBLE* constants, UBOOL bToSU3, const DOUBLE* detphase);

    static void ProjectRationalApproximation(class CFieldGauge* pGauge, class CFieldGauge* pGaugeBeforeProj, gaugetype** forcePointers,
        const Real* rationalCoeffs, UINT uiOrder, UBOOL bToSU3, CLGComplex* det);

    static void ProjectRationalApproximationForce(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pf0, 
        const gaugetype* const* forcePointers, const Real* rationalCoeffs, UINT uiOrder, UBOOL bToSU3, const CLGComplex* det);
};

__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingASQTADSU3)
class CLGAPI CGaugeSmearingASQTADSU3 : public CGaugeSmearingASQTAD<deviceSU3, 3>
{
    __CLGDECLARE_CLASS(CGaugeSmearingASQTADSU3)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
    UBOOL UseCaylayHamilton() const override  { return this->m_bUseCaylayHamilton; }
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGASQTAD_H_

//=============================================================================
// END OF FILE
//=============================================================================