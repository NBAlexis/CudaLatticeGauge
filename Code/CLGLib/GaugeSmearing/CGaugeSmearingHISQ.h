//=============================================================================
// FILENAME : CGaugeSmearingHISQ.h
// 
// DESCRIPTION:
// 
// Upto fat 7, it is:
// 
// 1/8 + 6/16 + (6*4)/64 + (6*4*2)/384
// 
// Note: level1 has no lepage
// 
// level two contains Lepage and Naik (Naik is implemented in Fermion)
// level2:
// one-link: 1+(epsi/8) see also 1004.0342
// Lepage: -1/8
// Naik:     -(1+epsi)/24
// see: 0710.0737
// default use epsi=0
// for epsi, see:
// hep-lat/0610092
// 
// In QUDA and MILC, the coefficients of fat3, fat7 are negative, this is because staggered phase (\eta_{\mu}) is absorbed into gauge link
// 
// REVISION:
//  [mm/dd/yy]
//  [12/03/2024 nbale]
//=============================================================================
#include "CGaugeSmearingASQTAD.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"

#ifndef _CGAUGESMEARINGHISQ_H_
#define _CGAUGESMEARINGHISQ_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EHISQLinkCache,
    EHLC_Full,
    EHLC_Median,
    EHLC_None,
    );

template<typename gaugetype, INT matrixN>
class __DLL_EXPORT CGaugeSmearingHISQ : public CGaugeSmearing
{
public:
    CGaugeSmearingHISQ()
        : CGaugeSmearing()
        , m_fOriginalL1(F(0.125))
        , m_fFat3L1(F(0.0625))
        , m_fFat5L1(F(0.015625))
        , m_fFat7L1(F(0.00260416666666667))
        , m_fLepageL1(F(0.0))

        , m_fOriginalL2(F(1.0))
        , m_fFat3L2(F(0.0625))
        , m_fFat5L2(F(0.015625))
        , m_fFat7L2(F(0.00260416666666667))
        , m_fLepageL2(F(-0.125))

        , m_eCache(EHLC_Full)
        , m_pGaugeNotProjected(NULL)
        , m_pEffectiveGaugeL1(NULL)
        , m_pP3_1_L1(NULL)
        , m_pP3_2_L1(NULL)
        , m_pP3_3_L1(NULL)
        , m_pP5_1_1_L1(NULL)
        , m_pP5_2_1_L1(NULL)
        , m_pP5_3_1_L1(NULL)
        , m_pP5_1_2_L1(NULL)
        , m_pP5_2_2_L1(NULL)
        , m_pP5_3_2_L1(NULL)
        , m_pP3_1_L2(NULL)
        , m_pP3_2_L2(NULL)
        , m_pP3_3_L2(NULL)
        , m_pP5_1_1_L2(NULL)
        , m_pP5_2_1_L2(NULL)
        , m_pP5_3_1_L2(NULL)
        , m_pP5_1_2_L2(NULL)
        , m_pP5_2_2_L2(NULL)
        , m_pP5_3_2_L2(NULL)
        , m_pNaik(NULL)
        , m_pNaikLinkPhase(NULL)
        , m_iPhaseFieldId(-1)
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

    ~CGaugeSmearingHISQ();

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;

    void GaugeSmearing(class CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject = TRUE) override
    {
        if (EHLC_Median == m_eCache)
        {
            GaugeSmearingMedian(pGauge, pOrignalGauge, pStaple, bProject);
        }
        else if (EHLC_None == m_eCache)
        {
            GaugeSmearingNone(pGauge, pOrignalGauge, pStaple, bProject);
        }
        else
        {
            GaugeSmearingFull(pGauge, pOrignalGauge, pStaple, bProject);
        }
    }

    const CFieldGauge* GetEffectiveGaugeLevel1() const override
    {
        return m_pEffectiveGaugeL1;
    }

    CCString GetInfos(const CCString& sTab) const override;

    const CFieldGauge* GetNaikLink() const override
    {
        return m_pNaik;
    }

    const CFieldGauge* GetNaikLinkPhase() const
    {
        return dynamic_cast<const CFieldGauge*>(m_pNaikLinkPhase);
    }

    /**
    * since the gauge is changed by gauge smearing, however, to calculate the force one need original gauge, we store it in "GaugeSmearing"
    * therefore we don't need to send in the pGauge
    */
    void DerivateOnU(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const override
    {
        if (EHLC_Median == m_eCache)
        {
            DerivateOnUMedian(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0);
        }
        else if (EHLC_None == m_eCache)
        {
            DerivateOnUNone(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0);
        }
        else
        {
            DerivateOnUFull(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0);
        }
    }

    static void CalculateNaikLink(gaugetype* naik, const gaugetype* effectivel1, BYTE byFieldId);
    static void CalculateU1NaikLink(class CFieldGaugeU1Real* naik, const class CFieldGaugeU1Real* externalfield);

    virtual EFieldType GetFieldType() const = 0;
    virtual UBOOL UseCaylayHamilton() const { return FALSE; }

    void OnPhaseChanged()
    {
        m_pNaikLinkPhase->Return();
        m_pNaikLinkPhase = NULL;
    }

    Real m_fOriginalL1;
    Real m_fFat3L1;
    Real m_fFat5L1;
    Real m_fFat7L1;
    Real m_fLepageL1;

    Real m_fOriginalL2;
    Real m_fFat3L2;
    Real m_fFat5L2;
    Real m_fFat7L2;
    Real m_fLepageL2;

    EHISQLinkCache m_eCache;

    //CFieldGauge* m_pOriginalGauge;
    CFieldGauge* m_pGaugeNotProjected;
    CFieldGauge* m_pEffectiveGaugeL1;
    CFieldGauge* m_pP3_1_L1;
    CFieldGauge* m_pP3_2_L1;
    CFieldGauge* m_pP3_3_L1;

    CFieldGauge* m_pP5_1_1_L1;
    CFieldGauge* m_pP5_2_1_L1;
    CFieldGauge* m_pP5_3_1_L1;
    CFieldGauge* m_pP5_1_2_L1;
    CFieldGauge* m_pP5_2_2_L1;
    CFieldGauge* m_pP5_3_2_L1;

    CFieldGauge* m_pP3_1_L2;
    CFieldGauge* m_pP3_2_L2;
    CFieldGauge* m_pP3_3_L2;

    CFieldGauge* m_pP5_1_1_L2;
    CFieldGauge* m_pP5_2_1_L2;
    CFieldGauge* m_pP5_3_1_L2;
    CFieldGauge* m_pP5_1_2_L2;
    CFieldGauge* m_pP5_2_2_L2;
    CFieldGauge* m_pP5_3_2_L2;

    CFieldGauge* m_pNaik;

    //Phase cache supports only static field, use another gauge smearing otherwise
    class CFieldGaugeU1Real* m_pNaikLinkPhase;
    INT m_iPhaseFieldId;

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

protected:
    void GaugeSmearingFull(CFieldGauge* pGauge, const CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject);
    void GaugeSmearingMedian(CFieldGauge* pGauge, const CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject);
    void GaugeSmearingNone(CFieldGauge* pGauge, const CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject);

    void DerivateOnUFull(const CFieldGauge* pEffectiveGauge, const CFieldGauge* pOrignalGauge, const CFieldGauge* pNaikForce, CFieldGauge* pf0) const;
    void DerivateOnUMedian(const CFieldGauge* pEffectiveGauge, const CFieldGauge* pOrignalGauge, const CFieldGauge* pNaikForce, CFieldGauge* pf0) const;
    void DerivateOnUNone(const CFieldGauge* pEffectiveGauge, const CFieldGauge* pOrignalGauge, const CFieldGauge* pNaikForce, CFieldGauge* pf0) const;
};


__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingHISQSU3)
class CLGAPI CGaugeSmearingHISQSU3 : public CGaugeSmearingHISQ<deviceSU3, 3>
{
    __CLGDECLARE_CLASS(CGaugeSmearingHISQSU3)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
    UBOOL UseCaylayHamilton() const override { return this->m_bUseCaylayHamilton; }
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGHISQ_H_

//=============================================================================
// END OF FILE
//=============================================================================