//=============================================================================
// FILENAME : CGaugeSmearingStout.h
// 
// DESCRIPTION:
// 
// The implementation use exp(iQ) of hep-lat/0311018, therefore supports only SU3 gauge field!
// I'm not sure whether to support other gauge in the future, the problem might be in the calculation of derivative of exp(iQ)
// If Q is small enough, we can calculate the derivative of exp(iQ) by Taylor expansion?
// 
// REVISION:
//  [mm/dd/yy]
//  [07/25/2025 nbale]
//=============================================================================
#ifndef _CGAUGESMEARINGSTOUT_H_
#define _CGAUGESMEARINGSTOUT_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EStoutLinkCache,
    ESLC_Full,
    ESLC_Median,
    ESLC_Small,
    );

__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingStoutSU3)
class CLGAPI CGaugeSmearingStoutSU3 : public CGaugeSmearing
{
    __CLGDECLARE_CLASS(CGaugeSmearingStoutSU3)
public:
    CGaugeSmearingStoutSU3();
    ~CGaugeSmearingStoutSU3();

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;

    void GaugeSmearing(class CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject = TRUE) override
    {
        if (ESLC_Median == m_eCache)
        {
            GaugeSmearingMedian(pGauge, pOrignal, pStaple, bProject);
        }
        else
        {
            GaugeSmearingFull(pGauge, pOrignal, pStaple, bProject);
        }
    }

    CCString GetInfos(const CCString& sTab) const override;

    /**
    * since the gauge is changed by gauge smearing, however, to calculate the force one need original gauge, we store it in "GaugeSmearing"
    * therefore we don't need to send in the pGauge
    */
    void DerivateOnU(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const override
    {
        if (ESLC_Median == m_eCache)
        {
            DerivateOnUMedian(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0);
        }
        else
        {
            DerivateOnUFull(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0);
        }
    }

    Real m_fRhojk;
    Real m_fRho4mu;

    //CFieldGauge* m_pOriginalGauge;

    //intermediate results can be used
    //(1) of hep-lat/0311018
    CFieldGauge* m_pC;
    //(2) of hep-lat/0311018
    CFieldGauge* m_pQ;
    CFieldGauge* m_pQ2;
    CFieldGauge* m_pExpIQ;
    CFieldGauge* m_pLambda;

    //5 complex parameters
    cuDoubleComplex* m_f012_exp2u_expm1u;

    //8 parameters
    DOUBLE* m_u_w2_cosw_xi0_denorm;
    //0: normal, 1: minus, 2: zero Q, not smearing
    BYTE* m_byMinus;

    EStoutLinkCache m_eCache;

protected:

    void GaugeSmearingFull(class CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject);
    void DerivateOnUFull(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const;

    void GaugeSmearingMedian(class CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject);
    void DerivateOnUMedian(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGSTOUT_H_

//=============================================================================
// END OF FILE
//=============================================================================