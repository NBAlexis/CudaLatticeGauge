//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Gamma.h
// 
// DESCRIPTION:
// If m_bExpGamma is turned on, it is applied as
// gamma4.Exp(gamma)
//
// REVISION:
//  [08/27/2023 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3GAMMA_H_
#define _CFIELDFERMIONWILSONSQUARESU3GAMMA_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3Gamma)

class CLGAPI CFieldFermionWilsonSquareSU3Gamma : public CFieldFermionWilsonSquareSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3Gamma)

public:

    CFieldFermionWilsonSquareSU3Gamma();
    ~CFieldFermionWilsonSquareSU3Gamma();

    void InitialOtherParameters(CParameters& params) override;
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;
    CCString GetInfos(const CCString& tab) const override;
    
    /**
    * If m_bExpGamma is turned on, it is applied as
    * gamma4.Exp(gamma)
    */
    UBOOL m_bExpGamma;

    Real m_fCoeffGamma1;
    Real m_fCoeffGamma2;
    Real m_fCoeffGamma3;
    Real m_fCoeffGamma4;
    Real m_fCoeffGamma5;
    Real m_fCoeffGamma51;
    Real m_fCoeffGamma52;
    Real m_fCoeffGamma53;
    Real m_fCoeffGamma54;
    Real m_fCoeffSigma12;
    Real m_fCoeffSigma13;
    Real m_fCoeffSigma14;
    Real m_fCoeffSigma23;
    Real m_fCoeffSigma24;
    Real m_fCoeffSigma34;

    /**
    * apply gamma4.Gamma using gamma4.exp(iGamma)
    * This is for simulation, so the gamma matrix is multiplied by an "I", if it is not gamma5-Hermitian
    */
    static void appApplyGammaExp(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
        const SIndex* pGaugeMove,
        const SIndex* pFermionMove,
        EGammaMatrix eGamma,
        UBOOL bDagger,
        UBOOL bExp,
        Real fGammaCoeff,
        Real fKappa,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * 
     */
    static void GammaForceExp(
        void* pForce,
        const void* pGaugeBuffer,
        const void* InverseD,
        const void* InverseDDagger,
        const SIndex* pFermionMove,
        UBOOL bExp,
        Real fGammaCoeff,
        Real fKappa,
        EGammaMatrix eGamma,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3GAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================