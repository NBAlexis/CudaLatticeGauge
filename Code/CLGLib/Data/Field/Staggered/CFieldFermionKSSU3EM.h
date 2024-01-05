//=============================================================================
// FILENAME : CFieldFermionKSSU3EM.h
// 
// DESCRIPTION:
// This is the class for external electro-magnetic
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3EM_H_
#define _CFIELDFERMIONKSSU3EM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3EM)

class CLGAPI CFieldFermionKSSU3EM : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3EM)

public:

    CFieldFermionKSSU3EM() : CFieldFermionKSSU3()
        , m_fQ(F(0.0))
    {
        
    }

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;
    void ApplyGammaKS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

protected:

    Real m_fQ;

public:

    Real GetQ() const { return m_fQ; }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3EM_H_

//=============================================================================
// END OF FILE
//=============================================================================