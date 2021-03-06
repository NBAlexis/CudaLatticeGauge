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
        , m_fa2Ez(F(0.0))
        , m_fa2Bz(F(0.0))
    {
        
    }

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fa2Ez;
    Real m_fa2Bz;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3EM_H_

//=============================================================================
// END OF FILE
//=============================================================================