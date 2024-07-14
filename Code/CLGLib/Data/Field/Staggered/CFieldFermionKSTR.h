//=============================================================================
// FILENAME : CFieldFermionKSTR.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [07/14/2020 nbale]
//=============================================================================
#include "CFieldFermionKST.h"

#ifndef _CFIELDFERMIONKSTR_H_
#define _CFIELDFERMIONKSTR_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKSTR : public CFieldFermionKST<deviceVector, deviceGauge, vectorN>
{
public:

    CFieldFermionKSTR() : CFieldFermionKST<deviceVector, deviceGauge, vectorN>()
        , m_bRealRotation(FALSE)
        , m_fOmega(F(0.0))
    {

    }

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;
    void CopyTo(CField* u) const override;

    UBOOL m_bRealRotation;
    Real m_fOmega;
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1R)
class CLGAPI CFieldFermionKSU1R : public CFieldFermionKSTR<CLGComplex, CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSU1R)
public:
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredU1; }
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSTR_H_

//=============================================================================
// END OF FILE
//=============================================================================