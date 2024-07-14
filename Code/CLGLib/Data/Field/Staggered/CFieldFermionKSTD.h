//=============================================================================
// FILENAME : CFieldFermionKSTD.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [07/14/2024 nbale]
//=============================================================================
#include "CFieldFermionKST.h"

#ifndef _CFIELDFERMIONKSTD_H_
#define _CFIELDFERMIONKSTD_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKSTD : public CFieldFermionKST<deviceVector, deviceGauge, vectorN>
{
protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override;
    void FixBoundary() override;
    

    CCString GetInfos(const CCString& tab) const override;

};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSTD_H_

//=============================================================================
// END OF FILE
//=============================================================================