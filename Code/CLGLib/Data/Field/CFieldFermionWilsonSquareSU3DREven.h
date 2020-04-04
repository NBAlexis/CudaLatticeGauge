//=============================================================================
// FILENAME : CFieldFermionWilsonSU3DREven.h
// 
// DESCRIPTION:
//   
//   Do not directly use it, copy a CFieldFermionWilsonSquareSU3 or CFieldFermionWilsonSquareSU3D instead
//   Assume Nx * Ny is even, for convinient to decompse threads
//
//   Note that, for DR, in fact, we only use naive discretization and exponential chemical potential.
//   In that case, Doo = Dee = 1, which is much simpler, so we only implement this situation.
//
// REVISION:
//  [03/31/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSU3DREVEN_H_
#define _CFIELDFERMIONWILSONSU3DREVEN_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSU3DREven)

class CLGAPI CFieldFermionWilsonSU3DREven : public CFieldFermionWilsonSU3DEven
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSU3DREven)

public:

    CFieldFermionWilsonSU3DREven() : CFieldFermionWilsonSU3DEven() {}

    //1 - kappa^2 Doo^(-1) Deo Dee^(-1)Doe
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    CCString GetInfos(const CCString& tab) const override;

    void WriteEvenSites(const CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) override;
    void WriteBackEvenSites(CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) const override;
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSU3DREVEN_H_

//=============================================================================
// END OF FILE
//=============================================================================