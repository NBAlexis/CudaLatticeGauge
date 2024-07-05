//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3D.h
// 
// DESCRIPTION:
//
// To make sure the usual case is faster, write specified code for Dirichlet
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================
#include "CFieldFermionWilsonSquareSU3.h"

#ifndef _CFIELDFERMIONWILSONSQUARESU3D_H_
#define _CFIELDFERMIONWILSONSQUARESU3D_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3D)

class CLGAPI CFieldFermionWilsonSquareSU3D : public CFieldFermionWilsonSquareSU3
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3D)

public:

    CFieldFermionWilsonSquareSU3D() : CFieldFermionWilsonSquareSU3() {}

    void FixBoundary() override;

    /**
    * For test
    */
    void PrepareForHMCOnlyRandomize();
    void PrepareForHMCNotRandomize(const CFieldGauge* pGauge);

    CCString GetInfos(const CCString& tab) const override;

protected:

    void PrepareForHMCS(const CFieldGauge* pGauge) override;
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;

    

    //No need to rewrite them, once D operator and Random Gaussian make sure boundary
    //All fields come from other operators (D, D+, DD+, etc), or Linear Algebra will keep the boundary
    //However, one can rewrite those functions for optimization (not to calculate the sites on surface).

    //virtual void AxpyPlus(const CField* x);
    //virtual void AxpyMinus(const CField* x);
    //virtual void Axpy(Real a, const CField* x);
    //virtual void Axpy(const CLGComplex& a, const CField* x);
    //virtual void ScalarMultply(const CLGComplex& a);
    //virtual void ScalarMultply(Real a);
    //virtual CLGComplex Dot(const CField* other) const;
};



__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================