//=============================================================================
// FILENAME : CFieldFermionWilsonSU3EvenOdd.h
// 
// DESCRIPTION:
//
// 
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSU3DEVENODD_H_
#define _CFIELDFERMIONWILSONSU3DEVENODD_H_

__BEGIN_NAMESPACE

#if FURTURE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSU3EvenOdd)

class CLGAPI CFieldFermionWilsonSU3EvenOdd : public CFieldFermionWilsonSquareSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSU3EvenOdd)

public:

    CFieldFermionWilsonSU3EvenOdd() : CFieldFermionWilsonSquareSU3D() {}


    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;
    CCString GetInfos(const CCString& tab) const override;

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

    //Dee = Doo = 1
    virtual void Doe() {}
    virtual void Deo() {}
    virtual void Deooe() {}
    virtual void InverseDoe() {}
    virtual void InverseDeo() {}
};

#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSU3DEVENODD_H_

//=============================================================================
// END OF FILE
//=============================================================================