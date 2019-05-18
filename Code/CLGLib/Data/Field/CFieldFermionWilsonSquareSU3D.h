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

#ifndef _CFIELDFERMIONWILSONSQUARESU3D_H_
#define _CFIELDFERMIONWILSONSQUARESU3D_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3D)

class CLGAPI CFieldFermionWilsonSquareSU3D : public CFieldFermionWilsonSquareSU3
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3D)

public:

    CFieldFermionWilsonSquareSU3D() : CFieldFermionWilsonSquareSU3() {}

    virtual void FixBoundary();
    virtual void PrepareForHMC(const CFieldGauge* pGauge);
    virtual void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const;
    virtual void DerivateDOperator(void* pForce, void* pCacheForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const;
    virtual CCString GetInfos(const CCString &tab) const;

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


#pragma region device functions

static __device__ __inline__ deviceWilsonVectorSU3 _deviceGetFermionBCWilsonSU3(
    const deviceWilsonVectorSU3* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    return idx.IsDirichlet() ?
        deviceWilsonVectorSU3::makeZeroWilsonVectorSU3()
        : pBuffer[idx.m_uiSiteIndex];
}

static __device__ __inline__ deviceWilsonVectorSU3 _deviceGetFermionBCWilsonSU3T(
    const deviceWilsonVectorSU3* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    return pBuffer[idx.m_uiSiteIndex];
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================