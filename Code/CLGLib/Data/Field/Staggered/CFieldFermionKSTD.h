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

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0_D(
            this->m_pDeviceData, 
            this->m_byFieldId, 
            (deviceGauge*)pForce, 
            (const deviceGauge*)pGaugeBuffer,
            byGaugeFieldId, 
            this->m_pRationalFieldPointers, 
            this->m_pMDNumerator, 
            this->m_rMD.m_uiDegree);
    }

    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_D(
            this->m_bEachSiteEta, (deviceVector*)pTargetBuffer, 
            (const deviceVector*)pBuffer, 
            (const deviceGauge*)pGaugeBuffer, 
            this->m_byFieldId, 
            byGaugeFieldId, 
            f2am, 
            bDagger, 
            eOCT, 
            fRealCoeff, 
            cCmpCoeff);
    }

public:

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        this->InitialField(EFIT_RandomGaussian);
        FixBoundary();
        this->D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
        FixBoundary();

        if (NULL != appGetFermionSolver(this->m_byFieldId) && !appGetFermionSolver(this->m_byFieldId)->IsAbsoluteAccuracy())
        {
            this->m_fLength = this->Dot(this).x;
        }

    }

    void FixBoundary() override
    {
        CCommonKernelSite<deviceVector>::FixBoundary(this->m_pDeviceData, this->m_byFieldId);
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GetInfos(tab);

        SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
        sRet = sRet + tab + appToString(boundary) + _T("\n");
        return sRet;
    }

};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1D)
class CLGAPI CFieldFermionKSU1D : public CFieldFermionKSTD<CLGComplex, CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSU1D)
public:
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredU1; }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSTD_H_

//=============================================================================
// END OF FILE
//=============================================================================