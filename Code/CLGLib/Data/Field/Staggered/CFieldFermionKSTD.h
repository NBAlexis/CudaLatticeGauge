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
//  [mm/dd/yy]
//  [07/14/2024 nbale]
//=============================================================================
#include "CFieldFermionKST.h"

#ifndef _CFIELDFERMIONKSTD_H_
#define _CFIELDFERMIONKSTD_H_

__BEGIN_NAMESPACE

template <class CFieldKS>
class __DLL_EXPORT CFieldFermionKSTD : public CFieldKS
{
protected:
    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        CFieldKS::_FermionKernel::DerivateD0_D(
            this->m_pDeviceData,
            this->m_byFieldId,
            (typename CFieldKS::_Gauge*)pForce,
            (const typename CFieldKS::_Gauge*)pGaugeBuffer,
            byGaugeFieldId,
            CRationalFieldPointer::GetInstance()->GetRationPoint<typename CFieldKS::_Vector>(this->m_byRationFieldPointerBufferLength),
            GRASet.m_pRASet[this->m_iMDIndex]->m_pDeviceData,
            GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree);
    }

    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldKS::_FermionKernel::DOperatorKS_D(
            this->m_bEachSiteEta, 
            (typename CFieldKS::_Vector*)pTargetBuffer,
            (const typename CFieldKS::_Vector*)pBuffer,
            (const typename CFieldKS::_Gauge*)pGaugeBuffer,
            this->m_byFieldId,
            byGaugeFieldId,
            f2am,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldKS::_FermionKernel::DOperatorKSOnEvenOrOdd_D(
            (typename CFieldKS::_Vector*)pTargetBuffer,
            (const typename CFieldKS::_Gauge*)pGaugeBuffer,
            this->m_byFieldId, 
            byGaugeFieldId, 
            bEven, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

public:

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override
    {
        if (this->m_bEvenPseudofermion)
        {
            appParanoiac(_T("PrepareForHMC with Even field\n"));
            CFieldKS::_SiteKernel::InitialBufferEvenOdd(this->m_pDeviceData, this->m_byFieldId, TRUE, TRUE, EFIT_RandomGaussian);
            FixBoundary(EFB_Momentum);
            //this->DebugPrintMe();
            //Note that, using D0_MC, A=D0D0 + m^2, A^{-1/2} should be expand with m^2 considered.
            this->D_MC(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
            FixBoundary(EFB_Momentum);
        }
        else
        {
            CFieldKS::_SiteKernel::InitialBuffer(this->m_pDeviceData, this->m_byFieldId, EFIT_RandomGaussian);
            FixBoundary(EFB_Momentum);
            this->D_MC(gaugeNum, bosonNum, 0, gaugeFields, pBoson, NULL);
            FixBoundary(EFB_Momentum);
        }

        if (NULL != appGetFermionSolver(this->m_byFieldId) && !appGetFermionSolver(this->m_byFieldId)->IsAbsoluteAccuracy())
        {
            this->m_fLength = this->Dot(this).x;
        }
    }

    void FixBoundary(EFixBoundary eType) override
    {
        CFieldKS::_SiteKernel::FixBoundary(this->m_pDeviceData, this->m_byFieldId);
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldKS::GetInfos(tab);
        SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
        sRet = sRet + tab + appToString(boundary) + _T("\n");
        return sRet;
    }

    UBOOL IsDirichlet() const override
    {
        return TRUE;
    }
};


#define _DEFINE_KS_Dirichlet(N) \
__CLG_REGISTER_HELPER_HEADER(CFieldFermionKS##N##D) \
class CLGAPI CFieldFermionKS##N##D : public CFieldFermionKSTD<CFieldFermionKS##N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKS##N##D) \
};

_DEFINE_KS_Dirichlet(U1)
_DEFINE_KS_Dirichlet(SU3)

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSTD_H_

//=============================================================================
// END OF FILE
//=============================================================================