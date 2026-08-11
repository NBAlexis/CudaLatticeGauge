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
//  [mm/dd/yy]
//  [07/14/2020 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKST.h"
#include "Data/Field/Staggered/CFieldFermionKSTD.h"
#include "Update/CStapleCache.h"

#ifndef _CFIELDFERMIONKSTR_H_
#define _CFIELDFERMIONKSTR_H_

__BEGIN_NAMESPACE

inline class CStapleCache* appGetStapleCache(BYTE byFieldId);

template<class CFieldKS>
class __DLL_EXPORT CFieldFermionKSTR : public CFieldKS
{
public:

    CFieldFermionKSTR() : CFieldKS()
        , m_bRealRotation(FALSE)
        , m_bShiftHalfCoord(FALSE)
        , m_fOmega(0.0)
        , m_byUseCachedGauge(0)
    {
    }

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        CFieldKS::DerivateD0(pForce, pGaugeBuffer, byGaugeFieldId);

        if (m_bRealRotation)
        {
            appCrucial(_T("DerivateD0 is not supported for real rotation!\n"));
        }

        CFieldKS::_FermionKernelR::DerivateD0_R(
            m_fOmega,
            m_bShiftHalfCoord,
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
        CFieldKS::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        if (m_bRealRotation)
        {
            CFieldKS::_FermionKernelR::DOperatorKS_R_RealRotation(
                m_fOmega,
                m_bShiftHalfCoord,
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
        else
        {
            if (1 == m_byUseCachedGauge && NULL != appGetStapleCache(byGaugeFieldId) && NULL != appGetStapleCache(byGaugeFieldId)->GetRotationBuffer())
            {
                CFieldKS::_FermionKernelR::DOperatorKS_R_ImaginaryRotation_Cached(
                    m_fOmega,
                    m_bShiftHalfCoord,
                    (typename CFieldKS::_Vector*)pTargetBuffer,
                    (const typename CFieldKS::_Vector*)pBuffer,
                    (const typename CFieldKS::_Gauge*)appGetStapleCache(byGaugeFieldId)->GetRotationBuffer(),
                    this->m_byFieldId,
                    byGaugeFieldId,
                    f2am,
                    bDagger,
                    eOCT,
                    fRealCoeff,
                    cCmpCoeff);
            }
            else
            {
                CFieldKS::_FermionKernelR::DOperatorKS_R_ImaginaryRotation(
                    m_fOmega,
                    m_bShiftHalfCoord,
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
        }
    }

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        //The rotation terms connect even and odd sites only, so even-odd works for torus and Dirichlet boundary
        if (m_bRealRotation)
        {
            appCrucial(_T("DOperatorKSOnEvenOrOdd is not supported for real rotation!\n"));
        }
        CFieldKS::DOperatorKSOnEvenOrOdd(pTargetBuffer, pGaugeBuffer, byGaugeFieldId, bEven, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        if (1 == m_byUseCachedGauge && NULL != appGetStapleCache(byGaugeFieldId) && NULL != appGetStapleCache(byGaugeFieldId)->GetRotationBuffer())
        {
            CFieldKS::_FermionKernelR::DOperatorKSOnEvenOrOdd_R_ImaginaryRotation_Cached(
                m_fOmega,
                m_bShiftHalfCoord,
                (typename CFieldKS::_Vector*)pTargetBuffer,
                (const typename CFieldKS::_Gauge*)appGetStapleCache(byGaugeFieldId)->GetRotationBuffer(),
                this->m_byFieldId,
                byGaugeFieldId,
                bEven,
                bDagger,
                eOCT,
                fRealCoeff,
                cCmpCoeff);
        }
        else
        {
            CFieldKS::_FermionKernelR::DOperatorKSOnEvenOrOdd_R_ImaginaryRotation(
                m_fOmega,
                m_bShiftHalfCoord,
                (typename CFieldKS::_Vector*)pTargetBuffer,
                (const typename CFieldKS::_Gauge*)pGaugeBuffer,
                this->m_byFieldId,
                byGaugeFieldId,
                bEven,
                bDagger,
                eOCT,
                fRealCoeff,
                cCmpCoeff);
        }
    }

    /**
    * Even-odd rotation force of one rational term
    * phi is the packed field with phi on even sites and D0 phi on odd sites
    * the force kernels need phi and D0 phi as two seperated fields,
    * so split the packed field into phi (odd zeroed) and D0 phi (even zeroed)
    */
    void CalculateForceEvenOddS_SingleTermOfRationalR(const CFieldGauge* pGauge, const CFieldFermionKS* phi, CFieldGauge* pForce, INT indexOfRationalTerm) const
    {
        if (m_bRealRotation)
        {
            appCrucial(_T("DerivateD0 is not supported for real rotation!\n"));
        }
        CField* pPhiOnly = appGetLattice()->GetPooledCopy(phi, _T(__FILE__), __LINE__);
        CField* pD0Only = appGetLattice()->GetPooledCopy(phi, _T(__FILE__), __LINE__);
        pPhiOnly->ZeroOnEvenOdd(FALSE);
        pD0Only->ZeroOnEvenOdd(TRUE);
        typename CFieldKS::_Vector* hostPointers[2] = {
            (typename CFieldKS::_Vector*)(pPhiOnly->GetData()),
            (typename CFieldKS::_Vector*)(pD0Only->GetData())
        };
        typename CFieldKS::_Vector** pPointers = CRationalFieldPointer::GetInstance()->GetRationPoint<typename CFieldKS::_Vector>(2);
        appSimpleCopyHD(pPointers, hostPointers, sizeof(typename CFieldKS::_Vector*) * 2);
        CFieldKS::_FermionKernelR::DerivateD0_ROnEvenOdd(
            m_fOmega,
            m_bShiftHalfCoord,
            this->m_pDeviceData,
            this->m_byFieldId,
            (typename CFieldKS::_Gauge*)(pForce->GetData()),
            (const typename CFieldKS::_Gauge*)(pGauge->GetData()),
            pGauge->m_byFieldId,
            pPointers,
            GRASet.m_pRASet[this->m_iMDIndex]->m_pDeviceData + indexOfRationalTerm,
            1);
        pPhiOnly->Return();
        pD0Only->Return();
    }

    void CalculateForceEvenOddS_SingleTermOfRational(const CFieldGauge* pGauge, const CFieldFermionKS* phi, CFieldGauge* pForce, Real fCoef, INT i) const override
    {
        CFieldKS::CalculateForceEvenOddS_SingleTermOfRational(pGauge, phi, pForce, fCoef, i);
        CalculateForceEvenOddS_SingleTermOfRationalR(pGauge, phi, pForce, i);
    }

public:

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldKS::InitialOtherParameters(params);

        // Keep the old EachSiteEta-implies-shift behavior when ShiftCoord is absent.
        INT iShiftCoord = this->m_bEachSiteEta ? 1 : 0;
        if (!params.FetchValueINT(_T("ShiftCoord"), iShiftCoord))
        {
            params.FetchValueINT(_T("ShiftCenter"), iShiftCoord);
        }
        m_bShiftHalfCoord = (0 != iShiftCoord);
        if (this->m_bEachSiteEta && !m_bShiftHalfCoord)
        {
            appCrucial(_T("CFieldFermionKSTR: EachSiteEta=1 requires ShiftCoord=1.\n"));
        }

        INT iReal = 0;
        if (params.FetchValueINT(_T("RealRotation"), iReal))
        {
            m_bRealRotation = (0 != iReal);
        }

        DOUBLE fValue = 0.1;
        if (params.FetchValueDOUBLE(_T("Omega"), fValue))
        {
            m_fOmega = fValue;
        }

        iReal = 0;
        if (params.FetchValueINT(_T("CachedGauge"), iReal))
        {
            m_byUseCachedGauge = static_cast<BYTE>(iReal);
        }
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldKS::GetInfos(tab);
        sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
        sRet = sRet + tab + _T("RealRotation : ") + appToString(m_bRealRotation) + _T("\n");
        sRet = sRet + tab + _T("ShiftCoord : ") + appToString(m_bShiftHalfCoord) + _T("\n");
        sRet = sRet + tab + _T("CachedGauge : ") + appToString(m_byUseCachedGauge) + _T("\n");
        return sRet;
    }

    void CopyParamTo(CField* U) const override
    {
        CFieldKS::CopyParamTo(U);
        CFieldFermionKSTR<CFieldKS>* pOther = dynamic_cast<CFieldFermionKSTR<CFieldKS>*>(U);
        if (NULL != pOther)
        {
            pOther->m_bRealRotation = m_bRealRotation;
            pOther->m_bShiftHalfCoord = m_bShiftHalfCoord;
            pOther->m_fOmega = m_fOmega;
            pOther->m_byUseCachedGauge = m_byUseCachedGauge;
        }
    }

    UBOOL m_bRealRotation;
    UBOOL m_bShiftHalfCoord;

    void SetFermionOmega(DOUBLE fOmega)
    {
        m_fOmega = fOmega;
        //this->UpdatePooledParamters();
    }

    DOUBLE GetOmega() const { return m_fOmega; }

protected:

    DOUBLE m_fOmega;
    BYTE m_byUseCachedGauge;
};

template<class CFieldKS>
class __DLL_EXPORT CFieldFermionKSTDR : public CFieldFermionKSTR<CFieldKS>
{
protected:
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        CFieldFermionKSTR<CFieldKS>::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        CFieldKS::_SiteKernel::FixBoundary(
            (typename CFieldKS::_Vector *)pTargetBuffer, 
            this->m_byFieldId);
    }
};


__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1R)
class CLGAPI CFieldFermionKSU1R : public CFieldFermionKSTR<CFieldFermionKSU1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSU1R)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3R)
class CLGAPI CFieldFermionKSSU3R : public CFieldFermionKSTR<CFieldFermionKSSU3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSSU3R)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3DR)
class CLGAPI CFieldFermionKSSU3DR : public CFieldFermionKSTDR<CFieldFermionKSSU3D>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSSU3DR)
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSTR_H_

//=============================================================================
// END OF FILE
//=============================================================================
