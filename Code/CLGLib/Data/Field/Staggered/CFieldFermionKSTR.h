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

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DerivateD0(pForce, pGaugeBuffer, byGaugeFieldId);

        if (m_bRealRotation)
        {
            appCrucial(_T("DerivateD0 is not supported for real rotation!\n"));
        }

        CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0_R(
            CCommonData::m_fOmega,
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
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        if (m_bRealRotation)
        {
            CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_RealRotation(
                CCommonData::m_fOmega,
                this->m_bEachSiteEta, 
                (deviceVector*)pTargetBuffer,
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
        else
        {
            CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_ImaginaryRotation(
                CCommonData::m_fOmega,
                this->m_bEachSiteEta, 
                (deviceVector*)pTargetBuffer,
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
    }

public:

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialOtherParameters(params);
        this->m_bEachSiteEta = TRUE;

        INT iReal = 0;
        if (params.FetchValueINT(_T("RealRotation"), iReal))
        {
            m_bRealRotation = (0 != iReal);
        }

        Real fValue = F(0.1);
        if (params.FetchValueReal(_T("Omega"), fValue))
        {
            m_fOmega = fValue;
        }
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GetInfos(tab);
        sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
        sRet = sRet + tab + _T("RealRotation : ") + appToString(m_bRealRotation) + _T("\n");
        return sRet;
    }

    void CopyTo(CField* U) const override
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyTo(U);
        CFieldFermionKSTR<deviceVector, deviceGauge, vectorN>* pOther = dynamic_cast<CFieldFermionKSTR<deviceVector, deviceGauge, vectorN>*>(U);
        if (NULL != pOther)
        {
            pOther->m_bRealRotation = m_bRealRotation;
            pOther->m_fOmega = m_fOmega;
        }
    }

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