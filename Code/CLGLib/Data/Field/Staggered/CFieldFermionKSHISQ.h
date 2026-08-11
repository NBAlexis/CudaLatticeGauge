//=============================================================================
// FILENAME : CFieldFermionKSHISQ.h
// 
// DESCRIPTION:
//  Only support one gauge
//  Only support even-odd in evaluation
//  
// level2:
// one-link: (1+epsi)/8
// Lepage: -1/8
// Naik:     -(1+epsi)/24
// see: 0710.0737
// default use epsi=0
// for epsi, see:
// hep-lat/0610092
// 
// 
// 
// REVISION:
//  [12/04/2022 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKST.h"
#include "CFieldFermionKSTR.h"

#ifndef _CFIELDFERMIONKSHISQ_H_
#define _CFIELDFERMIONKSHISQ_H_

__BEGIN_NAMESPACE
inline class CMultiShiftSolver* appGetMultiShiftSolver(BYTE byFieldId);

template<class CFieldKS>
class __DLL_EXPORT CFieldFermionHISQT : public CFieldKS
{
public:
    CFieldFermionHISQT()
        : CFieldKS()
        , m_fNaik(-F(0.04166666666666666666667))
        , m_fEpsilon(F(0.0))
    {

    }

protected:

    /**
    * This is only used for inverse, for example the calculation of meson correlator
    */
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        _RECORD(CFieldFermionHISQT::DOperatorKS);
        //appParanoiac(_T("CFieldFermionHISQSU3::DOperatorKS\n"));
        //This is just CFieldFermionKST::DOperatorKS
        //put effecitive gauge in
        {
            _RECORD2(CFieldFermionHISQT::DOperatorKS::DOperatorKS, a);
            CFieldKS::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        }

        {
            _RECORD2(CFieldFermionHISQT::DOperatorKS::_kernelDFermionKSNaikCached, a);
            CFieldKS::_FermionKernel::DOperatorNaik(
                this->m_bEachSiteEta,
                (typename CFieldKS::_Vector*)pTargetBuffer,
                (const typename CFieldKS::_Vector*)pBuffer,
                this->m_byFieldId,
                byGaugeFieldId,
                m_fNaik,
                m_fEpsilon * F(0.125),
                bDagger,
                eOCT,
                fRealCoeff,
                cCmpCoeff);
        }
    }

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        _RECORD(CFieldFermionHISQT::DOperatorKSOnEvenOrOdd);
        {
            _RECORD2(CFieldFermionHISQT::DOperatorKSOnEvenOrOdd::DOperatorKSOnEvenOrOdd, a);
            CFieldKS::DOperatorKSOnEvenOrOdd(pTargetBuffer, pGaugeBuffer, byGaugeFieldId, bEven, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
        }
        {
            _RECORD2(CFieldFermionHISQT::DOperatorKSOnEvenOrOdd::_kernelDFermionKSNaikCachedEvenOdd, a);
            CFieldKS::_FermionKernel::DOperatorNaikOnEvenOrOdd(
                (typename CFieldKS::_Vector*)pTargetBuffer,
                this->m_byFieldId, 
                byGaugeFieldId, 
                bEven, 
                m_fNaik,
                m_fEpsilon * F(0.125),
                bDagger, 
                eOCT,
                fRealCoeff, 
                cCmpCoeff);
        }
    }

    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override
    {
        appCrucial(_T("should never call CalculateForceEvenOddS of HISQ, because the Naik link is very special, see CActionFermionKSImprove or CActionFermionHISQCombined\n"));
    }

public:

    void InitialOtherParameters(CParameters& params) override
    {
        CFieldKS::InitialOtherParameters(params);
        //m_fNaik = F(0.0);
        UBOOL bHasNaik = params.FetchValueReal(_T("Naik"), m_fNaik);
        UBOOL bHasEpsilon = params.FetchValueReal(_T("Epsilon"), m_fEpsilon);
        if (bHasEpsilon && !bHasNaik)
        {
            m_fNaik = -(1 + m_fEpsilon) / 24;
        }
        if (!bHasEpsilon && !bHasNaik)
        {
            //hep-lat/0610092
            //Note that, m_f2am maybe set to zero and absorbed into denorminators of rational approximations
            //-(27/40/4)
            //+(327/1120/16)
            //-(15607/268800/64)
            //-(73697/3942400/256)
            m_fEpsilon = -F(0.16875) * this->m_f2am * this->m_f2am
                + F(0.018247767857142858) * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am
                - F(0.0009072149367559524) * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am
                - F(0.00007302123230773134) * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am * this->m_f2am;

            m_fNaik = -(1 + m_fEpsilon) / 24;
        }
    }

    void CopyParamTo(CField* f) const override
    {
        CFieldKS::CopyParamTo(f);
        CFieldFermionHISQT<CFieldKS>* target = dynamic_cast<CFieldFermionHISQT<CFieldKS>*>(f);
        if (NULL != target)
        {
            target->m_fNaik = m_fNaik;
            target->m_fEpsilon = m_fEpsilon;
        }
    }

    //pGauge should be effective gauge
    void CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const override
    {
        _RECORD(CFieldFermionHISQT::CalculateF0AndNaik);
        appParanoiac(_T("CFieldFermionHISQT CalculateF0AndNaik\n"));

        if (this->m_bEvenPseudofermion)
        {

            TArray<CField*> shiftsolutions;
            {
                _RECORD2(CFieldFermionHISQT::CalculateF0AndNaik::RationalApproximationPooled, a);
                this->RationalApproximationPooled(EFO_F_DDdagger, 1, 0, 0, &pGauge, NULL, NULL, this->m_iMDIndex, shiftsolutions);
            }
            for (INT i = 0; i < shiftsolutions.Num(); ++i)
            {
                CFieldFermionHISQT<CFieldKS>* fks = dynamic_cast<CFieldFermionHISQT<CFieldKS>*>(shiftsolutions[i]);
                fks->D0OnEvenOrOddS(pGauge, TRUE);

                this->CalculateForceEvenOddS_SingleTermOfRational(pGauge, fks, f0, GRASet.m_pRASet[this->m_iMDIndex]->m_lstA[i], i);
                if (abs(m_fEpsilon) > _CLG_FLT_EPSILON && NULL != pepsilonTerm)
                {
                    //For epsilon term, always use naive K-S force
                    pepsilonTerm->AddConnection(fks, GRASet.m_pRASet[this->m_iMDIndex]->m_lstA[i] * m_fEpsilon * F(0.125));
                }

                if (NULL != naik)
                {
                    CFieldKS::_FermionKernel::NaikConnection(GRASet.m_pRASet[this->m_iMDIndex]->m_lstA[i] * m_fNaik,
                        (const typename CFieldKS::_Vector*)(fks->GetData()),
                        (typename CFieldKS::_Gauge*)(naik->GetData()),
                        this->m_byFieldId);
                }

                fks->Return();
            }
        }
        else
        {
            TArray<CField*> phii;
            TArray<CFieldFermionHISQT<CFieldKS>*> phiid;
            for (UINT i = 0; i < GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree; ++i)
            {
                CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(this->m_byFieldId, _T(__FILE__), __LINE__));
                phii.AddItem(pPhi_i);
                CFieldFermionHISQT<CFieldKS>* pPhi_id = dynamic_cast<CFieldFermionHISQT<CFieldKS>*>(appGetLattice()->GetPooledFieldById(this->m_byFieldId, _T(__FILE__), __LINE__));
                phiid.AddItem(pPhi_id);
            }

            CMultiShiftSolver* solver = appGetMultiShiftSolver(this->m_byFieldId);
            if (NULL == solver)
            {
                appCrucial(_T("muitl shift solver not set! field id:%d\n"), this->m_byFieldId);
                _FAIL_EXIT;
            }
            TArray<CLGComplex> shifts;
            for (UINT i = 0; i < GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree; ++i)
            {
                shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[this->m_iMDIndex]->m_lstB[i], F(0.0)));
            }
            solver->Solve(phii, shifts, this, 1, 0, 0, &pGauge, NULL, NULL, EFO_F_DDdagger);
            //RationalApproximationPooled will allocate new fields
            //RationalApproximationPooled(EFO_F_DDdagger, 1, 0, &pGauge, NULL, m_iMDIndex, phii);

            const UINT uiBufferSize = sizeof(typename CFieldKS::_Vector*) * 2 * GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree;
            typename CFieldKS::_Vector** hostPointers = (typename CFieldKS::_Vector**)appAlloca(uiBufferSize);
            for (UINT i = 0; i < GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree; ++i)
            {
                CFieldFermionHISQT<CFieldKS>* phi_ks = dynamic_cast<CFieldFermionHISQT<CFieldKS>*>(phii[i]);
                phi_ks->FixBoundary(EFB_Field);
                phi_ks->CopyTo(phiid[i]);
                phiid[i]->D0S(pGauge);
                phiid[i]->FixBoundary(EFB_Field);
                hostPointers[i] = phi_ks->m_pDeviceData;
                hostPointers[i + GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree] = phiid[i]->m_pDeviceData;

                //do for the Naik terms
                if (NULL != naik)
                {
                    CFieldKS::_FermionKernel::NaikConnection2(GRASet.m_pRASet[this->m_iMDIndex]->m_lstA[i] * m_fNaik,
                        (const typename CFieldKS::_Vector*)(phi_ks->GetData()),
                        (const typename CFieldKS::_Vector*)(phiid[i]->GetData()),
                        (typename CFieldKS::_Gauge*)(naik->GetData()),
                        this->m_byFieldId);
                }
                if (abs(m_fEpsilon) > _CLG_FLT_EPSILON && NULL != pepsilonTerm)
                {
                    //For epsilon term, always use naive K-S force
                    CFieldKS::_FermionKernel::AddConnection2(GRASet.m_pRASet[this->m_iMDIndex]->m_lstA[i] * m_fEpsilon * F(0.125),
                        (const typename CFieldKS::_Vector*)(phi_ks->GetData()),
                        (const typename CFieldKS::_Vector*)(phiid[i]->GetData()),
                        (typename CFieldKS::_Gauge*)(pepsilonTerm->GetData()),
                        this->m_byFieldId);
                }
            }

            appSimpleCopyHD(CRationalFieldPointer::GetInstance()->GetRationPoint<typename CFieldKS::_Vector>(this->m_byRationFieldPointerBufferLength), hostPointers, uiBufferSize);

            this->DerivateD0(f0->GetData(), pGauge->GetData(), pGauge->m_byFieldId);

            for (UINT i = 0; i < GRASet.m_pRASet[this->m_iMDIndex]->m_uiDegree; ++i)
            {
                phii[i]->Return();
                phiid[i]->Return();
            }
        }
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldKS::GetInfos(tab);
        sRet = sRet + tab + _T("Naik = -(1+epsi)/24 : ") + appToString(m_fNaik) + _T("\n");
        sRet = sRet + tab + _T("Epsilon : ") + appToString(m_fEpsilon) + _T("\n");
        return sRet;
    }

    Real m_fNaik;
    Real m_fEpsilon;
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionHISQSU3)
class CLGAPI CFieldFermionHISQSU3 : public CFieldFermionHISQT<CFieldFermionKSSU3>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionHISQSU3)
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionHISQSU3R)
class CLGAPI CFieldFermionHISQSU3R : public CFieldFermionHISQT<CFieldFermionKSSU3R>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionHISQSU3R)
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSHISQ_H_

//=============================================================================
// END OF FILE
//=============================================================================