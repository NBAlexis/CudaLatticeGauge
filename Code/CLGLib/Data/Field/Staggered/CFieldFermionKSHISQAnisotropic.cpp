//=============================================================================
// FILENAME : CFieldFermionKSHISQAnisotropic.cpp
//
// DESCRIPTION:
// The anisotropic HISQ (aHISQ) staggered fermion field
//
// REVISION:
//  [mm/dd/yy]
//  [07/25/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldFermionKSHISQAnisotropic.h"
#include "CFieldFermionKSHISQAnisotropicKernel.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionHISQSU3Anisotropic)

void CFieldFermionHISQSU3Anisotropic::InitialOtherParameters(CParameters& params)
{
    CFieldFermionHISQSU3::InitialOtherParameters(params);

    params.FetchValueReal(_T("XiF"), m_fXiF);
    if (m_fXiF < _CLG_FLT_EPSILON)
    {
        appCrucial(_T("CFieldFermionHISQSU3Anisotropic: XiF (%f) must be positive!\n"), m_fXiF);
        _FAIL_EXIT;
    }
}

void CFieldFermionHISQSU3Anisotropic::CopyParamTo(CField* f) const
{
    CFieldFermionHISQSU3::CopyParamTo(f);
    CFieldFermionHISQSU3Anisotropic* target = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(f);
    if (NULL != target)
    {
        //every pooled field and rational shifted solution has the same operator
        target->m_fXiF = m_fXiF;
    }
}

CCString CFieldFermionHISQSU3Anisotropic::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionHISQSU3::GetInfos(tab);
    sRet = sRet + tab + _T("XiF (bare fermion anisotropy) : ") + appToString(m_fXiF) + _T("\n");
    return sRet;
}

void CFieldFermionHISQSU3Anisotropic::DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    _RECORD(CFieldFermionHISQSU3Anisotropic::DOperatorKS);
    //same structure as the parent, the temporal contribution is weighted by m_fXiF
    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorKS(m_bEachSiteEta,
        (deviceSU3Vector*)pTargetBuffer,
        (const deviceSU3Vector*)pBuffer,
        (const deviceSU3*)pGaugeBuffer,
        m_byFieldId,
        byGaugeFieldId,
        f2am,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        m_fXiF);

    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorNaik(
        m_bEachSiteEta,
        (deviceSU3Vector*)pTargetBuffer,
        (const deviceSU3Vector*)pBuffer,
        m_byFieldId,
        byGaugeFieldId,
        m_fNaik,
        m_fEpsilon * F(0.125),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        m_fXiF);
}

void CFieldFermionHISQSU3Anisotropic::DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    _RECORD(CFieldFermionHISQSU3Anisotropic::DOperatorKSOnEvenOrOdd);
    //same structure as the parent, the temporal contribution is weighted by m_fXiF
    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorKSOnEvenOrOdd(
        (deviceSU3Vector*)pTargetBuffer,
        (const deviceSU3*)pGaugeBuffer,
        m_byFieldId,
        byGaugeFieldId,
        bEven,
        f2am,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        m_fXiF);

    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorNaikOnEvenOrOdd(
        (deviceSU3Vector*)pTargetBuffer,
        m_byFieldId,
        byGaugeFieldId,
        bEven,
        m_fNaik,
        m_fEpsilon * F(0.125),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        m_fXiF);
}

void CFieldFermionHISQSU3Anisotropic::DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    //same as CFieldFermionKST::DerivateD0, the temporal contribution is weighted by m_fXiF
    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::DerivateD0(m_pDeviceData, m_byFieldId, (deviceSU3*)pForce, (const deviceSU3*)pGaugeBuffer, byGaugeFieldId,
        CRationalFieldPointer::GetInstance()->GetRationPoint<deviceSU3Vector>(m_byRationFieldPointerBufferLength),
        GRASet.m_pRASet[m_iMDIndex]->m_pDeviceData, GRASet.m_pRASet[m_iMDIndex]->m_uiDegree,
        m_fXiF);
}

void CFieldFermionHISQSU3Anisotropic::ConnectionSelf(void* res, Real fCoeff) const
{
    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::ConnectionOneFieldStaggered(
        m_pDeviceData, (deviceSU3*)res, m_byFieldId, fCoeff, m_fXiF);
}

void CFieldFermionHISQSU3Anisotropic::AddConnectionSelf(void* res, Real fCoeff) const
{
    CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::AddConnectionOneFieldStaggered(
        m_pDeviceData, (deviceSU3*)res, m_byFieldId, fCoeff, m_fXiF);
}

//pGauge should be effective gauge
//NOTE: keep in sync with CFieldFermionHISQT::CalculateF0AndNaik, the only difference is
//that NaikConnection, NaikConnection2 and AddConnection2 get m_fXiF for the temporal weight
void CFieldFermionHISQSU3Anisotropic::CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const
{
    _RECORD(CFieldFermionHISQSU3Anisotropic::CalculateF0AndNaik);
    appParanoiac(_T("CFieldFermionHISQSU3Anisotropic CalculateF0AndNaik\n"));

    if (m_bEvenPseudofermion)
    {

        TArray<CField*> shiftsolutions;
        {
            _RECORD2(CFieldFermionHISQSU3Anisotropic::CalculateF0AndNaik::RationalApproximationPooled, a);
            RationalApproximationPooled(EFO_F_DDdagger, 1, 0, 0, &pGauge, NULL, NULL, m_iMDIndex, shiftsolutions);
        }
        for (INT i = 0; i < shiftsolutions.Num(); ++i)
        {
            CFieldFermionHISQSU3Anisotropic* fks = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(shiftsolutions[i]);
            fks->D0OnEvenOrOddS(pGauge, TRUE);

            CalculateForceEvenOddS_SingleTermOfRational(pGauge, fks, f0, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i], i);
            if (abs(m_fEpsilon) > _CLG_FLT_EPSILON && NULL != pepsilonTerm)
            {
                //For epsilon term, always use naive K-S force
                pepsilonTerm->AddConnection(fks, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i] * m_fEpsilon * F(0.125));
            }

            if (NULL != naik)
            {
                CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::NaikConnection(GRASet.m_pRASet[m_iMDIndex]->m_lstA[i] * m_fNaik,
                    (const deviceSU3Vector*)(fks->GetData()),
                    (deviceSU3*)(naik->GetData()),
                    m_byFieldId,
                    m_fXiF);
            }

            fks->Return();
        }
    }
    else
    {
        TArray<CField*> phii;
        TArray<CFieldFermionHISQSU3Anisotropic*> phiid;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phii.AddItem(pPhi_i);
            CFieldFermionHISQSU3Anisotropic* pPhi_id = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
            phiid.AddItem(pPhi_id);
        }

        CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
        if (NULL == solver)
        {
            appCrucial(_T("muitl shift solver not set! field id:%d\n"), m_byFieldId);
            _FAIL_EXIT;
        }
        TArray<CLGComplex> shifts;
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[m_iMDIndex]->m_lstB[i], F(0.0)));
        }
        solver->Solve(phii, shifts, this, 1, 0, 0, &pGauge, NULL, NULL, EFO_F_DDdagger);

        const UINT uiBufferSize = sizeof(deviceSU3Vector*) * 2 * GRASet.m_pRASet[m_iMDIndex]->m_uiDegree;
        deviceSU3Vector** hostPointers = (deviceSU3Vector**)appAlloca(uiBufferSize);
        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            CFieldFermionHISQSU3Anisotropic* phi_ks = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(phii[i]);
            phi_ks->FixBoundary(EFB_Field);
            phi_ks->CopyTo(phiid[i]);
            phiid[i]->D0S(pGauge);
            phiid[i]->FixBoundary(EFB_Field);
            hostPointers[i] = phi_ks->m_pDeviceData;
            hostPointers[i + GRASet.m_pRASet[m_iMDIndex]->m_uiDegree] = phiid[i]->m_pDeviceData;

            //do for the Naik terms
            if (NULL != naik)
            {
                CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::NaikConnection2(GRASet.m_pRASet[m_iMDIndex]->m_lstA[i] * m_fNaik,
                    (const deviceSU3Vector*)(phi_ks->GetData()),
                    (const deviceSU3Vector*)(phiid[i]->GetData()),
                    (deviceSU3*)(naik->GetData()),
                    m_byFieldId,
                    m_fXiF);
            }
            if (abs(m_fEpsilon) > _CLG_FLT_EPSILON && NULL != pepsilonTerm)
            {
                //For epsilon term, always use naive K-S force
                CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>::AddConnection2(GRASet.m_pRASet[m_iMDIndex]->m_lstA[i] * m_fEpsilon * F(0.125),
                    (const deviceSU3Vector*)(phi_ks->GetData()),
                    (const deviceSU3Vector*)(phiid[i]->GetData()),
                    (deviceSU3*)(pepsilonTerm->GetData()),
                    m_byFieldId,
                    m_fXiF);
            }
        }

        appSimpleCopyHD(CRationalFieldPointer::GetInstance()->GetRationPoint<deviceSU3Vector>(m_byRationFieldPointerBufferLength), hostPointers, uiBufferSize);

        DerivateD0(f0->GetData(), pGauge->GetData(), pGauge->m_byFieldId);

        for (UINT i = 0; i < GRASet.m_pRASet[m_iMDIndex]->m_uiDegree; ++i)
        {
            phii[i]->Return();
            phiid[i]->Return();
        }
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
