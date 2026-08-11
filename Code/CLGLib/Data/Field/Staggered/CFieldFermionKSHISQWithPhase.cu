//=============================================================================
// FILENAME : CFieldFermionKSHISQWithPhase.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [12/30/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CFieldFermionKSHISQWithPhase.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "GaugeSmearing/CGaugeSmearingHISQ.h"

__BEGIN_NAMESPACE

#define _AddEachDir(res) \
if (0 == dir) \
{ \
    _add(res[uiSiteIndex], result); \
} \
__syncthreads(); \
if (1 == dir) \
{ \
    _add(res[uiSiteIndex], result); \
} \
__syncthreads(); \
if (2 == dir) \
{ \
    _add(res[uiSiteIndex], result); \
} \
__syncthreads(); \
if (3 == dir) \
{ \
    _add(res[uiSiteIndex], result); \
} 

__CLGIMPLEMENT_CLASS(CFieldFermionHISQWithPhaseSU3)

#pragma region kernels

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikWithCharge(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const Real * __restrict__ pPhaseNaikLink,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pNaikMove,
    deviceVector* pResultData,
    Real fCharge,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalDir;

    const INT etamu = (pEtaTable[uiSiteIndex] >> dir) & 1;

    SIndex sn = pNaikMove[2 * uiLinkIndex];
    deviceGauge vn = pNaikLink[uiLinkIndex];
    Real fPhase = fCharge * pPhaseNaikLink[uiLinkIndex];
    _mul(vn, _make_cuComplex(_sin(fPhase), _cos(fPhase)));

    deviceVector result = _mulVec(vn, pDeviceData[sn.m_uiSiteIndex]);
    if (_UBOOLXOR(sn.NeedToOpposite(), etamu & 1))
    {
        _oppo(result);
    }

    sn = pNaikMove[2 * uiLinkIndex + 1];
    const UINT uiLinkIndex2 = _deviceGetLinkIndex(sn.m_uiSiteIndex, sn.m_byDir);
    vn = pNaikLink[uiLinkIndex2];
    fPhase = fCharge * pPhaseNaikLink[uiLinkIndex2];
    _mul(vn, _make_cuComplex(_sin(fPhase), _cos(fPhase)));
    if (_UBOOLXOR(sn.NeedToOpposite(), etamu & 1))
    {
        _add(result, _dagmulVec(vn, pDeviceData[sn.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(vn, pDeviceData[sn.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _AddEachDir(pResultData);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikEvenOddWithCharge(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const Real* __restrict__ pPhaseNaikLink,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pNaikMove,
    Real fCharge,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalDirEO;

    const INT etamu = (pEtaTable[uiSiteIndex] >> dir) & 1;

    SIndex sn = pNaikMove[2 * uiLinkIndex];
    deviceGauge vn = pNaikLink[uiLinkIndex];
    Real fPhase = fCharge * pPhaseNaikLink[uiLinkIndex];
    _mul(vn, _make_cuComplex(_cos(fPhase), _sin(fPhase)));

    deviceVector result = _mulVec(vn, pDeviceData[sn.m_uiSiteIndex]);
    if (_UBOOLXOR(sn.NeedToOpposite(), etamu))
    {
        _oppo(result);
    }

    sn = pNaikMove[2 * uiLinkIndex + 1];
    const UINT uiLinkIndex2 = _deviceGetLinkIndex(sn.m_uiSiteIndex, sn.m_byDir);
    vn = pNaikLink[uiLinkIndex2];
    fPhase = fCharge * pPhaseNaikLink[uiLinkIndex2];
    _mul(vn, _make_cuComplex(_cos(fPhase), _sin(fPhase)));
    if (_UBOOLXOR(sn.NeedToOpposite(), etamu))
    {
        _add(result, _dagmulVec(vn, pDeviceData[sn.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(vn, pDeviceData[sn.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    //_add(pDeviceData[uiSiteIndex], result);
    _AddEachDir(pDeviceData);
}

/**
*  |
* ||
*  |
* is like
* 
* ____
* |  |
* |__|
* 
* so just calculate all "staples" with one of the links replaced by f0
* 
* 
*   /|
*  / |
* |  |
* |  |
* x  |
*    |
*    o
* 
* 
*   /|
*  / |
* x  |
*    |
*    |
* o  |
*  \ |
*   \|
*
*
*    x
*    |
* o  |
* |  |
* |  |
*  \ |
*   \|
* 
* 
* We loop for every "long" link, and add force for its three contributions
* 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNaikForceWithCharge(
    const deviceGauge* __restrict__ f0,
    const deviceGauge* __restrict__ gauge,
    const Real* __restrict__ u1naik,
    deviceGauge* force1,
    deviceGauge* force2,
    deviceGauge* force3,
    Real fCharge,
    BYTE byGaugeFieldId)
{
    intokernalDirInt4;

    const SSmallInt4 n2 = _deviceSmallInt4OffsetC(sSite4, __fwd(dir));
    const SSmallInt4 n3 = _deviceSmallInt4OffsetC(n2, __fwd(dir));

    //const SIndex& n1__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(sSite4) + dir];
    const SIndex& n2__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n2) + dir];
    const SIndex& n3__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n3) + dir];
    const UINT linkIndex2 = _deviceGetLinkIndex(n2__mu.m_uiSiteIndex, dir);
    const UINT linkIndex3 = _deviceGetLinkIndex(n3__mu.m_uiSiteIndex, dir);

    const Real fPhase = -u1naik[uiLinkIndex] * fCharge;
    deviceGauge f0v = f0[uiLinkIndex];
    _mul(f0v, _make_cuComplex(_cos(fPhase), _sin(fPhase)));

    deviceGauge f1 = _muldagC(f0v, gauge[linkIndex3]);
    _muldag(f1, gauge[linkIndex2]);
    //_mul(f1, uphase);

    deviceGauge f2 = _dagmulC(gauge[uiLinkIndex], f0v);
    _muldag(f2, gauge[linkIndex3]);
    //_mul(f2, uphase);

    deviceGauge f3 = _mulC(gauge[uiLinkIndex], gauge[linkIndex2]);
    _dagmul(f3, f0v);
    //_mul(f3, uphase);

    force1[uiLinkIndex] = f1;
    force2[linkIndex2] = f2;
    force3[linkIndex3] = f3;
    //if (!n1__mu.IsDirichlet())
    //{
    //    _add(force[uiLinkIndex], f1);
    //}
    //__syncthreads();

    //if (!n2__mu.IsDirichlet())
    //{
    //    _add(force[linkIndex2], f2);
    //}
    //__syncthreads();

    //if (!n3__mu.IsDirichlet())
    //{
    //    _add(force[linkIndex3], f3);
    //}
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNaikForce2(
    const deviceGauge* __restrict__ f0,
    const deviceGauge* __restrict__ gauge,
    deviceGauge* force1,
    deviceGauge* force2,
    deviceGauge* force3,
    BYTE byGaugeFieldId)
{
    intokernalDirInt4;

    const BYTE fwddir = dir + 1;
    const SSmallInt4 n2 = _deviceSmallInt4OffsetC(sSite4, fwddir);
    const SSmallInt4 n3 = _deviceSmallInt4OffsetC(n2, fwddir);

    //const SIndex& n1__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(sSite4) + dir];
    const SIndex n2__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n2) + dir];
    const SIndex n3__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n3) + dir];
    const UINT linkIndex2 = _deviceGetLinkIndex(n2__mu.m_uiSiteIndex, dir);
    const UINT linkIndex3 = _deviceGetLinkIndex(n3__mu.m_uiSiteIndex, dir);

    deviceGauge f1 = _muldagC(f0[uiLinkIndex], gauge[linkIndex3]);
    _muldag(f1, gauge[linkIndex2]);

    deviceGauge f2 = _dagmulC(gauge[uiLinkIndex], f0[uiLinkIndex]);
    _muldag(f2, gauge[linkIndex3]);

    deviceGauge f3 = _mulC(gauge[uiLinkIndex], gauge[linkIndex2]);
    _dagmul(f3, f0[uiLinkIndex]);

    force1[uiLinkIndex] = f1;

    force2[linkIndex2] = f2;

    force3[linkIndex3] = f3;
}

#pragma endregion

void CFieldFermionHISQWithPhaseSU3::CopyParamTo(CField* f) const
{
    CFieldFermionHISQSU3::CopyParamTo(f);
    CFieldFermionHISQWithPhaseSU3* target = dynamic_cast<CFieldFermionHISQWithPhaseSU3*>(f);
    if (NULL != target)
    {
        target->m_fCharge = m_fCharge;
        target->m_byU1FieldId = m_byU1FieldId;
    }
}

void CFieldFermionHISQWithPhaseSU3::InitialOtherParameters(CParameters& params)
{
    CFieldFermionHISQSU3::InitialOtherParameters(params);
    params.FetchValueReal(_T("Charge"), m_fCharge);
    INT iU1FieldId = 0;
    params.FetchValueINT(_T("U1FieldId"), iU1FieldId);
    m_byU1FieldId = static_cast<BYTE>(iU1FieldId);
}

void CFieldFermionHISQWithPhaseSU3::DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    //appParanoiac(_T("CFieldFermionHISQWithPhaseSU3::DOperatorKS\n"));
    //This is just CFieldFermionKST::DOperatorKS
    //put effecitive gauge in
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    if (NULL == pU1)
    {
        appCrucial(_T("Electric Magnetic field not set correctly!\n"));
        return;
    }

    CFieldFermionKSTKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorEM(
        (deviceSU3Vector*)pTargetBuffer,
        (const deviceSU3Vector*)pBuffer,
        (const deviceSU3*)pGaugeBuffer,
        pU1->m_pDeviceData, f2am, m_fCharge,
        m_bEachSiteEta, bDagger, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);

    //const deviceSU3* pEffectiveLevel1 = (const deviceSU3* )(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
    const CGaugeSmearingHISQSU3* pSmearing = dynamic_cast<const CGaugeSmearingHISQSU3*>(appGetGaugeSmearing(byGaugeFieldId));
    const deviceSU3* pNaikLink = (const deviceSU3*)(pSmearing->GetNaikLink()->GetData());
    const Real* pNaikLinkPhase = (const Real*)(pSmearing->GetNaikLinkPhase()->GetData());
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSNaikWithCharge TMPARG(deviceSU3Vector, deviceSU3), block, threads,
        (const deviceSU3Vector*)pBuffer,
        pNaikLink,
        pNaikLinkPhase,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[m_byFieldId],
        (deviceSU3Vector*)pTargetBuffer,
        m_fCharge,
        m_byFieldId,
        byGaugeFieldId,
        bDagger,
        m_fNaik,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

void CFieldFermionHISQWithPhaseSU3::DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    if (NULL == pU1)
    {
        appCrucial(_T("Electric Magnetic field not set correctly!\n"));
        return;
    }

    CFieldFermionKSTKernel<deviceSU3Vector, deviceSU3, 3>::DOperatorEMEvenOdd(
        (deviceSU3Vector*)pTargetBuffer,
        (const deviceSU3*)pGaugeBuffer,
        pU1->m_pDeviceData, f2am, m_fCharge,
        m_bEachSiteEta, bEven, bDagger, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);

    //const deviceSU3* pEffectiveLevel1 = (const deviceSU3*)(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
    const CGaugeSmearingHISQSU3* pSmearing = dynamic_cast<const CGaugeSmearingHISQSU3*>(appGetGaugeSmearing(byGaugeFieldId));
    const deviceSU3* pNaikLink = (const deviceSU3*)(pSmearing->GetNaikLink()->GetData());
    const Real* pNaikLinkPhase = (const Real*)(pSmearing->GetNaikLinkPhase()->GetData());
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelDFermionKSNaikEvenOddWithCharge TMPARG(deviceSU3Vector, deviceSU3), block, threads,
        (deviceSU3Vector*)pTargetBuffer,
        pNaikLink,
        pNaikLinkPhase,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[m_byFieldId],
        m_fCharge,
        m_byFieldId,
        byGaugeFieldId,
        bEven,
        bDagger,
        m_fNaik,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

void CFieldFermionHISQWithPhaseSU3::CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const
{
    appCrucial(_T("Should never call CFieldFermionHISQWithPhaseSU3::CalculateForceEvenOddS!, see CFieldFermionHISQ!\n"));
    //const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    //if (NULL == pU1)
    //{
    //    appCrucial(_T("Electric Magnetic field not set correctly!\n"));
    //    return;
    //}

    //CFieldGaugeSU3* pNaikf0 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
    //CFieldGaugeSU3* pNaikForce1 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
    //CFieldGaugeSU3* pNaikForce2 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
    //CFieldGaugeSU3* pNaikForce3 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
    //pNaikf0->Zero();
    //const deviceSU3* pEffectiveLevel1 = (const deviceSU3*)(appGetGaugeSmearing(pGauge->m_byFieldId)->GetEffectiveGaugeLevel1()->GetData());
    //appParanoiac(_T("CFieldFermionHISQWithPhaseSU3 CalculateForceEvenOddS\n"));
    //TArray<CField*> shiftsolutions;
    //RationalApproximationPooled(EFO_F_DDdagger, 1, 0, &pGauge, NULL, m_iMDIndex, shiftsolutions);
    //for (INT i = 0; i < shiftsolutions.Num(); ++i)
    //{
    //    CFieldFermionHISQWithPhaseSU3* fks = dynamic_cast<CFieldFermionHISQWithPhaseSU3*>(shiftsolutions[i]);
    //    fks->D0OnEvenOrOddS(pGauge, TRUE);
    //    pForce->AddConnection(fks, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i]);

    //    NaikConnection(GRASet.m_pRASet[m_iMDIndex]->m_lstA[i],
    //        (const deviceSU3Vector*)(fks->GetData()),
    //        (deviceSU3*)(pNaikf0->GetData()),
    //        m_byFieldId);

    //    fks->Return();
    //}

    //pForce->ApplyPhaseR(pU1, -m_fCharge);
    //const CGaugeSmearingHISQSU3* pSmearing = dynamic_cast<const CGaugeSmearingHISQSU3*>(appGetGaugeSmearing(pGauge->m_byFieldId));
    //const Real* pNaikLinkPhase = (const Real*)(pSmearing->GetNaikLinkPhase()->GetData());
    //preparethreadDir;
    //_LAUNCH_KERNEL(_kernelNaikForceWithCharge<deviceSU3>, block, threads,
    //    (const deviceSU3*)(pNaikf0->GetData()),
    //    pEffectiveLevel1,
    //    pNaikLinkPhase,
    //    (deviceSU3*)(pNaikForce1->GetData()),
    //    (deviceSU3*)(pNaikForce2->GetData()),
    //    (deviceSU3*)(pNaikForce3->GetData()),
    //    m_fCharge,
    //    pGauge->m_byFieldId
    //    );

    //pNaikForce1->AxpyPlus(pNaikForce2);
    //pNaikForce1->AxpyPlus(pNaikForce3);
    //pNaikForce2->Return();
    //pNaikForce3->Return();

    //pNaikf0->Return();
    //pNaikForce1->ScalarMultply(m_fNaik);
    //pSmearing->DerivateOnU(pGauge, pNaikForce1, pForce);
    //pNaikForce1->Return();
}

void CFieldFermionHISQWithPhaseSU3::CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const
{
    _RECORD(CFieldFermionHISQWithPhaseSU3::CalculateF0AndNaik);

    if (!m_bEvenPseudofermion)
    {
        appCrucial(_T("CFieldFermionHISQWithPhaseSU3 only support even pseudo fermion now!\n"));
        _FAIL_EXIT;
    }

    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    if (NULL == pU1)
    {
        appCrucial(_T("Electric Magnetic field not set correctly!\n"));
        return;
    }

    appParanoiac(_T("CFieldFermionHISQWithPhaseSU3 CalculateF0AndNaik\n"));
    CFieldGaugeSU3* pf0 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(f0->m_byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeSU3* pnaikf0 = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(naik->m_byFieldId, _T(__FILE__), __LINE__));
    //pf0->Zero();
    pnaikf0->Zero();
    TArray<CField*> shiftsolutions;
    {
        _RECORD2(CFieldFermionHISQWithPhaseSU3::CalculateForceEvenOddS::RationalApproximationPooled, a);
        RationalApproximationPooled(EFO_F_DDdagger, 1, 0, 0, &pGauge, NULL, NULL, m_iMDIndex, shiftsolutions);
    }
    for (INT i = 0; i < shiftsolutions.Num(); ++i)
    {
        CFieldFermionHISQWithPhaseSU3* fks = dynamic_cast<CFieldFermionHISQWithPhaseSU3*>(shiftsolutions[i]);
        fks->D0OnEvenOrOddS(pGauge, TRUE);
        if (0 == i)
        {
            pf0->SetAsConnection(fks, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i]);
        }
        else
        {
            pf0->AddConnection(fks, GRASet.m_pRASet[m_iMDIndex]->m_lstA[i]);
        }

        CFieldFermionKSTKernel<deviceSU3Vector, deviceSU3, 3>::NaikConnection(GRASet.m_pRASet[m_iMDIndex]->m_lstA[i],
            (const deviceSU3Vector*)(fks->GetData()),
            (deviceSU3*)(pnaikf0->GetData()),
            m_byFieldId);

        fks->Return();
    }

    pf0->ApplyPhaseR(pU1, -m_fCharge);
    f0->AxpyPlus(pf0);
    pf0->Return();

    const CGaugeSmearingHISQSU3* pSmearing = dynamic_cast<const CGaugeSmearingHISQSU3*>(appGetGaugeSmearing(pGauge->m_byFieldId));
    pnaikf0->ApplyPhaseR(pSmearing->GetNaikLinkPhase(), -m_fCharge);
    naik->Axpy(m_fNaik, pnaikf0);
    pnaikf0->Return();
}

CCString CFieldFermionHISQWithPhaseSU3::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionHISQSU3::GetInfos(tab);
    sRet = sRet + tab + _T("UIFieldID : ") + appToString(m_byU1FieldId) + _T("\n");
    sRet = sRet + tab + _T("Charge    : ") + appToString(m_fCharge) + _T("\n");
    return sRet;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================