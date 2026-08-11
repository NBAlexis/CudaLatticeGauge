//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3EM.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [05/01/2023 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "CFieldFermionWilsonSquareSU3EM.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3EM)

#pragma region kernels

/**
* almost copy of _kernelDFermionWilsonSquareSU3, except for the inclusion of the  U1 field
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3EM(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const Real* __restrict__ pU1Gauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fCharge,
    BYTE byFieldId,
    UBOOL bDDagger,
    UBOOL bOnlyHopping,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    //const SSmallInt4 test = __deviceSiteIndexToInt4(uiSiteIndex);

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    pResultData[uiSiteIndex] = bOnlyHopping ? result : pDeviceData[uiSiteIndex];
    if (bDDagger && !bOnlyHopping)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //=====================
        //check links
        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceWilsonVectorSU3 fermionelement = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];

        if (bDDagger)
        {
            fermionelement = gamma5.MulWilsonC(fermionelement);
        }
        //U(x,mu) phi(x+ mu)
        Real fphase = fCharge * pU1Gauge[linkIndex];
        fermionelement.MulComp(_make_cuComplex(_cos(fphase), _sin(fphase)));
        fermionelement = pGauge[linkIndex].MulWilsonVector(fermionelement);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            result.Sub(fermionelement);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(fermionelement));
        }
        else
        {
            result.Add(fermionelement);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(fermionelement));
        }

        fermionelement = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        linkIndex = _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir);
        fphase = fCharge * pU1Gauge[linkIndex];
        if (bDDagger)
        {
            fermionelement = gamma5.MulWilsonC(fermionelement);
        }
        //U^{dagger}(x-mu) phi(x-mu)
        if (x_m_mu_Gauge.NeedToDagger())
        {
            fermionelement.MulComp(_make_cuComplex(_cos(fphase), -_sin(fphase)));
            fermionelement = pGauge[linkIndex].DagMulWilsonVector(fermionelement);
        }
        else
        {
            fermionelement.MulComp(_make_cuComplex(_cos(fphase), _sin(fphase)));
            fermionelement = pGauge[linkIndex].MulWilsonVector(fermionelement);
        }
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(fermionelement);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(fermionelement));
        }
        else
        {
            result.Add(fermionelement);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(fermionelement));
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);

    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].MulReal(fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex].MulComp(cCoeff);
        break;
    default:
        break;
    }
}

/**
* if bEven, it is D on even to odd
*/
//__global__ void _CLG_LAUNCH_BOUND
//_kernelDFermionWD_DOnEvenOrOdd(
//    deviceWilsonVectorSU3* pDeviceData,
//    const deviceSU3* __restrict__ pGauge,
//    const Real* __restrict__ pU1Gauge,
//    const SIndex* __restrict__ pFermionMove,
//    const BYTE* __restrict__ pEtaTable,
//    Real fKai,
//    Real fCharge,
//    UBOOL bEven,
//    UBOOL bDDagger)
//{
//    UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) << 1U) | bEven;
//    if (uiSiteIndex >= _DC_Volume)
//    {
//        return;
//    }
//    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
//
//    // x o
//    // o x
//    // If bEven, we need uiSiteIndex - 1
//    // If not even, we need uiSiteIndex + 1
//    // So uiSiteIndex + 1 - 2 * bEven
//    const BYTE eta = pEtaTable[uiSiteIndex];
//    const BYTE mask = ((eta >> 4U) & 1U) ^ bEven;
//    uiSiteIndex += mask * (1 - (bEven << 1));
//    //eta = (mask) ? pEtaTable[uiSiteIndex] : eta;
//
//    deviceWilsonVectorSU3 result = _makeZero<deviceWilsonVectorSU3>();
//
//    UINT linkIndex = (uiSiteIndex << 2U);
//    //UINT dblinkIndex = (linkIndex << 1U);
//    //SIndex x_move_Fermion = pFermionMove[dblinkIndex];
//
//    #pragma unroll
//    for (BYTE idir = 0U; idir < 4U; ++idir)
//    {
//        const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];
//        UINT dblinkIndex = (linkIndex << 1U);
//        SIndex x_move_Fermion = pFermionMove[dblinkIndex];
//        Real fU1Phase = fCharge * pU1Gauge[linkIndex];
//        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
//        u_dagger_phi_x_m_m.MulComp(_make_cuComplex(_cos(fU1Phase), _sin(fU1Phase)));
//        if (bDDagger)
//        {
//            u_dagger_phi_x_m_m = gamma5.MulWilsonC(u_dagger_phi_x_m_m);
//        }
//        if (x_move_Fermion.NeedToOpposite())
//        {
//            _sub(result, u_dagger_phi_x_m_m);
//            _add(result, gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
//        }
//        else
//        {
//            _add(result, u_dagger_phi_x_m_m);
//            _sub(result, gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
//        }
//
//        x_move_Fermion = pFermionMove[dblinkIndex | 1];
//        dblinkIndex = _deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir);
//        fU1Phase = -fCharge * pU1Gauge[dblinkIndex];
//        u_dagger_phi_x_m_m = _dagmulVec(pGauge[dblinkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
//        u_dagger_phi_x_m_m.MulComp(_make_cuComplex(_cos(fU1Phase), _sin(fU1Phase)));
//        if (bDDagger)
//        {
//            u_dagger_phi_x_m_m = gamma5.MulWilsonC(u_dagger_phi_x_m_m);
//        }
//        if (x_move_Fermion.NeedToOpposite())
//        {
//            _sub(result, u_dagger_phi_x_m_m);
//            _sub(result, gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
//        }
//        else
//        {
//            _add(result, u_dagger_phi_x_m_m);
//            _add(result, gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
//        }
//        ++linkIndex;
//    }
//
//    if (bDDagger)
//    {
//        result = gamma5.MulWilsonC(result);
//    }
//
//    pDeviceData[uiSiteIndex].Add(result);
//}
//

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3EM(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const Real* __restrict__ pU1Gauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    Real fCharge,
    BYTE byFieldId)
{
    intokernalDir;

    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    const deviceWilsonVectorSU3& x_Right = pInverseD[uiSiteIndex];

    //idir = mu
        //Get Gamma mu
    const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + dir];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

    const deviceWilsonVectorSU3& x_p_mu_Right = pInverseD[x_p_mu_Fermion.m_uiSiteIndex];
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const Real fU1Phase = fCharge * pU1Gauge[linkIndex];
    const Real cs = _cos(fU1Phase);
    const Real sn = _sin(fU1Phase);

    deviceWilsonVectorSU3 right1(x_p_mu_Right);
    right1.Sub(gammaMu.MulWilsonC(right1));
    right1.MulComp(_make_cuComplex(cs, sn));
    deviceSU3 mid = deviceSU3::makeSU3Contract(right1, x_Left);

    deviceWilsonVectorSU3 right2(x_Right);
    right2.Add(gammaMu.MulWilsonC(right2));
    right2.MulComp(_make_cuComplex(cs, -sn));
    mid.Add(deviceSU3::makeSU3Contract(x_p_mu_Left, right2));

    _mul(mid, fKai * (1 - 2 * static_cast<INT>(x_p_mu_Fermion.NeedToOpposite())));
    pForce[linkIndex].Add(mid);
}

#pragma endregion


#pragma region DOperator

void CFieldFermionWilsonSquareSU3EM::DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    _RECORD(CFieldFermionWilsonSquareSU3EM::DOperator);
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const CFieldGaugeU1Real* pU1GaugeField = static_cast<const CFieldGaugeU1Real * >(appGetLattice()->GetFieldById(m_byU1FieldId));
    const Real* pGaugeU1 = (const Real*)pU1GaugeField->m_pDeviceData;

    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionWilsonSquareSU3EM, block, threads, 
        pSource,
        pGauge,
        pGaugeU1,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        pTarget,
        static_cast<Real>(m_fKai),
        m_fCharge,
        m_byFieldId,
        bDagger,
        FALSE,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionWilsonSquareSU3EM::DerivateDOperator(DOUBLE fCoeff, void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    _RECORD(CFieldFermionWilsonSquareSU3EM::DerivateDOperator);
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    //const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;
    const CFieldGaugeU1Real* pU1GaugeField = static_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    const Real* pGaugeU1 = (const Real*)pU1GaugeField->m_pDeviceData;

    //appCrucial(_T("need to change according to new implementation!\n"));

    preparethreadDir;
    _LAUNCH_KERNEL(_kernelDWilsonForceSU3EM, block, threads, 
        pDphiBuffer,
        pDDphiBuffer,
        //pGauge,
        pGaugeU1,
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        pForceSU3,
        static_cast<Real>(fCoeff),
        m_fCharge, 
        m_byFieldId);

}

//UBOOL CFieldFermionWilsonSquareSU3EM::CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const
//{
//    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
//        return FALSE;
//    }
//    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
//        return FALSE;
//    }
//
//    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
//    const CFieldGaugeU1Real* pGaugeU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
//
//    if (ER_NoRational == m_eRational)
//    {
//        CField* pDDaggerPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
//        CField* pDPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
//        CField* pCachedField = CCommonData::m_bStoreLastSolution ?
//            appGetLattice()->m_pFieldCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField)
//            : NULL;
//
//        if (NULL == pDDaggerPhi || EFT_FermionWilsonSquareSU3 != pDDaggerPhi->GetFieldType()
//            || NULL == pDPhi || EFT_FermionWilsonSquareSU3 != pDPhi->GetFieldType())
//        {
//            appCrucial(_T("Pooled field not found!\n"));
//            if (NULL != pDDaggerPhi)
//            {
//                pDDaggerPhi->Return();
//            }
//            if (NULL != pDPhi)
//            {
//                pDPhi->Return();
//            }
//            return FALSE;
//        }
//        CFieldFermionWilsonSquareSU3EM* pDDaggerPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3EM*>(pDDaggerPhi);
//        CFieldFermionWilsonSquareSU3EM* pDPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3EM*>(pDPhi);
//        //if (!pDDaggerPhiWilson->InverseDDdagger(pGaugeSU3))
//
//        //if (m_byEvenFieldId > 0)
//        //{
//        //    CopyTo(pDDaggerPhiWilson);
//        //    pDDaggerPhiWilson->InverseDDdagger(pGaugeSU3);
//        //}
//        //else
//        //{
//        TArray<const CFieldGauge*> gauge;
//        gauge.AddItem(pGaugeSU3);
//        if (!appGetFermionSolver(m_byFieldId)->Solve(
//            pDDaggerPhiWilson, this, 1, 0, gauge.GetData(), NULL,
//            EFO_F_DDdagger, ePhase, pCachedField))
//        {
//            appCrucial(_T("Sparse Linear Solver failed...\n"));
//            pDDaggerPhi->Return();
//            pDPhi->Return();
//            return FALSE;
//        }
//        //}
//
//        //phi 2 = D^{-1}phi = D+ (DD+)^{-1} phi
//        //It is faster to calcuate D+ phi2 then D^{-1} phi
//        pDDaggerPhiWilson->CopyTo(pDPhiWilson);
//        if (NULL != pCachedField)
//        {
//            //The gauge field is changing slowly, and D depends only on gauge, also change slowly
//            //Use the last solution as start point will accelerate the solver, so we cache it
//            pDDaggerPhiWilson->CopyTo(pCachedField);
//        }
//        pDPhiWilson->DdaggerS(pGaugeSU3);
//
//        CFieldGaugeSU3* pToAddForce = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGaugeSU3->m_byFieldId));
//
//        DerivateDOperator(
//            pToAddForce->m_pDeviceData,
//            pDPhiWilson->m_pDeviceData,
//            pDDaggerPhiWilson->m_pDeviceData,
//            pGaugeSU3->m_pDeviceData,
//            pGaugeSU3->m_byFieldId);
//
//        pToAddForce->ApplyPhaseR(pGaugeU1, m_fCharge);
//        pToAddForce->LeftMul(pGaugeSU3);
//        pToAddForce->TA();
//        pForce->Axpy(m_fKai, pToAddForce);
//
//        pToAddForce->Return();
//        pDDaggerPhi->Return();
//        pDPhi->Return();
//    }
//    else if (ER_WDMDRational == m_eRational)
//    {
//        appCrucial(_T("not supported yet"));
//
//    }
//    else //all rational mode
//    {
//        TArray<CField*> shiftsolutions;
//        CFieldFermionWilsonSquareSU3EM* lhs = dynamic_cast<CFieldFermionWilsonSquareSU3EM*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//        CFieldGaugeSU3* pToAddForceAll = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGaugeSU3->m_byFieldId));
//        CFieldGaugeSU3* pToAddForce = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetPooledFieldById(pGaugeSU3->m_byFieldId));
//        pToAddForceAll->Zero();
//
//        //this find [1/(D^+D+b)]phi
//        RationalApproximationPooled(EFO_F_DDdagger, 1, 0, &pGauge, NULL, m_iMDIndex, shiftsolutions);
//        for (INT i = 0; i < shiftsolutions.Num(); ++i)
//        {
//            CFieldFermionWilsonSquareSU3EM* fks = dynamic_cast<CFieldFermionWilsonSquareSU3EM*>(shiftsolutions[i]);
//            //Apply only d0
//            fks->CopyTo(lhs);
//            lhs->DdaggerS(pGaugeSU3);
//
//            DerivateDOperator(
//                pToAddForce->m_pDeviceData,
//                lhs->m_pDeviceData,
//                fks->m_pDeviceData,
//                pGaugeSU3->m_pDeviceData,
//                pGaugeSU3->m_byFieldId);
//
//            //preparethread;
//            //_LAUNCH_KERNEL(_kernelDWilsonForceSU3, block, threads, 
//            //    fks->m_pDeviceData,
//            //    rhs->m_pDeviceData,
//            //    //pGauge,
//            //    appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
//            //    pToAddForce->m_pDeviceData,
//            //    m_byFieldId);
//
//            pToAddForceAll->Axpy(m_fKai * GRASet.m_pRASet[m_iMDIndex]->m_lstA[i], pToAddForce);
//
//            fks->Return();
//        }
//        lhs->Return();
//        pToAddForce->Return();
//
//        pToAddForceAll->ApplyPhaseR(pGaugeU1, -m_fCharge);
//        pToAddForceAll->LeftMul(pGaugeSU3);
//        pToAddForceAll->TA();
//        pForce->AxpyPlus(pToAddForceAll);
//        pToAddForceAll->Return();
//    }
//
//    return TRUE;
//}

/*
void CFieldFermionWilsonSquareSU3EM::DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) 
{
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    //if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    //{
    //    appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
    //    return;
    //}
    //_RECORD(CFieldFermionKST::DDdaggerS);

    //const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGaugeBuffer = (const deviceSU3*)pFieldSU3->m_pDeviceData;
    const CFieldGaugeU1Real* pU1GaugeField = static_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    const Real* pGaugeU1 = (const Real*)pU1GaugeField->m_pDeviceData;

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    //TODO make even
    if (TRUE)
    {
        {
            //Add is used, so zero on odd
            _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::ZeroOnEvenOdd, a);
            ZeroOnEvenOdd(FALSE);
        }
        //appParanoiac(_T("DDdaggerS on even\n"));
        {
            _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::DOperatorKSOnEvenOrOdd, a);
            _kernelDFermionWD_DOnEvenOrOdd(m_pDeviceData, pGaugeBuffer, pGaugeU1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId], 
                appGetLattice()->m_pIndexCache->m_pEtaMu, 
                m_fKai, F(0.0), TRUE, FALSE);
        }
        {
            _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::DOperatorKSOnEvenOrOdd, a);
            _kernelDFermionWD_DOnEvenOrOdd(m_pDeviceData, pGaugeBuffer, pGaugeU1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                m_fKai, F(0.0), FALSE, TRUE);
        }
        {
            _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::ZeroOnEvenOdd, a);
            ZeroOnEvenOdd(FALSE);
        }
    }
    //else
    //{
    //    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    //    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
    //        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //    //why only apply coeff in the next step?
    //    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
    //        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    //    pPooled->Return();
    //}
}
*/

#pragma endregion

CFieldFermionWilsonSquareSU3EM::CFieldFermionWilsonSquareSU3EM()
    : CFieldFermionWilsonSquareSU3()
    , m_fCharge(F(0.0))
    , m_byU1FieldId(0)
{
    
}

CFieldFermionWilsonSquareSU3EM::~CFieldFermionWilsonSquareSU3EM()
{
    
}

void CFieldFermionWilsonSquareSU3EM::InitialOtherParameters(CParameters & params)
{
    CFieldFermionWilsonSquareSU3::InitialOtherParameters(params);

    Real fValue = F(0.0);
    if (params.FetchValueReal(_T("Charge"), fValue))
    {
        m_fCharge = fValue;
    }

    INT iU1FieldId = 0;
    if (params.FetchValueINT(_T("EMFieldID"), iU1FieldId))
    {
        m_byU1FieldId = static_cast<BYTE>(iU1FieldId);
    }
}

//void CFieldFermionWilsonSquareSU3EM::PrepareForHMCS(const CFieldGauge* pGauge)
//{
//    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
//        return;
//    }
//    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
//    const CFieldGaugeU1Real* pU1GaugeField = static_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
//    const Real* pGaugeU1 = (const Real*)pU1GaugeField->m_pDeviceData;
//
//    preparethread;
//    InitialField(EFIT_RandomGaussian);
//
//    {
//        _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::DOperatorKSOnEvenOrOdd, a);
//        _kernelDFermionWD_DOnEvenOrOdd(m_pDeviceData, pFieldSU3->m_pDeviceData, pGaugeU1,
//            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
//            appGetLattice()->m_pIndexCache->m_pEtaMu,
//            m_fKai, F(0.0), FALSE, TRUE);
//    }
//    {
//        _RECORD2(CFieldFermionWilsonSquareSU3EM::DDdaggerS::ZeroOnEvenOdd, a);
//        ZeroOnEvenOdd(FALSE);
//    }
//
//
//    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
//    {
//        m_fLength = Dot(this).x;
//    }
//
//    CCommonData::m_bStoreLastSolution = FALSE;
//}

void CFieldFermionWilsonSquareSU3EM::CopyParamTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3::CopyParamTo(U);
    CFieldFermionWilsonSquareSU3EM* pOther = dynamic_cast<CFieldFermionWilsonSquareSU3EM*>(U);

    pOther->m_fCharge = m_fCharge;
    pOther->m_byU1FieldId = m_byU1FieldId;
}

CCString CFieldFermionWilsonSquareSU3EM::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionWilsonSquareSU3::GetInfos(tab);
    sRet = sRet + tab + _T("Charge : ") + appToString(m_fCharge) + _T("\n");
    sRet = sRet + tab + _T("EMFieldID : ") + appToString(m_byU1FieldId) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================