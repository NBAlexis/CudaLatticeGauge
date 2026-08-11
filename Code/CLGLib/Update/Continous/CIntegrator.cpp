//=============================================================================
// FILENAME : CIntegrator.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [mm/dd/yy]
//  [12/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CIntegratorLeapFrog.h"
#include "Update/CStapleCache.h"

__BEGIN_NAMESPACE

CIntegrator::~CIntegrator()
{
    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        appSafeDelete(m_pGaugeField[i]);
        appSafeDelete(m_pForceField[i]);
        appSafeDelete(m_pMomentumField[i]);
    }

    //for (INT i = 0; i < m_pStapleField.Num(); ++i)
    //{
    //    appSafeDelete(m_pStapleField[i]);
    //}

    for (INT i = 0; i < m_pBosonFields.Num(); ++i)
    {
        appSafeDelete(m_pBosonFields[i]);
        appSafeDelete(m_pBosonForceFields[i]);
        appSafeDelete(m_pBosonMomentumFields[i]);
    }

    for (INT i = 0; i < m_pTensor2Field.Num(); ++i)
    {
        appSafeDelete(m_pTensor2Field[i]);
    }

    for (INT i = 0; i < m_pUPrime.Num(); ++i)
    {
        appSafeDelete(m_pUPrime[i]);
    }

    for (INT i = 0; i < m_pPhiPrime.Num(); ++i)
    {
        appSafeDelete(m_pPhiPrime[i]);
    }
}

/**
* Create the fields here
*/
void CIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    m_pOwner = pOwner;
    m_pLattice = pLattice;
    m_lstActions = pLattice->m_pActionList;

    //m_bStapleCached = FALSE;

    INT iStepCount = 50;
    params.FetchValueINT(_T("IntegratorStep"), iStepCount);
#if _CLG_DEBUG
    params.FetchValueINT(_T("IntegratorStepDebug"), iStepCount);
#endif

    Real fStepLength = F(1.0);
    params.FetchValueReal(_T("IntegratorStepLength"), fStepLength);
    m_uiStepCount = static_cast<UINT>(iStepCount);
    m_uiStepCountMetropolis = m_uiStepCount;
    m_fEStep = fStepLength / m_uiStepCount;
    iStepCount = 0;
    params.FetchValueINT(_T("IntegratorStepWarmup"), iStepCount);
    m_uiStepCountWarmup = static_cast<UINT>(iStepCount);

    INT iDebugForce = 0;
    params.FetchValueINT(_T("DebugForce"), iDebugForce);
    m_bDebugForce = (0 != iDebugForce);

    INT iBindDir = 0;
    params.FetchValueINT(_T("BindDir"), iBindDir);
    m_byBindDir = static_cast<BYTE>(iBindDir);

    for (INT i = 0; i < pLattice->m_pGaugeField.Num(); ++i)
    {
        if (pLattice->m_pGaugeField[i]->IsDynamic())
        {
            m_pGaugeField.AddItem(dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField[i]->GetCopy()));
            m_pGaugeField[i]->InitialField(EFIT_Zero);

            m_pForceField.AddItem(dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField[i]->GetCopy()));
            m_pForceField[i]->InitialField(EFIT_Zero);

            m_pMomentumField.AddItem(dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField[i]->GetCopy()));
            m_pMomentumField[i]->InitialField(EFIT_Zero);

            //if (CCommonData::m_bStoreStaple)
            //{
            //    m_pStapleField.AddItem(dynamic_cast<CFieldGauge*>(pLattice->m_pGaugeField[i]->GetCopy()));
            //    m_pStapleField[i]->InitialField(EFIT_Zero);
            //}
        }
        else
        {
            m_pGaugeField.AddItem(NULL);
            m_pForceField.AddItem(NULL);
            m_pMomentumField.AddItem(NULL);
            //if (CCommonData::m_bStoreStaple)
            //{
            //    m_pStapleField.AddItem(NULL);
            //}
        }
    }

    for (INT i = 0; i < pLattice->m_pBosonField.Num(); ++i)
    {
        if (pLattice->m_pBosonField[i]->IsDynamic())
        {
            m_pBosonFields.AddItem(dynamic_cast<CFieldBoson*>(pLattice->m_pBosonField[i]->GetCopy()));
            m_pBosonFields[i]->InitialField(EFIT_Zero);

            m_pBosonForceFields.AddItem(dynamic_cast<CFieldBoson*>(pLattice->m_pBosonField[i]->GetCopy()));
            m_pBosonForceFields[i]->InitialField(EFIT_Zero);

            m_pBosonMomentumFields.AddItem(dynamic_cast<CFieldBoson*>(pLattice->m_pBosonField[i]->GetCopy()));
            m_pBosonMomentumFields[i]->InitialField(EFIT_Zero);
        }
        else
        {
            m_pBosonFields.AddItem(NULL);
            m_pBosonForceFields.AddItem(NULL);
            m_pBosonMomentumFields.AddItem(NULL);
        }
    }

    //the tensor2 fields are working fields, they can be changed by the actions in OnFinishTrajectory
    for (INT i = 0; i < pLattice->m_pTensor2Field.Num(); ++i)
    {
        if (pLattice->m_pTensor2Field[i]->IsDynamic())
        {
            m_pTensor2Field.AddItem(dynamic_cast<CFieldTensor2*>(pLattice->m_pTensor2Field[i]->GetCopy()));
            m_pTensor2Field[i]->InitialField(EFIT_Zero);
        }
        else
        {
            m_pTensor2Field.AddItem(NULL);
        }
    }

    CCString sBackupFieldType;
    params.FetchStringValue(_T("BackupFieldTypeName"), sBackupFieldType);
    if (0 == sBackupFieldType.CompareNoCase(_T("CFieldGaugeSU3_12")))
    {
        m_eBackupFieldType = EBFT_SU3_12;
    }
    else
    {
        m_eBackupFieldType = EBFT_Same;
    }

    OnGaugeChanged();
}

void CIntegrator::Prepare(UBOOL bLastAccepted, UINT uiStep)
{
    //we may not accept the evaluation, so we need to copy it first
    if (!bLastAccepted || 0 == uiStep)
    {
        for (INT i = 0; i < m_pLattice->m_pGaugeField.Num(); ++i)
        {
            if (NULL != m_pGaugeField[i])
            {
                m_pLattice->m_pGaugeField[i]->CopyTo(m_pGaugeField[i]);
                m_pGaugeField[i]->SetOneDirectionUnity(m_byBindDir);
            }
        }
        for (INT i = 0; i < m_pLattice->m_pBosonField.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pLattice->m_pBosonField[i]->CopyTo(m_pBosonFields[i]);
            }
        }
        for (INT i = 0; i < m_pLattice->m_pTensor2Field.Num(); ++i)
        {
            if (NULL != m_pTensor2Field[i])
            {
                m_pLattice->m_pTensor2Field[i]->CopyTo(m_pTensor2Field[i]);
            }
        }
        OnGaugeChanged();

        //m_bStapleCached = FALSE;
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }

    OnCacheAndSmearing(3);

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->PrepareForHMC(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), uiStep);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }

    //generate a random momentum field to start
    InitialMomentumNoise();
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CIntegrator::RequireGaugeSmearing() 
{
    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        if (NULL != m_pGaugeField[i])
        {
            CGaugeSmearing* smearing = appGetGaugeSmearing(m_pGaugeField[i]->m_byFieldId);
            if (NULL != smearing && smearing->CalledWhenUpdate())
            {
                smearing->GaugeSmearingC(m_pGaugeField[i]);
            }
        }
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CIntegrator::OnCacheStaple(ECacheCall eCall)
{
    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        if (NULL != m_pGaugeField[i])
        {
            CStapleCache* cache = appGetStapleCache(m_pGaugeField[i]->m_byFieldId);
            if (NULL != cache)
            {
                cache->Cache(m_pGaugeField[i], eCall);
            }
        }
    }
}

void CIntegrator::OnFinishTrajectory(UBOOL bAccepted)
{
    //give the actions a chance to modify the fields which will be accepted
    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->OnFinishTrajectory(bAccepted, m_pGaugeField.Num(), m_pBosonFields.Num(), m_pTensor2Field.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pTensor2Field.GetData());
    }
    if (bAccepted)
    {
        for (INT i = 0; i < m_pLattice->m_pGaugeField.Num(); ++i)
        {
            if (NULL != m_pGaugeField[i])
            {
                m_pGaugeField[i]->CopyTo(m_pLattice->m_pGaugeField[i]);
            }
        }
        for (INT i = 0; i < m_pLattice->m_pBosonField.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pBosonFields[i]->CopyTo(m_pLattice->m_pBosonField[i]);
            }
        }
        for (INT i = 0; i < m_pLattice->m_pTensor2Field.Num(); ++i)
        {
            if (NULL != m_pTensor2Field[i])
            {
                m_pTensor2Field[i]->CopyTo(m_pLattice->m_pTensor2Field[i]);
            }
        }
    }
    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->OnFinishTrajectory(bAccepted);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CIntegrator::UpdateU(Real fStep)
{
    _RECORD(CIntegrator::UpdateU);
    FixGaugeBondary(m_pMomentumField, EFB_Momentum);
    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        if (NULL != m_pGaugeField[i])
        {
            m_pMomentumField[i]->SetOneDirectionZero(m_byBindDir);

            //U(k) = exp (i e P) U(k-1)
            if (abs(_HC_GaugeMomentumFactor - F(1.0)) > _CLG_FLT_EPSILON)
            {
                m_pMomentumField[i]->ExpMult(fStep / _HC_GaugeMomentumFactor, m_pGaugeField[i]);
            }
            else
            {
                m_pMomentumField[i]->ExpMult(fStep, m_pGaugeField[i]);
            }

            m_pGaugeField[i]->SetOneDirectionUnity(m_byBindDir);
        }
    }
    FixGaugeBondary(m_pGaugeField, EFB_Field);

    FixBosonBondary(m_pBosonMomentumFields, EFB_Momentum);
    for (INT i = 0; i < m_pBosonFields.Num(); ++i)
    {
        if (NULL != m_pBosonFields[i])
        {
            //m_pBosonMomentumFields[i]->DebugPrintMe();
            m_pBosonFields[i]->Axpy(fStep, m_pBosonMomentumFields[i]);
        }
    }
    FixBosonBondary(m_pBosonFields, EFB_Field);

    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    OnGaugeChanged();
}

void CIntegrator::UpdateP(Real fStep, ESolverPhase ePhase)
{
    _RECORD(CIntegrator::UpdateP);
    // recalc force
    
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    SetOneDirOne(m_pGaugeField, m_byBindDir);

    CalcForceOfActions(m_lstActions, EFC_All, ePhase);
    //P = P + e F
    AddForce(fStep, TRUE);

    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CIntegrator::CalcForceOfActions(const TArray<CAction*>& actionlst, EForceCalc eMode, ESolverPhase ePhase)
{
    ZeroForce();
    if (EFC_All == eMode)
    {
        OnCacheAndSmearing(3);
        for (INT i = 0; i < actionlst.Num(); ++i)
        {
            actionlst[i]->CalculateForce(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pForceField.GetData(), m_pBosonForceFields.GetData(),
                NULL, ePhase);
            //checkCudaErrors(cudaDeviceSynchronize());
            //checkCudaErrors(cudaGetLastError());
            //_CHECKCUDA;
        }
        //m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    }
    else if (EFC_Fermion == eMode)
    {
        OnCacheAndSmearing(2);
        for (INT i = 0; i < actionlst.Num(); ++i)
        {
            if (actionlst[i]->IsFermion())
            {
                actionlst[i]->CalculateForce(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pForceField.GetData(), m_pBosonForceFields.GetData(), NULL, ePhase);
                //checkCudaErrors(cudaDeviceSynchronize());
                //checkCudaErrors(cudaGetLastError());
                //_CHECKCUDA;
            }
        }
    }
    else if (EFC_Gauge == eMode)
    {
        OnCacheAndSmearing(1);
        for (INT i = 0; i < actionlst.Num(); ++i)
        {
            if (!actionlst[i]->IsFermion())
            {
                actionlst[i]->CalculateForce(m_pGaugeField.Num(), m_pBosonFields.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pForceField.GetData(), m_pBosonForceFields.GetData(),
                    NULL, ESP_Once);
                //checkCudaErrors(cudaDeviceSynchronize());
                //checkCudaErrors(cudaGetLastError());
                _CHECKCUDA;
            }
        }
        //m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    }

    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        if (NULL != m_pGaugeField[i] && NULL != m_pForceField[i])
        {
            // f = U f^+
            m_pForceField[i]->LeftMul(m_pGaugeField[i], FALSE, TRUE);
            m_pForceField[i]->TA();
        }
    }
    _CHECKCUDA;
    SetOneDirZero(m_pForceField, m_byBindDir);
    _CHECKCUDA;
}

void CIntegrator::FinishEvaluate()
{
    for (INT i = 0; i < m_pGaugeField.Num(); ++i)
    {
        if (NULL != m_pGaugeField[i])
        {
            m_pGaugeField[i]->ElementNormalize();
            m_pGaugeField[i]->SetOneDirectionUnity(m_byBindDir);
        }
    }
    OnGaugeChanged();
}

DOUBLE CIntegrator::GetEnergy(UBOOL bBeforeEvolution, TArray<DOUBLE>& actions)
{
    SetOneDirZero(m_pMomentumField, m_byBindDir);
    SetOneDirOne(m_pGaugeField, m_byBindDir);

    DOUBLE retv = CalcMomentumEnery();

    CCString sLog = _T("");
    sLog.Format(_T("kin:%f, "), retv);
    actions.AddItem(retv);

    m_bUDirty[0] = TRUE;
    m_bUDirty[1] = TRUE;
    OnCacheAndSmearing(3);
    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        DOUBLE fActionEnergy = m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField.Num(), m_pBosonFields.Num(), m_pTensor2Field.Num(), m_pGaugeField.GetData(), m_pBosonFields.GetData(), m_pTensor2Field.GetData(), NULL);

        CCString sThisActionInfo = _T("");
        sThisActionInfo.Format(_T(" Action%d:%f, "), i + 1, fActionEnergy);
        sLog += sThisActionInfo;
        retv += fActionEnergy;
        actions.AddItem(fActionEnergy);
    }

    appGeneral(_T("H (%s) = %s \n"), bBeforeEvolution ? "before" : "after" , sLog.c_str());
    return retv;
}

void CNestedIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);

    INT iNestedStepCount = 3;
    params.FetchValueINT(_T("NestedStep"), iNestedStepCount);
    if (iNestedStepCount < 1)
    {
        appCrucial(_T("NestedStep must >= 1, but set to be %d!\n"), iNestedStepCount)
    }
    m_uiNestedStep = static_cast<UINT>(iNestedStepCount);
    m_fNestedStepLength = m_fEStep / m_uiNestedStep;

    INT iInnerLeapfrog = 0;
    params.FetchValueINT(_T("InnerLeapfrog"), iInnerLeapfrog);
    m_bInnerLeapFrog = (0 != iInnerLeapfrog);
}

CCString CNestedIntegrator::GetNestedInfo(const CCString & sTab) const
{
    return sTab + _T("Nested : ") + appToString(static_cast<INT>(m_uiNestedStep)) + _T("\n")
         + sTab + _T("InnerLeapFrog : ") + (m_bInnerLeapFrog ? _T("1") : _T("0")) + _T("\n");
}

void CNestedIntegrator::UpdatePF(Real fStep, ESolverPhase ePhase)
{
    _RECORD(CNestedIntegrator::UpdatePF);

    CalcForceOfActions(m_lstActions, EFC_Fermion, ePhase);
    _CHECKCUDA;

    //P = P + e F
    AddForce(fStep, FALSE);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    if (m_bDebugForce)
    {
        appGeneral(_T(" ------ Fermion Force= %f \n"), CalcForce());
    }
}

void CNestedIntegrator::UpdatePG(Real fStep)
{
    _RECORD(CNestedIntegrator::UpdatePG);

    CalcForceOfActions(m_lstActions, EFC_Gauge, ESP_Once);
    //P = P + e F
    AddForce(fStep, FALSE);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    if (m_bDebugForce)
    {
        appGeneral(_T(" ------ Gauge Force= %f \n"), CalcForce());
    }
}

void CNestedIntegrator::NestedEvaluateLeapfrog(UBOOL bLast)
{
    const Real fHalfPstep = F(0.5) * m_fNestedStepLength;
    UpdatePG(fHalfPstep);
    appParanoiac("  leap frog nested sub step 0\n");

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(m_fNestedStepLength);

        if (uiStep < m_uiNestedStep)
        {
            UpdatePG(m_fNestedStepLength);
            appParanoiac("  leap frog nested sub step %d\n", uiStep);
        }
        else
        {
            UpdatePG(fHalfPstep);
            appParanoiac("  leap frog nested last step %d\n", uiStep);
        }
    }
}

void CMultiLevelNestedIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);

    TArray<UINT> nestedSteps;
    params.FetchValueArrayUINT(_T("NestedSteps"), nestedSteps);
    if (nestedSteps.Num() < 1)
    {
        appCrucial(_T("NestedSteps.Num must >= 1, but set to be 0!\n"));
        _FAIL_EXIT;
    }

    for (INT i = 0; i < nestedSteps.Num(); ++i)
    {
        if (nestedSteps[i] < 1)
        {
            appCrucial(_T("NestedSteps must >= 1, but set to be 0!\n"));
            _FAIL_EXIT;
        }
    }

    m_uiNestedStep = nestedSteps;

    m_fNestedStepLengths.RemoveAll();
    m_fNestedStepLengths.AddItem(m_fEStep);
    Real fStep = m_fEStep;
    for (INT i = 0; i < nestedSteps.Num(); ++i)
    {
        fStep = fStep / nestedSteps[i];
        m_fNestedStepLengths.AddItem(fStep);
    }

    m_iNestedActionId.RemoveAll();
    for (INT i = 0; i <= nestedSteps.Num(); ++i)
    {
        TArray<UINT> actionlist;
        params.FetchValueArrayUINT(_T("NestedActionList") + appToString(i), actionlist);
        if (actionlist.Num() < 1)
        {
            appCrucial(_T("NestedActionList.Num must >= 1, but set to be 0!\n"));
            _FAIL_EXIT;
        }
        m_iNestedActionId.AddItem(actionlist);
        if (i == nestedSteps.Num())
        {
            if (actionlist.Num() > 2 || actionlist[0] != 0)
            {
                appGeneral(_T("The last nested action list should only have Gauge action, but not!\n"));
            }
        }
    }

    m_fTotalStepLength = F(1.0);
    params.FetchValueReal(_T("IntegratorStepLength"), m_fTotalStepLength);
    INT iInnerLeapfrog = 0;
    params.FetchValueINT(_T("InnerLeapfrog"), iInnerLeapfrog);
    m_bInnerLeapFrog = (0 != iInnerLeapfrog);
}

CCString CMultiLevelNestedIntegrator::GetNestedInfo(const CCString& sTab) const
{
    CCString sRet;
    for (INT i = 0; i <= m_uiNestedStep.Num(); ++i)
    {
        const INT iStep = (0 == i) ? static_cast<INT>(m_uiStepCount) : static_cast<INT>(m_uiNestedStep[i - 1]);
        CCString sActionList = _T("[");
        for (INT j = 0; j < m_iNestedActionId[i].Num(); ++j)
        {
            sActionList = sActionList + appToString(static_cast<INT>(m_iNestedActionId[i][j]));
            if (j != m_iNestedActionId[i].Num() - 1)
            {
                sActionList = sActionList + _T(", ");
            }
        }
        sActionList = sActionList + _T("]");

        CCString sStepDetail = _T("");
        if (0 != i)
        {
            sStepDetail = appToString(static_cast<INT>(m_uiStepCount));
            for (INT j = 1; j <= i; ++j)
            {
                sStepDetail = sStepDetail + _T(" x ") + appToString(static_cast<INT>(m_uiNestedStep[j - 1]));
            }
            sStepDetail = _T("(") + sStepDetail + _T(")");
        }

        CCString sThisLine;
        sThisLine.Format(_T("Level:%d, Step:%d%s, Tau:%f, Actions:%s"),
            i, iStep, sStepDetail.c_str(), m_fNestedStepLengths[i], sActionList.c_str());
        sRet = sRet + sTab + sThisLine + _T("\n");
    }
    return sRet;
}

void CMultiLevelNestedIntegrator::UpdateP(Real fStep, TArray<UINT> actionList, ESolverPhase ePhase, UBOOL bUpdateP)
{
    _RECORD(CMultiLevelNestedIntegrator::UpdateP);

    TArray<CAction*> lst;
    UBOOL bOnlyFermion = TRUE;
    UBOOL bOnlyGauge = TRUE;
    for (INT i = 0; i < actionList.Num(); ++i)
    {
        //this is accumulate
        CAction* pAction = m_lstActions[actionList[i]];
        lst.AddItem(pAction);
        if (pAction->IsFermion())
        {
            bOnlyGauge = FALSE;
        }
        else
        {
            bOnlyFermion = FALSE;
        }
    }
    CalcForceOfActions(lst, bOnlyGauge ? EFC_Gauge : (bOnlyFermion ? EFC_Fermion : EFC_All), ePhase);

    //P = P + e F
    if (bUpdateP)
    {
        AddForce(fStep, FALSE);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CMultiLevelNestedIntegrator::NestedEvaluateLeapfrog(INT iLevel, Real fNestedStepLength, UBOOL bFirst, UBOOL bLast)
{
    //It is h/2m for nested, so times 0.5
    const UINT uiStepAll = (0 == iLevel) ? m_uiStepCount : m_uiNestedStep[iLevel - 1];
    fNestedStepLength = fNestedStepLength / uiStepAll;
    const Real fHalfEstep = F(0.5) * fNestedStepLength;
    //appDetailed(_T("  level(inner leapfrog) %d step 0, -------- length : %f \n"), iLevel, fNestedStepLength);

    UpdateP(fHalfEstep, iLevel,
        bFirst ? ESP_StartTrajectory : ESP_InTrajectory,
        TRUE);

    if (m_bDebugForce)
    {
        appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, CalcForce());
    }

    for (UINT uiStep = 1; uiStep < uiStepAll + 1; ++uiStep)
    {
        // middle step, exp(h/2 T) or NextNest
        if (iLevel != m_uiNestedStep.Num())
        {
            //NextNest
            NestedEvaluateLeapfrog(iLevel + 1, fNestedStepLength, bFirst, uiStep == uiStepAll);
        }
        else
        {
            //exp(h/2 T)
            //appDetailed("  nested(inner leapfrog) level %d step %d U\n", iLevel, uiStep);
            UpdateU(fNestedStepLength);
        }

        if (uiStep < uiStepAll)
        {
            appDetailed("  nested(inner leapfrog) level %d step %d P\n", iLevel, uiStep);
            UpdateP(fNestedStepLength, iLevel, ESP_InTrajectory, TRUE);
        }
        else
        {
            appDetailed("  nested(inner leapfrog) level %d last step %d\n", iLevel, uiStep);
            UpdateP(fHalfEstep, iLevel, bLast ? ESP_EndTrajectory : ESP_InTrajectory, TRUE);
        }

        if (m_bDebugForce)
        {
            appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, CalcForce());
        }
    }
}

__CLGIMPLEMENT_CLASS(CIntegratorLeapFrog)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================