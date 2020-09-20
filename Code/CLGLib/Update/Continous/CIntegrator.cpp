//=============================================================================
// FILENAME : CIntegrator.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CIntegrator::~CIntegrator()
{
    appSafeDelete(m_pGaugeField);
    appSafeDelete(m_pForceField);
    appSafeDelete(m_pMomentumField);
    appSafeDelete(m_pStapleField);
}

/**
* Create the fields here
*/
void CIntegrator::Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params)
{
    m_pOwner = pOwner;
    m_pLattice = pLattice;
    m_lstActions = pLattice->m_pActionList;
    m_bStapleCached = FALSE;

    INT iStepCount = 50;
    params.FetchValueINT(_T("IntegratorStep"), iStepCount);
    Real fStepLength = F(1.0);
    params.FetchValueReal(_T("IntegratorStepLength"), fStepLength);
    m_uiStepCount = (UINT)iStepCount;
    m_fEStep = fStepLength / m_uiStepCount;

    INT iDebugForce = 0;
    params.FetchValueINT(_T("DebugForce"), iDebugForce);
    m_bDebugForce = (0 != iDebugForce);

    INT iBindDir = 0;
    params.FetchValueINT(_T("BindDir"), iBindDir);
    m_byBindDir = static_cast<BYTE>(iBindDir);

    m_pGaugeField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pGaugeField->m_pOwner = pLattice;
    m_pGaugeField->InitialField(EFIT_Zero);

    m_pForceField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pForceField->m_pOwner = pLattice;
    m_pForceField->InitialField(EFIT_Zero);

    m_pMomentumField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
    m_pMomentumField->m_pOwner = pLattice;
    m_pMomentumField->InitialField(EFIT_Zero);

    if (CCommonData::m_bStoreStaple)
    {
        m_pStapleField = dynamic_cast<CFieldGauge*>(appCreate(pLattice->m_pGaugeField->GetClass()->GetName()));
        m_pStapleField->m_pOwner = pLattice;
        m_pStapleField->InitialField(EFIT_Zero);
    }
    else
    {
        m_pStapleField = NULL;
    }
}

void CIntegrator::Prepare(UBOOL bLastAccepted, UINT uiStep)
{
    //we may not accept the evaluation, so we need to copy it first
    if (!bLastAccepted || 0 == uiStep)
    {
        m_pLattice->m_pGaugeField->CopyTo(m_pGaugeField);
        m_pGaugeField->SetOneDirectionUnity(m_byBindDir);
        m_bStapleCached = FALSE;
        checkCudaErrors(cudaDeviceSynchronize());
    }

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->PrepareForHMC(m_pGaugeField, uiStep);
    }

    //generate a random momentum field to start
    m_pMomentumField->MakeRandomGenerator();
    m_pMomentumField->SetOneDirectionZero(m_byBindDir);
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_pGaugeField->CopyTo(m_pLattice->m_pGaugeField);
    }
    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        m_lstActions[i]->OnFinishTrajectory(bAccepted);
    }
    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::UpdateU(Real fStep) const
{
    m_pMomentumField->SetOneDirectionZero(m_byBindDir);

    //U(k) = exp (i e P) U(k-1)
    m_pMomentumField->ExpMult(fStep, m_pGaugeField);

    m_pGaugeField->SetOneDirectionUnity(m_byBindDir);

    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::UpdateP(Real fStep, UBOOL bCacheStaple, ESolverPhase ePhase)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    m_pGaugeField->SetOneDirectionUnity(m_byBindDir);

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, (0 == i && bCacheStaple) ? m_pStapleField : NULL, ePhase);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    m_pForceField->SetOneDirectionZero(m_byBindDir);
    
    //P = P + e F
    m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    m_pMomentumField->Axpy(fStep, m_pForceField);
    m_pMomentumField->SetOneDirectionZero(m_byBindDir);

    checkCudaErrors(cudaDeviceSynchronize());
}

void CIntegrator::FinishEvaluate() const
{
    m_pGaugeField->ElementNormalize();
    m_pGaugeField->SetOneDirectionUnity(m_byBindDir);
}

Real CIntegrator::GetEnergy(UBOOL bBeforeEvolution) const
{
    m_pMomentumField->SetOneDirectionZero(m_byBindDir);
    m_pGaugeField->SetOneDirectionUnity(m_byBindDir);

    Real retv = m_pMomentumField->CalculateKinematicEnergy();

#if _CLG_DEBUG
    CCString sLog = _T("");
    sLog.Format(_T("kin:%f, "), retv);
#endif

    for (INT i = 0; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
#if _CLG_DEBUG
        Real fActionEnergy = m_bStapleCached ? m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField);
        CCString sThisActionInfo = _T("");
        sThisActionInfo.Format(_T(" Action%d:%f, "), i + 1, fActionEnergy);
        sLog += sThisActionInfo;
        retv += fActionEnergy;
#else
        retv += m_bStapleCached ? m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField, m_pStapleField) : m_lstActions[i]->Energy(bBeforeEvolution, m_pGaugeField);
#endif
    }
#if _CLG_DEBUG
    appDetailed(_T("H (%s) = %s \n"), bBeforeEvolution ? "before" : "after" , sLog.c_str());
#endif
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
    return sTab + _T("Nested : ") + appIntToString(static_cast<INT>(m_uiNestedStep)) + _T("\n")
         + sTab + _T("InnerLeapFrog : ") + (m_bInnerLeapFrog ? _T("1") : _T("0")) + _T("\n");
}

void CNestedIntegrator::UpdatePF(Real fStep, ESolverPhase ePhase)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    for (INT i = 1; i < m_lstActions.Num(); ++i)
    {
        //this is accumulate
        m_lstActions[i]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL, ePhase);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //P = P + e F
    m_pMomentumField->Axpy(fStep, m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());

    if (m_bDebugForce)
    {
        const CLGComplex force = m_pForceField->Dot(m_pForceField);
        appGeneral(_T(" ------ Fermion Force= %f \n"), force.x);
    }
}

void CNestedIntegrator::UpdatePG(Real fStep, UBOOL bCacheStaple)
{
    // recalc force
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    m_lstActions[0]->CalculateForceOnGauge(m_pGaugeField, m_pForceField, bCacheStaple ? m_pStapleField : NULL, ESP_Once);
    checkCudaErrors(cudaDeviceSynchronize());

    //P = P + e F
    m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
    m_pMomentumField->Axpy(fStep, m_pForceField);
    checkCudaErrors(cudaDeviceSynchronize());

    if (m_bDebugForce)
    {
        const CLGComplex force = m_pForceField->Dot(m_pForceField);
        appGeneral(_T(" ------ Gauge Force= %f \n"), force.x);
    }
}

void CNestedIntegrator::NestedEvaluateLeapfrog(UBOOL bLast)
{
    const Real fHalfPstep = F(0.5) * m_fNestedStepLength;
    UpdatePG(fHalfPstep, FALSE);
    appDetailed("  leap frog nested sub step 0\n");

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(m_fNestedStepLength);

        if (uiStep < m_uiNestedStep)
        {
            UpdatePG(m_fNestedStepLength, FALSE);
            appDetailed("  leap frog nested sub step %d\n", uiStep);
        }
        else
        {
            UpdatePG(fHalfPstep, bLast);
            appDetailed("  leap frog nested last step %d\n", uiStep);
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
        params.FetchValueArrayUINT(_T("NestedActionList") + appIntToString(i), actionlist);
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
            sActionList = sActionList + appIntToString(static_cast<INT>(m_iNestedActionId[i][j]));
            if (j != m_iNestedActionId[i].Num() - 1)
            {
                sActionList = sActionList + _T(", ");
            }
        }
        sActionList = sActionList + _T("]");

        CCString sStepDetail = _T("");
        if (0 != i)
        {
            sStepDetail = appIntToString(static_cast<INT>(m_uiStepCount));
            for (INT j = 1; j <= i; ++j)
            {
                sStepDetail = sStepDetail + _T(" x ") + appIntToString(static_cast<INT>(m_uiNestedStep[j - 1]));
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

void CMultiLevelNestedIntegrator::UpdateP(Real fStep, TArray<UINT> actionList, ESolverPhase ePhase, UBOOL bCacheStaple, UBOOL bUpdateP)
{
    m_pForceField->Zero();
    checkCudaErrors(cudaDeviceSynchronize());

    for (INT i = 0; i < actionList.Num(); ++i)
    {
        //this is accumulate
        const CAction* pAction = m_lstActions[actionList[i]];
        if (pAction->IsFermion())
        {
            pAction->CalculateForceOnGauge(m_pGaugeField, m_pForceField, NULL, ePhase);
        }
        else
        {
            pAction->CalculateForceOnGauge(m_pGaugeField, m_pForceField, bCacheStaple ? m_pStapleField : NULL, ESP_Once);
            m_bStapleCached = CCommonData::m_bStoreStaple && bCacheStaple;
        }
        
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //P = P + e F
    if (bUpdateP)
    {
        m_pMomentumField->Axpy(fStep, m_pForceField);
    }
    checkCudaErrors(cudaDeviceSynchronize());
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
        FALSE,
        TRUE);

    if (m_bDebugForce)
    {
        const CLGComplex force = m_pForceField->Dot(m_pForceField);
        appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, force.x);
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
            UpdateP(fNestedStepLength, iLevel, ESP_InTrajectory, FALSE, TRUE);
        }
        else
        {
            appDetailed("  nested(inner leapfrog) level %d last step %d\n", iLevel, uiStep);
            UpdateP(fHalfEstep, iLevel, bLast ? ESP_EndTrajectory : ESP_InTrajectory, bLast, TRUE);
        }

        if (m_bDebugForce)
        {
            const CLGComplex force = m_pForceField->Dot(m_pForceField);
            appGeneral(_T(" ------ Force (%d) = %f \n"), iLevel, force.x);
        }
    }
}

__CLGIMPLEMENT_CLASS(CIntegratorLeapFrog)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================