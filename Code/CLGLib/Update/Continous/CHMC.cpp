//=============================================================================
// FILENAME : CHMC.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [mm/dd/yy]
//  [12/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CHMC)


CHMC::~CHMC()
{
    appSafeDelete(m_pIntegrator);
}

void CHMC::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    m_iAcceptedConfigurationCount = 0;

    INT iMetro = 0;
    params.FetchValueINT(_T("Metropolis"), iMetro);
    m_bMetropolis = (0 != iMetro);

    INT iSave = 0;
    params.FetchValueINT(_T("SaveConfiguration"), iSave);
    m_bSaveConfigurations = (0 != iSave);

    INT iAdaptive = 0;
    params.FetchValueINT(_T("Adaptive"), iAdaptive);
    m_bAdaptiveUpdate = (0 != iAdaptive);

    INT iReport = 1;
    params.FetchValueINT(_T("ReportMeasure"), iReport);
    m_bReport = (0 != iReport);

    INT iSkip = 1;
    params.FetchValueINT(_T("Skip"), iSkip);
    m_uiSkip = static_cast<UINT>(iSkip);

    if (m_bSaveConfigurations)
    {
        m_sConfigurationPrefix = _T("Untitled");
        params.FetchStringValue(_T("ConfigurationFilePrefix"), m_sConfigurationPrefix);
        m_sConfigurationPrefix.Format(_T("%s_%d"), m_sConfigurationPrefix.c_str(), appGetTimeStamp());

        CCString sSaveType = _T("EFFT_CLGBin");
        if (params.FetchStringValue(_T("ConfigurationFileType"), sSaveType))
        {
            m_eSaveFieldType = __STRING_TO_ENUM(EFieldFileType, sSaveType);
        }
    }

    if (m_bAdaptiveUpdate)
    {
        TArray<INT> minMax;
        if (params.FetchValueArrayINT(_T("MinMaxStep"), minMax))
        {
            if (minMax.Num() > 1 && minMax[0] >= 1 && minMax[1] > minMax[0])
            {
                m_uiMinStep = static_cast<UINT>(minMax[0]);
                m_uiMaxStep = static_cast<UINT>(minMax[1]);
            }
        }

        TArray<Real> growReduceThreshold;
        if (params.FetchValueArrayReal(_T("GrowReduceThreshold"), growReduceThreshold))
        {
            if (minMax.Num() > 1 
             && growReduceThreshold[1] > F(0.000001) 
             && growReduceThreshold[0] < -F(0.000001) - growReduceThreshold[1])
            {
                m_fGrowStep = growReduceThreshold[0];
                m_fReduceStep = growReduceThreshold[1];
            }
        }
    }
}

UINT CHMC::Update(UINT iSteps, UBOOL bMeasure)
{
    _RECORD(CHMC::Update);
    ++m_uiUpdateCall;
    UINT uiCurrentStep = 0;
    UINT uiAccept = 0;
    UINT uiCurrentStepAccept = 0;
    UBOOL bAccepted = FALSE;
    CCString sAccept = appDressColor(EVC_GREEN, _T("Accept"));
    CCString sReject = appDressColor(EVC_RED, _T("Reject"));

    DOUBLE fEnergy = 0.0;
    DOUBLE fEnergyNew = 0.0;

    if (!m_bAdaptiveUpdate)
    {
        m_pIntegrator->FixStep(!m_bMetropolis);
    }

    TArray<DOUBLE> actions;
    for (UINT i = 0; i < iSteps; ++i)
    {
        if (0 == i)
        {
            m_pOwner->FixAllFieldBoundary();
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
        //m_pOwner->m_pGaugeField->DebugPrintMe();
        m_pIntegrator->Prepare(bAccepted, i);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
        m_pOwner->FixAllFieldBoundary();
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());

        actions.RemoveAll();
        if (m_bMetropolis || m_bTestHDiff)
        {
            fEnergy = m_pIntegrator->GetEnergy(TRUE, actions);
        }
        m_pIntegrator->Evaluate();
        if (m_bMetropolis || m_bTestHDiff)
        {
            fEnergyNew = m_pIntegrator->GetEnergy(FALSE, actions);
        }

        DOUBLE diff_H = 1.0;
        DOUBLE rand = 0.0;

        if (m_bMetropolis || m_bTestHDiff)
        {
            m_fLastHDiff = fEnergy - fEnergyNew;
            if (m_bTestHDiff)
            {
                m_lstHDiff.AddItem(m_fLastHDiff);
                m_lstH.AddItem(fEnergy);
            }

            if (m_bMetropolis)
            {
                diff_H = _hostexpd(m_fLastHDiff);  // Delta H (SA)
#if _CLG_MULTI_GPU
                if (NULL != appGetComm())
                {
                    //GetRandomReal() routes to CRandom::HostRandomF (mt19937 seeded
                    //by std::random_device in HostRandom.cpp) unless RandomType is
                    //ER_Schrage: the draw is NOT reproducible across processes, so
                    //a 1-vs-N accept-sequence comparison could never match even
                    //with the broadcast. The host Schrage stream is seeded by
                    //RandomSeed (rank-invariant) and nothing else consumes it, so
                    //drawing the Metropolis rand from it gives an identical
                    //sequence on every rank and every rank count. The broadcast
                    //stays as a guard so all ranks share the same verdict anyway.
                    rand = AMD * static_cast<DOUBLE>(appGetLattice()->m_pRandom->GetRandomUISchrage());
                    appGetComm()->BroadcastFromRoot(rand);
                }
                else
#endif
                {
                    rand = GetRandomReal();
                }
            }

            BYTE byUpdateChange = 0;
            if (m_bAdaptiveUpdate)
            {
                if (m_fLastHDiff < m_fGrowStep && m_pIntegrator->GetStepCount() < m_uiMaxStep)
                {
                    byUpdateChange = 1;
                    m_pIntegrator->ChangeStepCount(TRUE);
                }
                else if (appAbs(m_fLastHDiff) < m_fReduceStep && m_pIntegrator->GetStepCount() > m_uiMinStep)
                {
                    byUpdateChange = 2;
                    m_pIntegrator->ChangeStepCount(FALSE);
                }
            }

            //Metropolis
            appGeneral(_T(" HMC: step = %d, H_dff%s = %f (%f - %f)%s\n"),
                i + 1,
                m_bMetropolis ? _T("") : _T("(warmup)"),
                m_fLastHDiff,
                fEnergy,
                fEnergyNew,
                0 == byUpdateChange ? _T("") : (1 == byUpdateChange ? _T(", step++") : _T(", step--"))
                );
        }

        CCString sActionInfo;
        for (INT j = 0; j < actions.Num(); ++j)
        {
            if (j == 0)
            {
                sActionInfo += _T("\nkin old:") + appToString(actions[j]) + _T("\t");
            }
            else if (j == (actions.Num() / 2))
            {
                sActionInfo += _T("\nkin new:") + appToString(actions[j]) + _T("\t");
            }
            else if (j < (actions.Num() / 2))
            {
                sActionInfo += _T("act") + appToString(j) + _T(" old:") + appToString(actions[j]) + _T("\t");
            }
            else
            {
                sActionInfo += _T("act") + appToString(j - (actions.Num() / 2)) + _T(" new:") + appToString(actions[j]) + _T("\t");
            }
        }

        appGeneral(_T("%s\n"), sActionInfo.c_str());

        if (std::isnan(m_fLastHDiff) || is_nan_bitwise_robust(m_fLastHDiff) || std::isnan(fEnergy) || is_nan_bitwise_robust(fEnergy) || std::isnan(fEnergyNew) || is_nan_bitwise_robust(fEnergyNew))
        {
            //If we give up this trajectory, can we recover from this nan?
            appCrucial(_T("  Rejected because HDIff is nan (accepted:%d)\n"), uiAccept);
            bAccepted = FALSE;
        }
        else if (rand <= diff_H)
        {

            ++uiAccept;
            ++uiCurrentStepAccept;
            if (m_bMetropolis)
            {
                if (m_uiSkip > 1)
                {
                    appGeneral(_T(" random(0,1)=%f < exp(Hdff)=%f %s (total accepted:%d \ncurrent configuration: (accept/step/total)=%d/%d/%d )\n"), rand, diff_H, sAccept.c_str(), uiAccept,
                        uiCurrentStepAccept, uiCurrentStep + 1, m_uiSkip);
                }
                else
                {
                    appGeneral(_T(" random(0,1)=%f < exp(Hdff)=%f %s (accepted:%d)\n"), rand, diff_H, sAccept.c_str(), m_iAcceptedConfigurationCount + 1);
                }
            }
            else
            {
                appGeneral(_T(" Warmup %s (accepted:%d)\n"), sAccept.c_str(), uiAccept);
            }
            
            bAccepted = TRUE;
        }
        else
        {
            if (m_uiSkip > 1)
            {
                appGeneral(_T(" random(0,1)=%f > exp(Hdff)=%f %s (accepted:%d \ncurrent configuration: (accept/step/total)=%d/%d/%d)\n"), rand, diff_H, sReject.c_str(), uiAccept, 
                    uiCurrentStepAccept, uiCurrentStep + 1, m_uiSkip);
            }
            else
            {
                appGeneral(_T(" random(0,1)=%f > exp(Hdff)=%f %s (accepted:%d)\n"), rand, diff_H, sReject.c_str(), m_iAcceptedConfigurationCount);
            }
            
            bAccepted = FALSE;
        }
        m_pIntegrator->OnFinishTrajectory(bAccepted); //Here we copy the gauge field back
        checkCudaErrors(cudaGetLastError());

        //If rejected, just accept the old configuration and trigger the measure
        ++uiCurrentStep;

        if (m_bMetropolis)
        {
            if ((m_uiSkip < 2 && bAccepted)
             || (m_uiSkip >= 2 && uiCurrentStep >= m_uiSkip)
                )
            {
                uiCurrentStepAccept = 0;
                uiCurrentStep = 0;
                ++m_iAcceptedConfigurationCount;
                appGeneral(appDressColor(EVC_CYAN, _T("Get One Configuration")) + _T(" (configuration %d)\n"), m_iAcceptedConfigurationCount);

                if (bMeasure)
                {
                    //In 'OnFinishTrajectory', the field is already copy to 'CLatticeData'
                    TArray<const CFieldGauge*> gauges;
                    TArray<const CFieldBoson*> bosons;
                    TArray<const CFieldTensor2*> tensor2s;
                    for (INT j = 0; j < m_pIntegrator->m_pGaugeField.Num(); ++j)
                    {
                        gauges.AddItem(m_pIntegrator->m_pGaugeField[j]);
                    }
                    for (INT j = 0; j < m_pIntegrator->m_pBosonFields.Num(); ++j)
                    {
                        bosons.AddItem(m_pIntegrator->m_pBosonFields[j]);
                    }
                    for (INT j = 0; j < m_pIntegrator->m_pTensor2Field.Num(); ++j)
                    {
                        tensor2s.AddItem(m_pIntegrator->m_pTensor2Field[j]);
                    }

                    m_pOwner->FixAllFieldBoundary();
                    m_pOwner->OnUpdatorConfigurationAccepted(
                        gauges.Num(),
                        bosons.Num(),
                        tensor2s.Num(),
                        gauges.GetData(),
                        bosons.GetData(),
                        tensor2s.GetData(),
                        NULL);
                }

                if (m_bSaveConfigurations && !m_bTestHDiff)
                {
                    SaveConfiguration(i);
                }
            }
        }
        else
        {
            if (bMeasure)
            {
                //In 'OnFinishTrajectory', the field is already copy to 'CLatticeData'
                TArray<const CFieldGauge*> gauges;
                TArray<const CFieldBoson*> bosons;
                TArray<const CFieldTensor2*> tensor2s;
                for (INT j = 0; j < m_pIntegrator->m_pGaugeField.Num(); ++j)
                {
                    gauges.AddItem(m_pIntegrator->m_pGaugeField[j]);
                }
                for (INT j = 0; j < m_pIntegrator->m_pBosonFields.Num(); ++j)
                {
                    bosons.AddItem(m_pIntegrator->m_pBosonFields[j]);
                }
                for (INT j = 0; j < m_pIntegrator->m_pTensor2Field.Num(); ++j)
                {
                    tensor2s.AddItem(m_pIntegrator->m_pTensor2Field[j]);
                }

                m_pOwner->FixAllFieldBoundary();
                m_pOwner->OnUpdatorConfigurationAccepted(
                    gauges.Num(),
                    bosons.Num(),
                    tensor2s.Num(),
                    gauges.GetData(),
                    bosons.GetData(), 
                    tensor2s.GetData(),
                    NULL);
            }
        }
        appFlushLog();
    }

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    
    m_pOwner->OnUpdatorFinished(bMeasure, m_bReport);
#if !_CLG_DEBUG
    appFlushLog();
#endif
    return m_iAcceptedConfigurationCount;
}

void CHMC::UpdateUntileAccept(UINT iSteps, UBOOL bMeasure)
{
    if (m_uiSkip < 2)
    {
        m_iAcceptedConfigurationCount = 0;
        while (m_iAcceptedConfigurationCount < iSteps)
        {
            Update(1, bMeasure);
        }
    }
    else
    {
        m_iAcceptedConfigurationCount = 0;
        Update(iSteps * m_uiSkip, bMeasure);
    }
}

CCString CHMC::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : HMC\n");
    sRet = sRet + tab + _T("Integrator : \n");
    sRet = sRet + m_pIntegrator->GetInfos(tab + _T("    "));
    sRet = sRet + tab + _T("Metropolis : ") + (m_bMetropolis ? _T("1\n") : _T("0\n"));
    sRet = sRet + tab + _T("Skip : ") + appToString(m_uiSkip);
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================