//=============================================================================
// FILENAME : CHMC.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
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

    if (m_bSaveConfigurations)
    {
        m_sConfigurationPrefix = _T("Untitled");
        params.FetchStringValue(_T("ConfigurationFilePrefix"), m_sConfigurationPrefix);
        m_sConfigurationPrefix.Format(_T("%s_%d"), m_sConfigurationPrefix.c_str(), appGetTimeStamp());
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
    ++m_uiUpdateCall;
    UBOOL bAccepted = FALSE;

#if !_CLG_DOUBLEFLOAT
    DOUBLE fEnergy = 0.0;
    DOUBLE fEnergyNew = 0.0;
#else
    Real fEnergy = F(0.0);
    Real fEnergyNew = F(0.0);
#endif

    for (UINT i = 0; i < iSteps; ++i)
    {
        m_pIntegrator->Prepare(bAccepted, i);
        m_pOwner->FixAllFieldBoundary();
        if (m_bMetropolis || m_bTestHDiff)
        {
            fEnergy = m_pIntegrator->GetEnergy(TRUE);
        }
        m_pIntegrator->Evaluate();
        if (m_bMetropolis || m_bTestHDiff)
        {
            fEnergyNew = m_pIntegrator->GetEnergy(FALSE);
        }

#if !_CLG_DOUBLEFLOAT
        DOUBLE diff_H = 1.0;
        DOUBLE rand = 0.0;
#else
        Real diff_H = F(1.0);
        Real rand = F(0.0);
#endif

        if (m_bMetropolis || m_bTestHDiff)
        {
#if !_CLG_DOUBLEFLOAT
            m_fLastHDiff = fEnergy - fEnergyNew;
#else
            m_fLastHDiff = fEnergy - fEnergyNew;
#endif
            if (m_bTestHDiff)
            {
                m_lstHDiff.AddItem(m_fLastHDiff);
                m_lstH.AddItem(fEnergy);
            }

            if (!m_bTestHDiff)
            {
#if !_CLG_DOUBLEFLOAT
                diff_H = _hostexpd(m_fLastHDiff);  // Delta H (SA)
#else
                diff_H = _hostexp(m_fLastHDiff);  // Delta H (SA)
#endif
                rand = GetRandomReal();
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
            appGeneral(_T(" HMC: step = %d, H_dff = %f (%f - %f)%s\n"),
                i + 1,
                m_fLastHDiff,
                fEnergy,
                fEnergyNew,
                0 == byUpdateChange ? _T("") : (1 == byUpdateChange ? _T(", step++") : _T(", step--"))
                );
        }

        if (rand <= diff_H)
        {
            ++m_iAcceptedConfigurationCount;
            appGeneral(_T("  Accepted (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            bAccepted = TRUE;
        }
        else
        {
            appGeneral(_T("  Rejected (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            bAccepted = FALSE;
        }
        m_pIntegrator->OnFinishTrajectory(bAccepted); //Here we copy the gauge field back

        //If rejected, just accept the old configuration and trigger the measure
        if (bMeasure && bAccepted)
        {
            m_pOwner->FixAllFieldBoundary();
            m_pOwner->OnUpdatorConfigurationAccepted(
                appGetLattice()->m_pGaugeField,
                bAccepted ? m_pIntegrator->m_pStapleField : NULL);
        }

        if (m_bSaveConfigurations && (bAccepted || 0 == i))
        {
            SaveConfiguration(i + 1);
        }
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

CCString CHMC::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : HMC\n");
    sRet = sRet + tab + _T("Integrator : \n");
    sRet = sRet + m_pIntegrator->GetInfos(tab + _T("    "));
    sRet = sRet + tab + _T("Metropolis : ") + (m_bMetropolis ? _T("1\n") : _T("0\n"));
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================