//=============================================================================
// FILENAME : MatchingRho.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#include "MatchingRho.h"

int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("MatchingRho.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/MatchingRho.yaml"), params);
#endif

    appSetupLog(params);

    INT iVaule = 200;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 2500;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OnlyMeasure"), iVaule);
    const UBOOL bOnlyMeasure = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("DoMeasureFermion"), iVaule);
    const UBOOL bMeasureFermion = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("DoSmearing"), iVaule);
    const UBOOL bDoSmearing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SaveIndexStart"), iVaule);
    const UINT iSaveIndexStart = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    TArray<Real> pionCorrelator;
    TArray<Real> rhoCorrelator;
    //TArray<Real> potentialR;
    //TArray<CLGComplex> potentialC;

    //=========================================================
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    //CMeasurePolyakov* pPL = dynamic_cast<CMeasurePolyakov*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureMesonCorrelator* pMC = dynamic_cast<CMeasureMesonCorrelator*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CFieldGaugeSU3* pStaple = NULL;
    if (bDoSmearing)
    {
        pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField->GetCopy());
    }

    //Themalization
    appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
    //Initial of the gauge field depends on the .yaml file
    //appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    appGetLattice()->m_pMeasurements->Reset();
    UINT uiAccepCountBeforeE = 0;
    if (!bOnlyMeasure)
    {
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            const UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
            if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
            {
                uiAccepCountBeforeE = uiAccepCountBeforeE2;
            }
        }
    }

    //================================================================
    //Start working
    UINT uiAccepCountAfterE = 0;
    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    appGetLattice()->m_pMeasurements->Reset();
    while (
        (bOnlyMeasure && uiAccepCountAfterE < iEquib)
     || (!bOnlyMeasure && appGetLattice()->m_pUpdator->GetConfigurationCount() < iEquib)
        )
    {
        UINT uiAccepCountBeforeE2 = uiAccepCountAfterE;
        if (!bOnlyMeasure)
        {
            uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, TRUE);
            if (uiAccepCountAfterE != uiAccepCountBeforeE2)
            {
                uiAccepCountAfterE = uiAccepCountBeforeE2;

                //save measures
                for (UINT uiLt = 0; uiLt < _HC_Lt; ++uiLt)
                {
                    pionCorrelator.AddItem(pMC->m_lstResultsLastConf[0][uiLt]);
                    rhoCorrelator.AddItem((
                        pMC->m_lstResultsLastConf[1][uiLt]
                        + pMC->m_lstResultsLastConf[2][uiLt]
                        + pMC->m_lstResultsLastConf[3][uiLt]) / F(3.0));
                }

                //for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                //{
                //    if (0 == potentialR.Num())
                //    {
                //        potentialR.AddItem(_hostsqrt(pPL->m_lstR[i]));
                //    }

                //    potentialC.AddItem(pPL->m_lstC[(uiAccepCountAfterE - 1) * pPL->m_lstR.Num() + i]);
                //}

                //save configurations
                sFileName.Format(_T("Matching_%d"), uiAccepCountAfterE + iSaveIndexStart);
                sFileName = sSavePrefix + sFileName;
                //=================================
                //Save info
                appGetTimeNow(buff1, 256);
                appGetTimeUtc(buff2, 256);
                sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\n"),
                    appGetTimeStamp(),
                    buff1,
                    buff2);
                sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                //=================================
                //Save config
                appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));
            }
        }
        else
        {
            ++uiAccepCountAfterE;
            sFileName.Format(_T("Matching_%d.con"), uiAccepCountAfterE);
            sFileName = sSavePrefix + sFileName;
            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            if (bDoSmearing)
            {
                appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField, pStaple);
            }

            if (bMeasureFermion)
            {
                pMC->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                for (UINT uiLt = 0; uiLt < _HC_Lt; ++uiLt)
                {
                    pionCorrelator.AddItem(pMC->m_lstResultsLastConf[0][uiLt]);
                    rhoCorrelator.AddItem((
                        pMC->m_lstResultsLastConf[1][uiLt]
                      + pMC->m_lstResultsLastConf[2][uiLt]
                      + pMC->m_lstResultsLastConf[3][uiLt]) / F(3.0));
                }
            }
            else
            {
                pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            }

            appSetLogDate(FALSE);
            appGeneral(0 == uiAccepCountAfterE % 50 ? _T("\n=") : _T("="));
            appSetLogDate(TRUE);
        }
    }
    if (!bOnlyMeasure)
    {
        appGetLattice()->m_pMeasurements->Report();
    }
    else
    {
        if (!bMeasureFermion)
        {
            pPL->Report();

            appSetLogDate(FALSE);

            //extract result
            assert(static_cast<INT>(iEquib)* pPL->m_lstR.Num() == pPL->m_lstP.Num());
            UINT uiMaxL = (_HC_Lx + 1) / 2 - 1;
            uiMaxL = uiMaxL * uiMaxL;
            TArray<CCString> r_idx;

            appGeneral(_T("pr={"));

            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                appGeneral(_T("%2.12f%s"), _hostsqrt(static_cast<Real>(pPL->m_lstR[i])), (i == pPL->m_lstR.Num() - 1) ? _T("") : _T(", "));
            }

            appGeneral(_T("};\n"));

            appGeneral(_T("pl={\n"));

            for (UINT j = 0; j < iEquib; ++j)
            {
                appGeneral(_T("{"));
                for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                {
                    const CLGComplex cV = pPL->m_lstP[j * pPL->m_lstR.Num() + i];
                    appGeneral(_T("%2.12f %s %2.12f I%s"), cV.x, cV.y < F(0.0) ? _T("") : _T("+"), cV.y, (i == pPL->m_lstR.Num() - 1) ? _T("") : _T(", "));
                }
                appGeneral(_T("}%s\n"), (j == iEquib - 1) ? _T("") : _T(","));
            }

            appGeneral(_T("\n};\n"));

            appSetLogDate(TRUE);
        }
    }

    //=================================
    //report final result
    // we are satisfied with the report of Polyakov Loop, so only report the Meason correlator

    if (!bOnlyMeasure || bMeasureFermion)
    {
        appSetLogDate(FALSE);
        appGeneral(_T("\n ==================== Pion correlator C(p=0, nt) ==============\n\n"));
        appGeneral(_T("{\n"));
        for (UINT iConf = 0; iConf < uiAccepCountAfterE; ++iConf)
        {
            appGeneral(_T("{"));
            for (UINT iT = 0; iT < _HC_Lt; ++iT)
            {
                appGeneral(_T("%2.12f%s "), pionCorrelator[iConf * _HC_Lt + iT], (iT == _HC_Lt - 1) ? _T("") : _T(","));
            }
            appGeneral(_T("}%s\n"), (iConf == uiAccepCountAfterE - 1) ? _T("") : _T(","));
        }
        appGeneral(_T("}\n"));

        appGeneral(_T("\n ==================== Pho correlator C(p=0, nt) ==============\n\n"));
        appGeneral(_T("{\n"));
        for (UINT iConf = 0; iConf < uiAccepCountAfterE; ++iConf)
        {
            appGeneral(_T("{"));
            for (UINT iT = 0; iT < _HC_Lt; ++iT)
            {
                appGeneral(_T("%2.12f%s "), rhoCorrelator[iConf * _HC_Lt + iT], (iT == _HC_Lt - 1) ? _T("") : _T(","));
            }
            appGeneral(_T("}%s\n"), (iConf == uiAccepCountAfterE - 1) ? _T("") : _T(","));
        }
        appGeneral(_T("}\n"));

        appSetLogDate(TRUE);
    }

    if (bDoSmearing)
    {
        delete pStaple;
    }

    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
