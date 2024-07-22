//=============================================================================
// FILENAME : MatchingRho.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#include "MatchingRho.h"

#define __Show_Correlator(name) \
appGeneral(_T("\n ==================== %s correlator C(p=0, nt) ==============\n\n"), _T(#name)); \
appGeneral(_T("{\n")); \
for (UINT iConf = 0; iConf < uiAccepCountAfterE; ++iConf) \
{ \
    appGeneral(_T("{")); \
    for (UINT iT = 0; iT < _HC_Lt; ++iT) \
    { \
        appGeneral(_T("%2.20f%s "), name##Correlator[iConf * _HC_Lt + iT], (iT == _HC_Lt - 1) ? _T("") : _T(",")); \
    } \
    appGeneral(_T("}%s\n"), (iConf == uiAccepCountAfterE - 1) ? _T("") : _T(",")); \
} \
appGeneral(_T("}\n")); \


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
    params.FetchValueINT(_T("DoMeasureCondensation"), iVaule);
    const UBOOL bMeasureCondensation = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("DoMeasureTrace"), iVaule);
    const UBOOL bMeasureTrace = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    const UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("CompressedFile"), iVaule);
    const UBOOL bCompressedFile = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("DoSmearing"), iVaule);
    const UBOOL bDoSmearing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("TestZero"), iVaule);
    const UBOOL bTestZero = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SaveIndexStart"), iVaule);
    const UINT iSaveIndexStart = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVPrefix;
    params.FetchStringValue(_T("CSVPrefix"), sCSVPrefix);
    appGeneral(_T("CSV prefix: %s\n"), sCSVPrefix.c_str());

#if !_CLG_DOUBLEFLOAT
    TArray<DOUBLE> pionCorrelator;
    TArray<DOUBLE> rhoCorrelator;

    TArray<DOUBLE> rho0Correlator;
    TArray<DOUBLE> rho1Correlator;
    TArray<DOUBLE> rho2Correlator;
    TArray<DOUBLE> rho3Correlator;
#else
    TArray<Real> pionCorrelator;
    TArray<Real> rhoCorrelator;

    TArray<Real> rho0Correlator;
    TArray<Real> rho1Correlator;
    TArray<Real> rho2Correlator;
    TArray<Real> rho3Correlator;
#endif
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
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    CMeasurePandChiralTalor* pTalor = dynamic_cast<CMeasurePandChiralTalor*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    pTalor->m_fBetaOverN = appGetLattice()->m_pActionList[0]->GetBetaOverN();

    CFieldGaugeSU3* pStaple = NULL;
    if (bDoSmearing)
    {
        pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
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
            //If do some simulation
            uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, TRUE);
            if (uiAccepCountAfterE != uiAccepCountBeforeE2)
            {
                uiAccepCountAfterE = uiAccepCountBeforeE2;

                if (NULL != pMC)
                {
                    //save measures
                    for (UINT uiLt = 0; uiLt < _HC_Lt; ++uiLt)
                    {
                        pionCorrelator.AddItem(pMC->m_lstResultsLastConf[0][uiLt]);
                        rho1Correlator.AddItem(pMC->m_lstResultsLastConf[1][uiLt]);
                        rho2Correlator.AddItem(pMC->m_lstResultsLastConf[2][uiLt]);
                        rho3Correlator.AddItem(pMC->m_lstResultsLastConf[3][uiLt]);

#if !_CLG_DOUBLEFLOAT
                        rhoCorrelator.AddItem((
                            pMC->m_lstResultsLastConf[1][uiLt]
                            + pMC->m_lstResultsLastConf[2][uiLt]
                            + pMC->m_lstResultsLastConf[3][uiLt]) / 3.0);
#else
                        rhoCorrelator.AddItem((
                            pMC->m_lstResultsLastConf[1][uiLt]
                            + pMC->m_lstResultsLastConf[2][uiLt]
                            + pMC->m_lstResultsLastConf[3][uiLt]) / F(3.0));
#endif
                        rho0Correlator.AddItem(pMC->m_lstResultsLastConf[4][uiLt]);
                    }
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
                //Save config
                CCString sMD5;
                if (bCompressedFile)
                {
                    sMD5 = appGetLattice()->m_pGaugeField[0]->SaveToCompressedFile(sFileName + _T(".cco"));
                }
                else
                {
                    sMD5 = appGetLattice()->m_pGaugeField[0]->SaveToFile(sFileName + _T(".con"));
                }

                //=================================
                //Save info
                appGetTimeNow(buff1, 256);
                appGetTimeUtc(buff2, 256);
                sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\nMD5 : %s"),
                    appGetTimeStamp(),
                    buff1,
                    buff2,
                    sMD5.c_str());
                sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

            }
        }
        else
        {
            //only load configuration and measure
            ++uiAccepCountAfterE;
            sFileName.Format(_T("Matching_%d"), uiAccepCountAfterE + iSaveIndexStart);
            sFileName = sSavePrefix + sFileName;

            if (bCompressedFile)
            {
                appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName + _T(".cco"), EFFT_CLGBinCompressed);
            }
            else
            {
                appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName + _T(".con"), EFFT_CLGBin);
            }
            
            if (bDoSmearing)
            {
                appGetLattice()->m_pGaugeField[0]->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField[0], pStaple);
            }

            if (bMeasureFermion)
            {
                pMC->OnConfigurationAccepted(_FIELDS, NULL);
                for (UINT uiLt = 0; uiLt < _HC_Lt; ++uiLt)
                {
                    pionCorrelator.AddItem(pMC->m_lstResultsLastConf[0][uiLt]);
                    rho1Correlator.AddItem(pMC->m_lstResultsLastConf[1][uiLt]);
                    rho2Correlator.AddItem(pMC->m_lstResultsLastConf[2][uiLt]);
                    rho3Correlator.AddItem(pMC->m_lstResultsLastConf[3][uiLt]);
#if !_CLG_DOUBLEFLOAT
                    rhoCorrelator.AddItem((
                        pMC->m_lstResultsLastConf[1][uiLt]
                        + pMC->m_lstResultsLastConf[2][uiLt]
                        + pMC->m_lstResultsLastConf[3][uiLt]) / 3.0);
#else
                    rhoCorrelator.AddItem((
                        pMC->m_lstResultsLastConf[1][uiLt]
                      + pMC->m_lstResultsLastConf[2][uiLt]
                      + pMC->m_lstResultsLastConf[3][uiLt]) / F(3.0));
#endif
                    rho0Correlator.AddItem(pMC->m_lstResultsLastConf[4][uiLt]);
                }
            }
            else if (bMeasureCondensation)
            {
                if (NULL != pCC)
                {
                    if (bTestZero)
                    {
                        appGetLattice()->m_pGaugeField[0]->InitialField(EFIT_Identity);
                    }

                    CFieldFermionWilsonSquareSU3*  pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
                    CFieldFermionWilsonSquareSU3*  pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
                    const UINT iFieldCount = pCC->GetFieldCount();
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1->InitialField(EFIT_RandomGaussian);
                        }
                        pF1->FixBoundary();
                        pF1->CopyTo(pF2);
                        pF1->InverseD(_FIELDS);
                        pF1->FixBoundary();

                        pCC->OnConfigurationAcceptedZ4(
                            _FIELDS,
                            NULL,
                            pF2,
                            pF1,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                    pF1->Return();
                    pF2->Return();
                }
                else
                {
                    appGeneral(_T("Measurement is not found!\n"));
                }
            }
            else
            {
                if (NULL != pPL && !bMeasureTrace)
                {
                    pPL->OnConfigurationAccepted(_FIELDS, NULL);
                }
                
                if (bMeasureTrace && NULL != pTalor)
                {
                    pTalor->OnConfigurationAccepted(_FIELDS, NULL);
                    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
                    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
                    const UINT iFieldCount = pTalor->GetFieldCount();
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1->InitialField(EFIT_RandomGaussian);
                        }
                        //pF1->FixBoundary();
                        pF1->CopyTo(pF2);
                        pF1->InverseD(_FIELDS);
                        //pF1->FixBoundary();

                        pTalor->OnConfigurationAcceptedZ4(
                            _FIELDS,
                            NULL,
                            pF2,
                            pF1,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                    pF1->Return();
                    pF2->Return();
                }
            }

            if (0 == uiAccepCountAfterE % 100)
            {
                appGeneral(_T("\n %d ="), uiAccepCountAfterE);
            }
            else
            {
                appPushLogDate(FALSE);
                appGeneral(0 == uiAccepCountAfterE % 50 ? _T("\n=") : _T("="));
                appPopLogDate();
            }
        }
    }
    if (!bOnlyMeasure)
    {
        appGetLattice()->m_pMeasurements->Report();
    }
    else
    {
        if (!bMeasureFermion && !bMeasureTrace)
        {
            pPL->Report();

            appPushLogDate(FALSE);

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

            appPopLogDate();
        }
    }

    //=================================
    //report final result
    // we are satisfied with the report of Polyakov Loop, so only report the Meason correlator

    if ((!bOnlyMeasure || bMeasureFermion) && NULL != pMC)
    {
        appPushLogDate(FALSE);

        __Show_Correlator(pion);
        __Show_Correlator(rho);
        __Show_Correlator(rho1);
        __Show_Correlator(rho2);
        __Show_Correlator(rho3);
        __Show_Correlator(rho0);

        appPopLogDate();
    }

    if (bOnlyMeasure && bMeasureCondensation && NULL != pCC)
    {
        appPushLogDate(FALSE);

        appGeneral(_T("\n ==================== condensation ==============\n\n")); 
        appGeneral(_T("{\n")); 
        for (UINT iConf = 0; iConf < uiAccepCountAfterE; ++iConf) 
        { 
            appGeneral(_T("{"));
            for (UINT iType = 0; iType < CMeasureChiralCondensate::_kCondMeasureCount; ++iType)
            {
                appGeneral(_T("%2.18f %s %2.18f I%s"), 
                    pCC->m_lstCond[iType][iConf].x,
                    (pCC->m_lstCond[iType][iConf].y < 0) ? _T("-") : _T("+"),
                    appAbs(pCC->m_lstCond[iType][iConf].y),
                    ((iType + 1) == CMeasureChiralCondensate::_kCondMeasureCount) ? _T("") : _T(",")
                    );

            }
            appGeneral(_T("}%s\n"), (iConf == uiAccepCountAfterE - 1) ? _T("") : _T(",")); 
        }
        appGeneral(_T("}\n"));

        for (INT i = 0; i < 6; ++i)
        {
            TArray<CLGComplex> thisConfiguration; 
            for (UINT j = 0; j < uiAccepCountAfterE; ++j)
            { 
                thisConfiguration.AddItem(pCC->m_lstCond[0][j * pCC->m_lstR.Num() + i]); 
            }
            CCString sCSVFileName;
            sCSVFileName.Format(_T("%s_Chiral_r%d.csv"), sCSVPrefix.c_str(), i);
            WriteStringFileComplexArray(sCSVFileName, thisConfiguration);
        } 

        appPopLogDate();

        for (UINT iType = 0; iType < 1/*CMeasureChiralCondensate::_kCondMeasureCount*/; ++iType)
        {
            CCString sCSVFileName;
            sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), __ENUM_TO_STRING(ECondList, (ECondList)iType).c_str());
            WriteStringFileComplexArray(sCSVFileName, pCC->m_lstCondAll[iType]);
        }
    }

    if (bOnlyMeasure && bMeasureTrace && NULL != pTalor)
    {
        CCString sCSVFileName;
        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("D"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstTraceRes[ECPCTTT_D]);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("MD"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstTraceRes[ECPCTTT_MD]);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("DMD"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstTraceRes[ECPCTTT_DMD]);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("MDMD"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstTraceRes[ECPCTTT_MDMD]);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("DMDMD"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstTraceRes[ECPCTTT_DMDMD]);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("PL"));
        WriteStringFileComplexArray(sCSVFileName, pTalor->m_lstPolyakov);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("GO"));
        WriteStringFileRealArray(sCSVFileName, pTalor->m_lstPolyakovSOmega);

        sCSVFileName.Format(_T("%s_%s.csv"), sCSVPrefix.c_str(), _T("GOS"));
        WriteStringFileRealArray(sCSVFileName, pTalor->m_lstPolyakovSOmegaSq);
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
