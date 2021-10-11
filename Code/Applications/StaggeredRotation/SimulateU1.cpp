//=============================================================================
// FILENAME : SimulateU1.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [10/05/2021 nbale]
//=============================================================================

#include "StaggeredRotation.h"

INT SimulateStaggeredRotationU1(CParameters& params)
{
#pragma region Parameters

    appSetupLog(params);

    INT iVaule = 99;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("EquvibSkip"), iVaule);
    const UINT iEquibSkip = static_cast<UINT>(iVaule);

    iVaule = 250;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    const UINT iMinNt = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    const UINT iMaxNt = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    const UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OmegaStart"), iVaule);
    const UINT iOmegaStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    const UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_plaquttes;
    params.FetchValueArrayReal(_T("Plaquttes"), old_plaquttes);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    TArray<CCString> sOldFileNames;
    TArray<Real> fOldFilePolyakov;

    for (UINT i = iMinNt; i <= iMaxNt; ++i)
    {
        CCString sFileName;
        Real fPlaqutte = F(0.0);
        CCString sFileKeyName;
        CCString sPlaqKeyName;
        sFileKeyName.Format(_T("Nt%dFileName"), i);
        sPlaqKeyName.Format(_T("Nt%dPlaqutte"), i);
        params.FetchStringValue(sFileKeyName, sFileName);
        params.FetchValueReal(sPlaqKeyName, fPlaqutte);
        sOldFileNames.AddItem(sFileName);
        fOldFilePolyakov.AddItem(fPlaqutte);

        appGeneral(_T("file: %s, plaq : %f\n"), sFileName.c_str(), fPlaqutte);
    }

    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;
#pragma endregion

    for (UINT uiNt = iMinNt; uiNt <= iMaxNt; ++uiNt)
    {
        TArray<INT> latticeDecomp;
        params.FetchValueArrayINT(_T("LatticeLength"), latticeDecomp);
        latticeDecomp[3] = uiNt;
        TArray<CCString> sLatticeDecomp;
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[0]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[1]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[2]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[3]));
        params.SetStringVectorVaule(_T("LatticeLength"), sLatticeDecomp);

        if (!appInitialCLG(params))
        {
            appCrucial(_T("Initial Failed!\n"));
            return 1;
        }

        CMeasurePlaqutteEnergy* pPE = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
        TArray<Real> plaqutteList;

        CActionGaugePlaquetteRotatingU1* pGaugeRotation = dynamic_cast<CActionGaugePlaquetteRotatingU1*>(appGetLattice()->GetActionById(1));
        CCString sHeader;
        sHeader.Format(_T("Nt%d"), uiNt);
        appSetLogHeader(sHeader);
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);

        //=============== Check oldfiles ==================
        UBOOL bNeedBake = TRUE;
        if (!bAdditive && !sOldFileNames[uiNt - iMinNt].IsEmpty())
        {
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sOldFileNames[uiNt - iMinNt], EFFT_CLGBin);
            pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            const Real fError = appAbs(pPE->m_lstData[0] - fOldFilePolyakov[uiNt - iMinNt]);
#if _CLG_DOUBLEFLOAT
            if (fError < F(1E-07))
#else
            if (fError < F(1E-05))
#endif
            {
                appGeneral(_T("\n ================ Bake using old file =================\n"));
                bNeedBake = FALSE;
            }
            else
            {
                appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                    pPE->m_lstData[0], fOldFilePolyakov[uiNt - iMinNt], fError);
            }
        }

        if (bAdditive)
        {
            bNeedBake = FALSE;
        }

        if (bNeedBake && iBeforeEquib > 0)
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
            pGaugeRotation->SetOmega(F(0.0));

            appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
            appGetLattice()->m_pMeasurements->Reset();
            UINT uiAccepCountBeforeE = 0;
            while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
            {
                const UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
                if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
                {
                    uiAccepCountBeforeE = uiAccepCountBeforeE2;
                    pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
            }
            assert(pPE->m_lstData.Num() == static_cast<INT>(iBeforeEquib));
            appSetLogDate(FALSE);
            appGeneral(_T("\n Plaq ={\n"));
            for (INT i = 0; i < pPE->m_lstData.Num(); ++i)
            {
                appGeneral(_T("%f, "), pPE->m_lstData[i]);
            }
            appGeneral(_T("}\n"));
            appSetLogDate(TRUE);
        }
        else
        {
            appGeneral(_T("Not Baked\n"));
        }

        UINT uiOmega = iOmegaStart;
        const Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            sHeader.Format(_T("Nt%dO%d"), uiNt, uiOmega);
            appSetLogHeader(sHeader);
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);

            pGaugeRotation->SetOmega(fSep * uiOmega);

            if (bAdditive)
            {
                if (old_plaquttes.Num() <= static_cast<INT>(uiOmega))
                {
                    appGeneral(_T("\n ================ not have the initial value===========\n"));
                    appFailQuitCLG();
                    return 1;
                }
                const Real fPlaqOld = old_plaquttes[uiOmega];

                appGetLattice()->m_pMeasurements->Reset();
                sFileName.Format(_T("%sR_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, iSaveStartIndex);
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
                pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                const Real fError = appAbs(pPE->m_lstData[0] - fPlaqOld);
#if _CLG_DOUBLEFLOAT
                if (fError < F(1E-07))
#else
                if (fError < F(1E-05))
#endif
                {
                    appGeneral(_T("\n ================ using old file start from %d =================\n"), iSaveStartIndex);
                }
                else
                {
                    appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                        pPE->m_lstData[0], fPlaqOld, fError);
                    appFailQuitCLG();
                    return 1;
                }
            }

            UINT iConfigNumberNow = 0;
            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);

            while (iConfigNumberNow < iEquib)
            {
                appGetLattice()->m_pUpdator->Update(1, (iConfigNumberNow < iEquibSkip) ? FALSE : TRUE);
                const UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
                if (uiAcce != iConfigNumberNow)
                {
                    sFileName.Format(_T("R_Nt%d_O%d_%d"), uiNt, uiOmega, uiAcce + iSaveStartIndex);
                    sFileName = sSavePrefix + sFileName;

                    //=================================
                    //Save config
                    const CCString MD5 = appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));

                    //=================================
                    //Save info
                    appGetTimeNow(buff1, 256);
                    appGetTimeUtc(buff2, 256);
                    sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\nMD5 : %s\n"),
                        appGetTimeStamp(),
                        buff1,
                        buff2,
                        MD5.c_str());
                    sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                    appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                    iConfigNumberNow = uiAcce;
                }
            }

#pragma region gather measurements

            appGetLattice()->m_pMeasurements->Report();

            //===================== Polyakov loop =====================
            assert(pPE->m_lstData.Num() == static_cast<INT>(iEquib - iEquibSkip));

            //============= polyakov gather =============
            plaqutteList.AddItem(pPE->m_fLastRealResult);

#pragma endregion

            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), uiNt);
        appSetLogDate(FALSE);
        assert(plaqutteList.Num() == static_cast<INT>(iAfterEquib + 1));

        appGeneral(_T("Plaquttes={\n"));
        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n ") : _T("%2.10f,\n"),
                plaqutteList[i]);
        }
        appGeneral(_T("}\n\n"));

        appSetLogDate(TRUE);

        appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();

    }

    return 0;
}
