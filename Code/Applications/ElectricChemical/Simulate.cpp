//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/29/2022 nbale]
//=============================================================================
#include "ElectricChemical.h"

INT Simulate(CParameters& params)
{
#pragma region Parameters

    appSetupLog(params);

    INT iVaule = 99;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("EquvibSkip"), iVaule);
    UINT iEquibSkip = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix : %s\n"), sSavePrefix.c_str());

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const INT iListStart = iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const INT iListEnd = iVaule;

    TArray<Real> lstElectric;
    params.FetchValueArrayReal(_T("Electric"), lstElectric);

    TArray<Real> lstChemical;
    params.FetchValueArrayReal(_T("Chemical"), lstChemical);

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_polyakovs;
    params.FetchValueArrayReal(_T("Polyakovs"), old_polyakovs);

    CCString sOldFileName;
    const UBOOL bHasOldFile = params.FetchStringValue(_T("OldFileName"), sOldFileName);

    Real fOldPolyakov = F(0.0);
    if (bHasOldFile)
    {
        params.FetchValueReal(_T("OldPolyakov"), fOldPolyakov);
    }

    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;

#pragma endregion

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    TArray<TArray<CLGComplex>> polykovX_nx;
    TArray<CLGComplex> polykov;
    TArray<Real> polykovphase;
    for (UINT uiX = 0; uiX < static_cast<UINT>(CCommonData::m_sCenter.x); ++uiX)
    {
        TArray<CLGComplex> a;
        TArray<Real> b;
        polykovX_nx.AddItem(a);
    }

    CFieldFermionKSSU3GammaEM* pU = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetFieldById(2));
    CFieldFermionKSSU3GammaEM* pD = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetFieldById(3));
    CFieldGaugeU1Real* pU1 = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(4));

    UBOOL bNeedBake = TRUE;
    if (!bAdditive && bHasOldFile)
    {
        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sOldFileName, EFFT_CLGBin);
        pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
        Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldPolyakov);
#if _CLG_DOUBLEFLOAT
        if (fError < F(1E-07))
#else
        if (fError < F(1E-05))
#endif
        {
            appGeneral(_T("\n ================ Bake using old file %s =================\n"), sOldFileName.c_str());
        }
        else
        {
            appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                _cuCabsf(pPL->m_lstLoop[0]), fOldPolyakov, fError);
        }
        bNeedBake = FALSE;
    }

    if (bAdditive)
    {
        bNeedBake = FALSE;
    }

    if (bNeedBake && iBeforeEquib > 0)
    {
        if (bHasOldFile)
        {
            appGeneral(_T("!!!! Has old file but still bake !!!!\n"));
        }
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pMeasurements->Reset();
        UINT uiAccepCountBeforeE = 0;
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
            if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
            {
                uiAccepCountBeforeE = uiAccepCountBeforeE2;
                pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            }
        }
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iBeforeEquib));
        appSetLogDate(FALSE);
        appGeneral(_T("\n|<P>|,arg<P>={\n"));
        for (INT i = 0; i < pPL->m_lstLoop.Num(); ++i)
        {
            appGeneral(_T("{%f, %f},\n"), _cuCabsf(pPL->m_lstLoop[i]), __cuCargf(pPL->m_lstLoop[i]));
        }
        appGeneral(_T("}\n"));
        appSetLogDate(TRUE);
    }
    else
    {
        appGeneral(_T("Not Baked\n"));
    }

    for (INT uiOmega = iListStart; uiOmega < lstElectric.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        CCString sHeader;
        sHeader.Format(_T("%d"), uiOmega);
        appSetLogHeader(sHeader);
        appGeneral(_T("\n========= Electric =%f Chemical = %f  ==========\n"), lstElectric[uiOmega], lstChemical[uiOmega]);

        pU->m_fCoeffGamma54 = lstChemical[uiOmega];
        pU->UpdatePooledParamters();
        pD->m_fCoeffGamma54 = lstChemical[uiOmega];
        pD->UpdatePooledParamters();
        pU1->InitialU1Real(EURT_None, EURT_E_t, EURT_None, F(0.0), lstElectric[uiOmega], F(0.0));

        //pU1->DebugPrintMe();

        if (bAdditive)
        {
            Real fPolyaOld = F(0.0);
            if (old_polyakovs.Num() <= static_cast<INT>(uiOmega))
            {
                appGeneral(_T("\n ================ not have the initial value===========\n"));
                appFailQuitCLG();
                return 1;
            }
            fPolyaOld = old_polyakovs[uiOmega];

            appGetLattice()->m_pMeasurements->Reset();
            sFileName.Format(_T("%sEC_%d_%d.con"), sSavePrefix.c_str(), uiOmega, iSaveStartIndex);
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fPolyaOld);
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
                    _cuCabsf(pPL->m_lstLoop[0]), fPolyaOld, fError);
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
                sFileName.Format(_T("%sEC_%d_%d"), sSavePrefix.c_str(), uiOmega, uiAcce + iSaveStartIndex);

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
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iEquib - iEquibSkip));
        assert(pPL->m_lstAverageLoopDensity.Num()
            == static_cast<INT>(CCommonData::m_sCenter.x));

        //============= polyakov gather =============
        polykov.AddItem(pPL->m_cAverageLoop);
        polykovphase.AddItem(__cuCargf(pPL->m_cAverageLoop));
        for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
        {
            polykovX_nx[iX].AddItem(pPL->m_lstAverageLoopDensity[iX]);
        }

#pragma endregion

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    }

    appGeneral(_T("\n========= Finished! ==========\n\n"));

    appQuitCLG(); 

    return 0;
}
