//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"

INT SimulateStaggeredEM(CParameters& params)
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

    iVaule = 0;
    params.FetchValueINT(_T("EMStart"), iVaule);
    UINT iEMStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("EMEnd"), iVaule);
    UINT iEMEnd = static_cast<UINT>(iVaule);

    TArray<Real> lstE;
    params.FetchValueArrayReal(_T("ElectricList"), lstE);

    TArray<Real> lstM;
    params.FetchValueArrayReal(_T("MagneticList"), lstM);

    if (0 == lstE.Num() && 0 != lstM.Num())
    {
        for (INT i = 0; i < lstM.Num(); ++i)
        {
            lstE.AddItem(F(0.0));
        }
    }

    if (0 == lstM.Num() && 0 != lstE.Num())
    {
        for (INT i = 0; i < lstE.Num(); ++i)
        {
            lstM.AddItem(F(0.0));
        }
    }

    if ((0 == lstE.Num() && 0 == lstM.Num()) || lstE.Num() != lstM.Num())
    {
        appCrucial(_T("The length of the list ElectricList or MagneticList is wrong!\n"));
        return 1;
    }

    if (iEMEnd < static_cast<UINT>(lstE.Num()))
    {
        iEMEnd = static_cast<UINT>(lstE.Num());
    }

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_polyakovs;
    params.FetchValueArrayReal(_T("Polyakovs"), old_polyakovs);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

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

    CCString sOldFileName;
    Real fOldPolyakov = F(0.0);
    CCString sOldFileKeyName;
    CCString sOldPolyaKeyName;
    sOldFileKeyName.Format(_T("Nt%dFileName"), _HC_Lt);
    sOldPolyaKeyName.Format(_T("Nt%dPolyakov"), _HC_Lt);
    params.FetchStringValue(sOldFileKeyName, sOldFileName);
    params.FetchValueReal(sOldPolyaKeyName, fOldPolyakov);

    appGeneral(_T("file: %s, |p| : %f\n"), sOldFileName.c_str(), fOldPolyakov);

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    //TArray<TArray<CLGComplex>> polykovX_nx;
    TArray<CLGComplex> polykov;
    TArray<Real> polykovphase;
    //for (UINT uiX = 0; uiX < static_cast<UINT>(CCommonData::m_sCenter.x); ++uiX)
    //{
    //    TArray<CLGComplex> a;
    //    TArray<Real> b;
    //    polykovX_nx.AddItem(a);
    //}

    //CActionGaugePlaquetteRotating* pGaugeRotation = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    //const BYTE uiFermionFieldCount = static_cast<BYTE>(appGetLattice()->m_pOtherFields.Num());

    CCString sHeader;
    sHeader.Format(_T("Nt%d"), _HC_Lt);
    appSetLogHeader(sHeader);
    appGeneral(_T("Run for Nt = %d, start baking.\n"), _HC_Lt);

    //=============== Check oldfiles ==================
    UBOOL bNeedBake = TRUE;
    if (!bAdditive && !sOldFileName.IsEmpty())
    {
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sOldFileName, EFFT_CLGBin);
        pPL->OnConfigurationAccepted(_FIELDS, NULL);
        Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldPolyakov);
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
                _cuCabsf(pPL->m_lstLoop[0]), fOldPolyakov, fError);
        }
    }

    if (bAdditive)
    {
        bNeedBake = FALSE;
    }

    if (bNeedBake && iBeforeEquib > 0)
    {
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        CCommonData::m_fBz = F(0.0);
        CCommonData::m_fEz = F(0.0);

        appGetLattice()->m_pGaugeField[0]->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pMeasurements->Reset();
        UINT uiAccepCountBeforeE = 0;
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
            if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
            {
                uiAccepCountBeforeE = uiAccepCountBeforeE2;
                pPL->OnConfigurationAccepted(_FIELDS, NULL);
            }
        }
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iBeforeEquib));
        appPushLogDate(FALSE);
        appGeneral(_T("\n|<P>|,arg<P>={\n"));
        for (INT i = 0; i < pPL->m_lstLoop.Num(); ++i)
        {
            appGeneral(_T("{%f, %f},\n"), _cuCabsf(pPL->m_lstLoop[i]), __cuCargf(pPL->m_lstLoop[i]));
        }
        appGeneral(_T("}\n"));
        appPopLogDate();
    }
    else
    {
        appGeneral(_T("Not Baked\n"));
    }

    UINT uiEM = iEMStart;
    while (uiEM < iEMEnd)
    {
        sHeader.Format(_T("Nt%dEM%d"), _HC_Lt, uiEM);
        appSetLogHeader(sHeader);
        appGeneral(_T("\n========= E=%f M=%f ==========\n"), lstE[uiEM], lstM[uiEM]);
        CCommonData::m_fBz = lstM[uiEM];
        CCommonData::m_fEz = lstE[uiEM];

        if (bAdditive)
        {
            Real fPolyaOld = F(0.0);
            if (old_polyakovs.Num() <= static_cast<INT>(uiEM))
            {
                appGeneral(_T("\n ================ not have the initial value===========\n"));
                appFailQuitCLG();
                return 1;
            }
            fPolyaOld = old_polyakovs[uiEM];

            appGetLattice()->m_pMeasurements->Reset();
            sFileName.Format(_T("%sR_Nt%d_EM%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiEM, iSaveStartIndex);
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pPL->OnConfigurationAccepted(_FIELDS, NULL);
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
            UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
            if (uiAcce != iConfigNumberNow)
            {
                sFileName.Format(_T("R_Nt%d_EM%d_%d"), _HC_Lt, uiEM, uiAcce + iSaveStartIndex);
                sFileName = sSavePrefix + sFileName;

                //=================================
                //Save config
                const CCString MD5 = appGetLattice()->m_pGaugeField[0]->SaveToFile(sFileName + _T(".con"));

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
        //assert(pPL->m_lstAverageLoopDensity.Num()
        //    == static_cast<INT>(CCommonData::m_sCenter.x));

        //============= polyakov gather =============
        polykov.AddItem(pPL->m_cAverageLoop);
        polykovphase.AddItem(__cuCargf(pPL->m_cAverageLoop));
        //for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
        //{
        //    polykovX_nx[iX].AddItem(pPL->m_lstAverageLoopDensity[iX]);
        //}

#pragma endregion

        ++uiEM;

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    }

    appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), _HC_Lt);
    appPushLogDate(FALSE);
    assert(polykov.Num() == static_cast<INT>(iEMEnd - iEMStart));

    appGeneral(_T("|Polyakov|={\n"));
    for (INT i = 0; i < static_cast<INT>(iEMEnd - iEMStart); ++i)
    {
        appGeneral((i == static_cast<INT>(iEMEnd - iEMStart - 1)) ? _T("%2.10f\n ") : _T("%2.10f,\n"),
            _cuCabsf(polykov[i]));
    }
    appGeneral(_T("}\n\narg(Polyakov)={\n"));

    for (INT i = 0; i < static_cast<INT>(iEMEnd - iEMStart); ++i)
    {
        appGeneral((i == static_cast<INT>(iEMEnd - iEMStart - 1)) ? _T("%2.10f\n ") : _T("%2.10f,\n"),
            polykovphase[i]);
    }

    //for (UINT x = 0; x < static_cast<UINT>(CCommonData::m_sCenter.x); ++x)
    //{
    //    appGeneral(_T("Polyakov[x=%d]={\n"), x);
    //    for (INT i = 0; i < static_cast<INT>(iEMEnd - iEMStart); ++i)
    //    {
    //        appGeneral((i == static_cast<INT>(iEMEnd - iEMStart - 1)) ? _T("%2.10f %s %2.10f I\n") : _T("%2.10f %s %2.10f I,\n"),
    //            polykovX_nx[x][i].x,
    //            polykovX_nx[x][i].y < F(0.0) ? _T("-") : _T("+"),
    //            appAbs(polykovX_nx[x][i].y)
    //        );
    //    }
    //    appGeneral(_T("}\n\n"));
    //}

    appPopLogDate();

    appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), _HC_Lt);
    appQuitCLG();

    return 0;
}
