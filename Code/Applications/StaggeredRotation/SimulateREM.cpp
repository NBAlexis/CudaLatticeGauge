//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"

INT SimulateStaggeredRotationEM(CParameters& params)
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
    UINT iEquibSkip = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const UINT iEMStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const UINT iEMEnd = static_cast<UINT>(iVaule);

    TArray<Real> lstMagnetic;
    params.FetchValueArrayReal(_T("MagneticList"), lstMagnetic);

    TArray<Real> lstOmega;
    params.FetchValueArrayReal(_T("OmegaList"), lstOmega);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    const UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    const UBOOL bAdditive = 0 != iVaule;

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

    //TArray<INT> latticeDecomp;
    //params.FetchValueArrayINT(_T("LatticeLength"), latticeDecomp);
    //latticeDecomp[3] = uiNt;
    //TArray<CCString> sLatticeDecomp;
    //sLatticeDecomp.AddItem(appIntToString(latticeDecomp[0]));
    //sLatticeDecomp.AddItem(appIntToString(latticeDecomp[1]));
    //sLatticeDecomp.AddItem(appIntToString(latticeDecomp[2]));
    //sLatticeDecomp.AddItem(appIntToString(latticeDecomp[3]));
    //params.SetStringVectorVaule(_T("LatticeLength"), sLatticeDecomp);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CCString sOldFileName;
    Real fOldFilePolyakov = F(0.0);
    CCString sFileKeyName;
    CCString sPolyaKeyName;
    sFileKeyName.Format(_T("Nt%dFileName"), _HC_Lt);
    sPolyaKeyName.Format(_T("Nt%dPolyakov"), _HC_Lt);
    params.FetchStringValue(sFileKeyName, sOldFileName);
    params.FetchValueReal(sPolyaKeyName, fOldFilePolyakov);

    appGeneral(_T("file: %s, |p| : %f\n"), sOldFileName.c_str(), fOldFilePolyakov);

    CFieldGaugeU1Real* pU1 = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(6));
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

    CActionGaugePlaquetteRotating* pGaugeRotation = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    CFieldFermionKSSU3REM* pF2 = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(2));
    CFieldFermionKSSU3REM* pF3 = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(3));
    CFieldFermionKSSU3REM* pF4 = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(4));
    CFieldFermionKSSU3REM* pF5 = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetFieldById(5));
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
        Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldFilePolyakov);
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
                _cuCabsf(pPL->m_lstLoop[0]), fOldFilePolyakov, fError);
        }
    }

    if (bAdditive)
    {
        bNeedBake = FALSE;
    }

    if (bNeedBake && iBeforeEquib > 0)
    {
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        pGaugeRotation->SetGaugeOmega(F(0.0));
        pF2->SetFermionOmega(F(0.0));
        pF3->SetFermionOmega(F(0.0));
        pF4->SetFermionOmega(F(0.0));
        pF5->SetFermionOmega(F(0.0));
        pU1->InitialU1Real(EURT_None, EURT_None, pU1->m_eB, F(0.0), F(0.0), F(0.0), TRUE);

        appGetLattice()->m_pGaugeField[0]->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pMeasurements->Reset();
        UINT uiAccepCountBeforeE = 0;
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            const UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
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

    UINT uiListIdx = iEMStart;
    while (uiListIdx < iEMEnd)
    {
        sHeader.Format(_T("Nt%dREM%d"), _HC_Lt, uiListIdx);
        appSetLogHeader(sHeader);
        appGeneral(_T("\n========= Omega=%f Magnetic=%f ==========\n"), lstOmega[uiListIdx], lstMagnetic[uiListIdx]);

        pGaugeRotation->SetGaugeOmega(lstOmega[uiListIdx]);
        pF2->SetFermionOmega(lstOmega[uiListIdx]);
        pF3->SetFermionOmega(lstOmega[uiListIdx]);
        pF4->SetFermionOmega(lstOmega[uiListIdx]);
        pF5->SetFermionOmega(lstOmega[uiListIdx]);
        pU1->InitialU1Real(EURT_None, EURT_None, pU1->m_eB, F(0.0), F(0.0), lstMagnetic[uiListIdx], TRUE);

        if (bAdditive)
        {
            if (old_polyakovs.Num() <= static_cast<INT>(uiListIdx))
            {
                appGeneral(_T("\n ================ not have the initial value===========\n"));
                appFailQuitCLG();
                return 1;
            }
            const Real fPolyaOld = old_polyakovs[uiListIdx];

            appGetLattice()->m_pMeasurements->Reset();
            sFileName.Format(_T("%sR_Nt%d_REM%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiListIdx, iSaveStartIndex);
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            const Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fPolyaOld);
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
                sFileName.Format(_T("R_Nt%d_REM%d_%d"), _HC_Lt, uiListIdx, uiAcce + iSaveStartIndex);
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

        ++uiListIdx;

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    }

    appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), _HC_Lt);
    appPushLogDate(FALSE);
    const INT uiRealRuned = static_cast<INT>(iEMEnd - iEMStart);
    assert(polykov.Num() == uiRealRuned);

    appGeneral(_T("|Polyakov|={\n"));
    for (INT i = 0; i <= uiRealRuned; ++i)
    {
        appGeneral(i == uiRealRuned ? _T("%2.10f\n ") : _T("%2.10f,\n"),
            _cuCabsf(polykov[i]));
    }
    appGeneral(_T("}\n\narg(Polyakov)={\n"));

    for (INT i = 0; i <= uiRealRuned; ++i)
    {
        appGeneral(i == uiRealRuned ? _T("%2.10f\n ") : _T("%2.10f,\n"),
            polykovphase[i]);
    }

    //for (INT x = 0; x < static_cast<INT>(CCommonData::m_sCenter.x); ++x)
    //{
    //    appGeneral(_T("Polyakov[x=%d]={\n"), x);
    //    for (INT i = 0; i <= uiRealRuned; ++i)
    //    {
    //        appGeneral(i == uiRealRuned ? _T("%2.10f %s %2.10f I\n") : _T("%2.10f %s %2.10f I,\n"),
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
