//=============================================================================
// FILENAME : RotationImprovedFermion.cpp
// 
// DESCRIPTION:
//
// REVISION:
// [mm/dd/yy]
// [08/22/2025 nbale]
//=============================================================================

#include "RotationImprovedFermion.h"

extern INT SimulateRotation(CParameters& params)
{
    appSetupLog(params);

    INT iVaule = 100;
    params.FetchValueINT(_T("Warmup"), iVaule);
    UINT iWarmUp = static_cast<UINT>(iVaule);

    iVaule = 1000;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 250;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    UINT iMinNt = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OmegaStart"), iVaule);
    UINT iOmegaStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    TArray<CCString> sOldFileNames;

    for (UINT i = iMinNt; i <= iMaxNt; ++i)
    {
        CCString sFileName;
        CCString sFileKeyName;
        sFileKeyName.Format(_T("Nt%dFileName"), i);
        params.FetchStringValue(sFileKeyName, sFileName);
        sOldFileNames.AddItem(sFileName);
        appGeneral(_T("file: %s \n"), sFileName.c_str());
    }

    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    CCString sFileName;

    for (UINT uiNt = iMinNt; uiNt <= iMaxNt; ++uiNt)
    {
        TArray<INT> latticeDecomp;
        params.FetchValueArrayINT(_T("LatticeLength"), latticeDecomp);
        latticeDecomp[3] = uiNt;
        TArray<CCString> sLatticeDecomp;
        sLatticeDecomp.AddItem(appToString(latticeDecomp[0]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[1]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[2]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[3]));
        params.SetStringVectorVaule(_T("LatticeLength"), sLatticeDecomp);

        if (!appInitialCLG(params))
        {
            appCrucial(_T("Initial Failed!\n"));
            return 1;
        }

        CActionGaugePlaquetteRotating* pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
        
        CCString sHeader;
        sHeader.Format(_T("Nt%d"), uiNt);
        appSetLogHeader(sHeader);
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);

        //=============== Check oldfiles ==================
        UBOOL bNeedBake = TRUE;
        if (!bAdditive && !sOldFileNames[uiNt - iMinNt].IsEmpty())
        {
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sOldFileNames[uiNt - iMinNt], EFFT_CLGBin);
            appGeneral(_T("\n ================ Bake using old file: %s =================\n"), sOldFileNames[uiNt - iMinNt].c_str());
            bNeedBake = FALSE;
        }

        if (bAdditive)
        {
            bNeedBake = FALSE;
        }

        if (bNeedBake && iWarmUp > 0)
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));

            appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
            appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
            appGetLattice()->m_pUpdator->Update(iWarmUp, FALSE);
        }
        else
        {
            appGeneral(_T("Not Baked\n"));
        }

        UINT uiOmega = iOmegaStart;
        Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            sHeader.Format(_T("Nt%dO%d"), uiNt, uiOmega);
            appSetLogHeader(sHeader);
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);

            if (NULL != pGauageAction)
            {
                pGauageAction->SetGaugeOmega(fSep * uiOmega);
                appGeneral(_T("\n========= Gauge Action set omega =%f  ==========\n"), fSep * uiOmega);               
            }

            INT iFermionFieldNum = 0;
            for (BYTE i = 2; i < 10; ++i)
            {
                //CFieldFermionWilsonSquareCloverSU3DR is inherent from CFieldFermionWilsonSquareSU3DR
                CFieldFermionWilsonSquareSU3DR* pFermion1 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(i));
                if (NULL != pFermion1)
                {
                    pFermion1->SetFermionOmega(fSep * uiOmega);
                    ++iFermionFieldNum;
                }
                //CFieldFermionHISQSU3R : public CFieldFermionHISQT<CFieldFermionKSSU3R>
                CFieldFermionKSSU3R* pFermion2 = dynamic_cast<CFieldFermionKSSU3R*>(appGetLattice()->GetFieldById(i));
                if (NULL != pFermion2)
                {
                    pFermion2->SetFermionOmega(fSep * uiOmega);
                    ++iFermionFieldNum;
                }
            }
            appGeneral(_T("\n========= %d fermion fields set omega =%f  ==========\n"), iFermionFieldNum, fSep * uiOmega);

            if (bAdditive)
            {
                appGetLattice()->m_pMeasurements->Reset();
                sFileName.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, iSaveStartIndex - 1);
                appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
                appGeneral(_T("\n ================ start from file =%s ===========\n"), sFileName.c_str());
            }

            appGetLattice()->m_pMeasurements->Reset();

            sFileName.Format(_T("%sRotate_Nt%d_O%d"), sSavePrefix.c_str(), uiNt, uiOmega);
            appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sFileName, static_cast<UINT>(iSaveStartIndex));
            appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
            //TODO SET FILE INDEX START
            appGetLattice()->m_pUpdator->UpdateUntileAccept(iEquib, FALSE);

            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), uiNt);
        appQuitCLG();
    }

    return 0;
}


//=============================================================================
// END OF FILE
//=============================================================================
