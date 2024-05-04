//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"

__DEFINE_ENUM(EDistributionJobKSREM,
    EDJKSR_Polyakov,
    EDJKSR_Chiral,
    EDJKSR_AngularMomentum,
    EDJKSR_ChiralAndFermionMomentum,
    EDJKSR_PlaqutteEnergy,
    EDJKSR_CheckMD5,
    EDJKSR_VR,
    EDJKSR_DoubleToFloat,
    )



#define __REM_MEASURE_FERMION(ftype, savetag) \
    if (0 == (savetag & uiLoadFermion)) \
    { \
        if (bZ4) \
        { \
            pF1##ftype->InitialField(EFIT_RandomZ4); \
        } \
        else \
        { \
            pF1##ftype->InitialField(EFIT_RandomGaussian); \
        } \
        pF1##ftype->FixBoundary(); \
        pF1##ftype->CopyTo(pF2##ftype); \
        pF1##ftype->InverseD(appGetLattice()->m_pGaugeField); \
        pF1##ftype->FixBoundary(); \
        if (bSaveFermion) \
        { \
            CCString sFermionFile = ""; \
            sFermionFile.Format(_T("%s_%s_Nt%d_%s%d_%d_F%d"), sFermionHead.c_str(), #ftype, _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, uiSaveFermionStart + i); \
            CCString sMD51 = pF1##ftype->SaveToFile(sFermionFile + _T("_F1.con")); \
            CCString sMD52 = pF2##ftype->SaveToFile(sFermionFile + _T("_F2.con")); \
            CCString sFileContent = ""; \
            sFileContent = _T("Stochastic Fermion File for ") + sFileName; \
            if (bZ4) \
            { \
                sFileContent = sFileContent + _T("\nZ4\n"); \
            } \
            else \
            { \
                sFileContent = sFileContent + _T("\nGaussian\n"); \
            } \
            sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n"); \
            sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n"); \
            sFileContent = sFileContent + _T("Beta: ") + appToString(CCommonData::m_fBeta) + _T("\n"); \
            sFileContent = sFileContent + _T("Omega: ") + appToString(CCommonData::m_fOmega) + _T("\n"); \
            sFileContent = sFileContent + _T("Magnetic: ") + appToString(pU1->m_feBz) + _T("\n"); \
            sFileContent = sFileContent + _T("MagneticType: ") + __ENUM_TO_STRING(EU1RealType, pU1->m_eB) + _T("\n"); \
            sFileContent = sFileContent + _T("Mass: ") + appToString(pF1##ftype->m_f2am) + _T("\n"); \
            sFileContent = sFileContent + _T("Chage: ") + appToString(pF1##ftype->m_fQ) + _T("\n"); \
            sFileContent = sFileContent + _T("ShiftCenter: ") + (pF1##ftype->m_bEachSiteEta ? _T("TRUE") : _T("FALSE")) + _T("\n"); \
            appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent); \
        } \
    } \
    else \
    { \
        CCString sF1FileName = ""; \
        CCString sF2FileName = ""; \
        sF1FileName.Format(_T("%s%d/%s/%s_%s_Nt%d_%s%d_%d_F%d_F1.con"), \
            sLoadFermionHead.c_str(), uiListIdx, #ftype, sLoadFermionFile.c_str(), #ftype, _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, i + 1); \
        sF2FileName.Format(_T("%s%d/%s/%s_%s_Nt%d_%s%d_%d_F%d_F2.con"), \
            sLoadFermionHead.c_str(), uiListIdx, #ftype, sLoadFermionFile.c_str(), #ftype, _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, i + 1); \
        pF1##ftype->InitialFieldWithFile(sF1FileName, EFFT_CLGBin); \
        pF2##ftype->InitialFieldWithFile(sF2FileName, EFFT_CLGBin); \
    } \
 \
pCC##ftype->OnConfigurationAcceptedZ4( \
    1, 0, gauge.GetData(), NULL, \
    NULL, \
    pF2##ftype, \
    pF1##ftype, \
    0 == i, \
    iFieldCount == i + 1); \
\
pFA##ftype->OnConfigurationAcceptedZ4( \
    1, 0, gauge.GetData(), NULL, \
    NULL, \
    pF2##ftype, \
    pF1##ftype, \
    0 == i, \
    iFieldCount == i + 1); \


INT MeasurementREM(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const UINT iListStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const UINT iListEnd = static_cast<UINT>(iVaule);

    TArray<Real> lstMagnetic;
    params.FetchValueArrayReal(_T("MagneticList"), lstMagnetic);

    TArray<Real> lstOmega;
    params.FetchValueArrayReal(_T("OmegaList"), lstOmega);

    iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("MeasureCCS"), iVaule);
    //UBOOL bMeasureCCS = (0 != iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sValue = _T("EDJKSR_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJobKSREM eJob = __STRING_TO_ENUM(EDistributionJobKSREM, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sMiddle;
    params.FetchStringValue(_T("Middle"), sMiddle);
    appGeneral(_T("middle: %s\n"), sMiddle.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    Real fBeta = F(0.0);
    params.FetchValueReal(_T("GaugeBate"), fBeta);

    iVaule = 0;
    params.FetchValueINT(_T("SaveFermionFile"), iVaule);
    UBOOL bSaveFermion = 0 != iVaule;

    iVaule = 1;
    params.FetchValueINT(_T("FermionFileIndexStart"), iVaule);
    UINT uiSaveFermionStart = static_cast<UINT>(iVaule);

    CCString sFermionHead;
    params.FetchStringValue(_T("FermionFileHead"), sFermionHead);
    appGeneral(_T("FermionFileHead: %s\n"), sFermionHead.c_str());

    iVaule = 0;
    params.FetchValueINT(_T("LoadFermion"), iVaule);
    const UINT uiLoadFermion = iVaule;

    CCString sLoadFermionFile;
    params.FetchStringValue(_T("LoadFermionFile"), sLoadFermionFile);
    appGeneral(_T("Load Fermion File Name: %s\n"), sLoadFermionFile.c_str());

    CCString sLoadFermionHead;
    params.FetchStringValue(_T("LoadFermionHead"), sLoadFermionHead);
    appGeneral(_T("Load Fermion File Head: %s\n"), sLoadFermionHead.c_str());

    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }
    appGeneral(_T("load type: %s\n"), __ENUM_TO_STRING(EFieldFileType, eLoadType).c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

#pragma endregion

    UINT uiMaxL = (_HC_Lx + 1) / 2 - 1;
    uiMaxL = uiMaxL * uiMaxL;
    TArray<TArray<CCString>> r_omega_idx;
    for (UINT i = 0; i < iListEnd * iListEnd * uiMaxL; ++i)
    {
        TArray<CCString> newlst;
        r_omega_idx.AddItem(newlst);
    }
    TArray<Real> lstR;
    TArray<TArray<CLGComplex>> lstPolyIn;
    TArray<TArray<CLGComplex>> lstPolyOut;
    TArray<TArray<CLGComplex>> lstPolyInZ;
    TArray<TArray<CLGComplex>> lstPolyOutZ;

   
    CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CFieldGaugeU1Real* pU1 = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(5));
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCCu = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureChiralCondensateKS* pCCd = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureChiralCondensateKS* pCCs = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));

    CMeasureAngularMomentumKSREM* pFAu = dynamic_cast<CMeasureAngularMomentumKSREM*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    CMeasureAngularMomentumKSREM* pFAd = dynamic_cast<CMeasureAngularMomentumKSREM*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    CMeasureAngularMomentumKSREM* pFAs = dynamic_cast<CMeasureAngularMomentumKSREM*>(appGetLattice()->m_pMeasurements->GetMeasureById(7));

    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(8));

    CActionGaugePlaquetteRotating* pAG = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);

    CFieldFermionKSSU3REM* pF1u = NULL;
    CFieldFermionKSSU3REM* pF2u = NULL;
    CFieldFermionKSSU3REM* pF1d = NULL;
    CFieldFermionKSSU3REM* pF2d = NULL;
    CFieldFermionKSSU3REM* pF1s = NULL;
    CFieldFermionKSSU3REM* pF2s = NULL;

    if (EDJKSR_ChiralAndFermionMomentum == eJob
        || EDJKSR_Chiral == eJob)
    {
        pF1u = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(2));
        pF2u = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(2));
        pF1d = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(3));
        pF2d = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(3));
        pF1s = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(4));
        pF2s = dynamic_cast<CFieldFermionKSSU3REM*>(appGetLattice()->GetPooledFieldById(4));
    }

    appPushLogDate(FALSE);

    for (UINT uiListIdx = iListStart; uiListIdx < iListEnd; ++uiListIdx)
    {
        CCommonData::m_fOmega = lstOmega[uiListIdx];
        pU1->InitialU1Real(EURT_None, EURT_None, pU1->m_eB, F(0.0), F(0.0), lstMagnetic[uiListIdx]);

        if (NULL != pAG)
        {
            pAG->SetOmega(CCommonData::m_fOmega);
        }
        appGeneral(_T("(* ==== Omega(%f) Magnetic(%f) ========= *)\n"), lstOmega[uiListIdx], lstMagnetic[uiListIdx]);
        pPL->Reset();
        pJG->Reset();
        pCCu->Reset();
        pCCd->Reset();
        pCCs->Reset();
        pFAu->Reset();
        pFAd->Reset();
        pFAs->Reset();

        pCCu->SetFieldCount(iFieldCount);
        pCCd->SetFieldCount(iFieldCount);
        pCCs->SetFieldCount(iFieldCount);
        pFAu->SetFieldCount(iFieldCount);
        pFAd->SetFieldCount(iFieldCount);
        pFAs->SetFieldCount(iFieldCount);

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s%d/%sR_Nt%d_%s%d_%d.con"), sSubFolderPrefix.c_str(), uiListIdx, sSavePrefix.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN);
                sTxtFileName.Format(_T("%s%d/%sR_Nt%d_%s%d_%d.txt"), sSubFolderPrefix.c_str(), uiListIdx, sSavePrefix.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN);
            }
            else
            {
                sFileName.Format(_T("%sR_Nt%d_%s%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN);
                sTxtFileName.Format(_T("%sR_Nt%d_%s%d_%d.txt"), sSavePrefix.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN);
            }
            //appGeneral(_T("checking %s ..."), sFileName);
            if (EDJKSR_CheckMD5 == eJob)
            {
                UINT uiSize = 0;
                BYTE* fileContent = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
                CCString sMD5 = "MD5 : " + CLGMD5Hash(fileContent, uiSize);
                CCString sMD5old = "MD5 : " + CLGMD5Hash_OLD(fileContent, uiSize);
                CCString sFileContent = appGetFileSystem()->ReadAllText(sTxtFileName);
                if (sFileContent.Find(sMD5) >= 0)
                {
                    if (EDJKSR_CheckMD5 == eJob)
                    {
                        appGeneral(_T("MD5 Found and good %s \n"), sFileName.c_str());
                    }
                    else
                    {
                        appGeneral(_T("-"));
                    }
                }
                else if (sFileContent.Find(sMD5old) >= 0)
                {
                    if (EDJKSR_CheckMD5 == eJob)
                    {
                        appGeneral(_T("MD5 Found and good but use old bad MD5 %s \n"), sFileName.c_str());
                    }
                    else
                    {
                        appGeneral(_T("-"));
                    }
                    sFileContent = sFileContent.Replace(sMD5old, sMD5);
                    appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
                }
                else if (sFileContent.Find("MD5 : ") >= 0)
                {
                    appCrucial(_T("MD5 Found and NOT good %s \n"), sFileName.c_str());
                }
                else
                {
                    if (EDJKSR_CheckMD5 == eJob)
                    {
                        appGeneral(_T("MD5 Not Found so add to it %s \n"), sFileName.c_str());
                    }
                    else
                    {
                        appGeneral(_T("+"));
                    }
                    sFileContent = sFileContent + "\n" + sMD5 + "\n";
                    appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
                }

                if (EDJKSR_CheckMD5 == eJob)
                {
                    continue;
                }
            }
            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, eLoadType);
            TArray<CFieldGauge*> gauge;
            gauge.AddItem(appGetLattice()->m_pGaugeField);

            switch (eJob)
            {
                case EDJKSR_Polyakov:
                {
                    pPL->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                }
                break;
                case EDJKSR_Chiral:
                {
                    appCrucial(_T("Do not calculate chiral solely.\n"));
                }
                break;
                case EDJKSR_AngularMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                }
                break;
                case EDJKSR_ChiralAndFermionMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        //if (0 == (1 & uiLoadFermion))
                        //{
                        //    if (bZ4)
                        //    {
                        //        pF1u->InitialField(EFIT_RandomZ4);
                        //    }
                        //    else
                        //    {
                        //        pF1u->InitialField(EFIT_RandomGaussian);
                        //    }
                        //    pF1u->FixBoundary();
                        //    pF1u->CopyTo(pF2u);
                        //    pF1u->InverseD(appGetLattice()->m_pGaugeField);
                        //    pF1u->FixBoundary();
                        //    if (bSaveFermion)
                        //    {
                        //        CCString sFermionFile = "";
                        //        sFermionFile.Format(_T("%s_u_Nt%d_%s%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, uiSaveFermionStart + i);
                        //        CCString sMD51 = pF1u->SaveToFile(sFermionFile + _T("_F1.con"));
                        //        CCString sMD52 = pF2u->SaveToFile(sFermionFile + _T("_F2.con"));
                        //        CCString sFileContent = "";
                        //        sFileContent = _T("Stochastic Fermion File for ") + sFileName;
                        //        if (bZ4)
                        //        {
                        //            sFileContent = sFileContent + _T("\nZ4\n");
                        //        }
                        //        else
                        //        {
                        //            sFileContent = sFileContent + _T("\nGaussian\n");
                        //        }
                        //        sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n");
                        //        sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n");
                        //        sFileContent = sFileContent + _T("Beta: ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
                        //        sFileContent = sFileContent + _T("Omega: ") + appFloatToString(CCommonData::m_fOmega) + _T("\n");
                        //        sFileContent = sFileContent + _T("Magnetic: ") + appFloatToString(CCommonData::m_fBz) + _T("\n");
                        //        sFileContent = sFileContent + _T("Mass: ") + appFloatToString(pF1u->m_f) + _T("\n");
                        //        sFileContent = sFileContent + _T("ShiftCenter: ") + (pF1u->m_bEachSiteEta ? _T("TRUE") : _T("FALSE")) + _T("\n");
                        //        appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                        //    }
                        //}
                        //else
                        //{
                        //    CCString sF1FileName = "";
                        //    CCString sF2FileName = "";
                        //    sF1FileName.Format(_T("%s%d/u/%s_Light_Nt%d_%s%d_%d_F%d_F1.con"),
                        //        sLoadFermionHead.c_str(), uiListIdx, sLoadFermionFile.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, i + 1);
                        //    sF2FileName.Format(_T("%s%d/u/%s_Light_Nt%d_%s%d_%d_F%d_F2.con"),
                        //        sLoadFermionHead.c_str(), uiListIdx, sLoadFermionFile.c_str(), _HC_Lt, sMiddle.c_str(), uiListIdx, uiN, i + 1);
                        //    pF1u->InitialFieldWithFile(sF1FileName, EFFT_CLGBin);
                        //    pF2u->InitialFieldWithFile(sF2FileName, EFFT_CLGBin);
                        //}

                        //pCCu->OnConfigurationAcceptedZ4(
                        //    appGetLattice()->m_pGaugeField,
                        //    NULL,
                        //    pF2u,
                        //    pF1u,
                        //    0 == i,
                        //    iFieldCount == i + 1);

                        //pFAu->OnConfigurationAcceptedZ4(
                        //    appGetLattice()->m_pGaugeField,
                        //    NULL,
                        //    pF2u,
                        //    pF1u,
                        //    0 == i,
                        //    iFieldCount == i + 1);
                        __REM_MEASURE_FERMION(u, 1);
                        __REM_MEASURE_FERMION(d, 2);
                        __REM_MEASURE_FERMION(s, 4);

                    }
                }
                break;
                default:
                    break;
            }

            if ((iEndN - uiN + 1) % uiNewLine == 0)
            {
                appPushLogDate(TRUE);
                appGeneral(_T("\n="));
                appPopLogDate();
            }
            else
            {
                appGeneral(_T("="));
            }
            
        }
        appGeneral(_T("\n*)\n"));

        if (EDJKSR_CheckMD5 == eJob)
        {
            continue;
        }

#pragma endregion

        switch (eJob)
        {
            case EDJKSR_Polyakov:
            {
                CCString sFileNameWrite1;
                CCString sFileNameWrite2;
                sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_R.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_REM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiListIdx);
                
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pPL->m_lstR.Num() == pPL->m_lstP.Num());
                
                if (uiListIdx == iListStart)
                {
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        lstR.AddItem(F(0.5)* _hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                    }
                    WriteStringFileRealArray(sFileNameWrite1, lstR);
                }

                TArray<TArray<CLGComplex>> polyakovOmgR;
                TArray<CLGComplex> polyIn;
                TArray<CLGComplex> polyOut;

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    TArray<CLGComplex> thisConfiguration;
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        thisConfiguration.AddItem(pPL->m_lstP[j * pPL->m_lstR.Num() + i]);
                    }
                    polyakovOmgR.AddItem(thisConfiguration);
                    polyIn.AddItem(pPL->m_lstLoopInner[j]);
                    polyOut.AddItem(pPL->m_lstLoop[j]);
                }
                lstPolyIn.AddItem(polyIn);
                lstPolyOut.AddItem(polyOut);
                WriteStringFileComplexArray2(sFileNameWrite2, polyakovOmgR);

                if (pPL->m_bMeasureLoopZ)
                {
                    CCString sFileNameWrite3;
                    sFileNameWrite3.Format(_T("%s_polyakovZ_Nt%d_REM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiListIdx);
                    polyakovOmgR.RemoveAll();
                    polyIn.RemoveAll();
                    polyOut.RemoveAll();

                    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                    {
                        TArray<CLGComplex> thisConfiguration;
                        for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                        {
                            thisConfiguration.AddItem(pPL->m_lstPZ[j * pPL->m_lstR.Num() + i]);
                        }
                        polyakovOmgR.AddItem(thisConfiguration);
                        polyIn.AddItem(pPL->m_lstLoopZInner[j]);
                        polyOut.AddItem(pPL->m_lstLoopZ[j]);
                    }
                    lstPolyInZ.AddItem(polyIn);
                    lstPolyOutZ.AddItem(polyOut);
                    WriteStringFileComplexArray2(sFileNameWrite3, polyakovOmgR);
                }
            }
            break;
            case EDJKSR_Chiral:
            {

            }
            break;
            case EDJKSR_AngularMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGS, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGChen, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGPot, uiListIdx, REM);
            }
            break;
            case EDJKSR_ChiralAndFermionMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGS, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGChen, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf, uiListIdx, REM);
                _CLG_EXPORT_ANGULAR(pJG, JGPot, uiListIdx, REM);

                _CLG_EXPORT_CHIRAL(pCCu, ChiralKS, uiListIdx, REM);
                if (pCCu->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCu, ConnectSusp, uiListIdx, REM);
                }


                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma1, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma2, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma3, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma4, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma5, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma51, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma52, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma53, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSGamma54, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSSigma13, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSSigma14, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSSigma23, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSSigma24, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCu, CMTKSSigma34, uiListIdx, REM);

                _CLG_EXPORT_CHIRAL(pCCd, ChiralKS, uiListIdx, REM);
                if (pCCd->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCd, ConnectSusp, uiListIdx, REM);
                }
                
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma1, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma2, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma3, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma4, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma5, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma51, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma52, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma53, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSGamma54, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSSigma13, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSSigma14, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSSigma23, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSSigma24, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCd, CMTKSSigma34, uiListIdx, REM);

                _CLG_EXPORT_CHIRAL(pCCs, ChiralKS, uiListIdx, REM);
                if (pCCs->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCs, ConnectSusp, uiListIdx, REM);
                }

                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma1, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma2, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma3, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma4, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma5, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma51, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma52, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma53, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSGamma54, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSSigma13, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSSigma14, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSSigma23, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSSigma24, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pCCs, CMTKSSigma34, uiListIdx, REM);

                _CLG_EXPORT_CHIRAL(pFAu, OrbitalKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAu, SpinKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAu, PotentialKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAd, OrbitalKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAd, SpinKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAd, PotentialKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAs, OrbitalKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAs, SpinKS, uiListIdx, REM);
                _CLG_EXPORT_CHIRAL(pFAs, PotentialKS, uiListIdx, REM);

                if (uiListIdx == iListStart)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCCu->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(F(0.5) * _hostsqrt(static_cast<Real>(pCCu->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFileRealArray(sRadiousFile, lstRadius);
                }

                //if (bMeasureCCS)
                //{
                //    CCString sFileNameWriteCCS;
                //    sFileNameWriteCCS.Format(_T("%s_lightCCS_Nt%d_REM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiListIdx);
                //    WriteStringFileComplexArray(sFileNameWriteCCS, pCCSLight->m_lstResults);
                //    sFileNameWriteCCS.Format(_T("%s_heavyCCS_Nt%d_REM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiListIdx);
                //    WriteStringFileComplexArray(sFileNameWriteCCS, pCCSHeavy->m_lstResults);
                //}
            }
            break;
            //case EDJKS_PlaqutteEnergy:
            //{
            //    //CCString sFileName;
            //    //sFileName.Format(_T("%s_plaqutte.csv"), sCSVSavePrefix.c_str());
            //    //WriteStringFileRealArray(sFileName, pPE->m_lstData);
            //}
            //break;
            //case EDJKS_VR:
            //    {
            //        CCString sCSVFile;
            //        sCSVFile.Format(_T("%s_VR_Nt%d_REM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiListIdx);
            //        TArray<TArray<CLGComplex>> vrs;
            //        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            //        {
            //            TArray<CLGComplex> thisConfiguration;
            //            for (INT i = 0; i < pWilson->m_lstR.Num(); ++i)
            //            {
            //                for (UINT t = 0; t < (_HC_Lt - 1); ++t)
            //                {
            //                    thisConfiguration.AddItem(pWilson->m_lstC[j][i][t]);
            //                }
            //            }
            //            vrs.AddItem(thisConfiguration);
            //        }
            //        WriteStringFileComplexArray2(sCSVFile, vrs);
            //    }
            //break;
            //case EDJKS_DoubleToFloat:
            //    {
            //        //do nothing
            //    }
            //    break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    //if (NULL != pStaple)
    //{
    //    appSafeDelete(pStaple);
    //}

    if (EDJKSR_CheckMD5 == eJob)
    {
        if (NULL != pF1u)
        {
            pF1u->Return();
            pF2u->Return();
            pF1d->Return();
            pF2d->Return();
            pF1s->Return();
            pF2s->Return();
        }

        appQuitCLG();
        return 0;
    }

    switch (eJob)
    {
        case EDJKSR_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyIn);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOut);

            sFileNameWrite1.Format(_T("%s_polyakovZ_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakovZ_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyInZ);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOutZ);
        }
        break;
        default:
            break;
    }

    appGeneral(_T("\n(*"));
    appPopLogDate();

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));
    if (NULL != pF1u)
    {
        pF1u->Return();
        pF2u->Return();
        pF1d->Return();
        pF2d->Return();
        pF1s->Return();
        pF2s->Return();
    }

    appQuitCLG();

    return 0;
}


