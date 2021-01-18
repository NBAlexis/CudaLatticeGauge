//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"

#define _CLG_EXPORT_CHIRAL(measureName, lstName) \
CCString sFileNameWrite##measureName##lstName = _T("%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##All = _T("%s_condensate"); \
sFileNameWrite##measureName##lstName = sFileNameWrite##measureName##lstName + _T(#measureName) + _T(#lstName) + _T("_Nt%d_O%d.csv"); \
sFileNameWrite##measureName##lstName##All = sFileNameWrite##measureName##lstName##All + _T(#measureName) + _T(#lstName) + _T("_Nt%d_All_O%d.csv"); \
sFileNameWrite##measureName##lstName.Format(sFileNameWrite##measureName##lstName, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
sFileNameWrite##measureName##lstName##All.Format(sFileNameWrite##measureName##lstName##All, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
TArray<TArray<CLGComplex>> lstName##measureName##OverR; \
TArray<CLGComplex> lstName##measureName##All; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<CLGComplex> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lstCond[lstName][j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##measureName##OverR.AddItem(thisConfiguration); \
    lstName##measureName##All.AddItem(measureName->m_lstCondAll[lstName][j]); \
} \
WriteStringFileComplexArray2(sFileNameWrite##measureName##lstName, lstName##measureName##OverR); \
WriteStringFileComplexArray(sFileNameWrite##measureName##lstName##All, lstName##measureName##All);


__DEFINE_ENUM(EDistributionJobKS,
    EDJKS_Polyakov,
    EDJKS_Chiral,
    EDJKS_AngularMomentum,
    EDJKS_ChiralAndFermionMomentum,
    EDJKS_PlaqutteEnergy,
    )


enum { kExportDigital = 20, };

#if !_CLG_WIN

void _gcvt_s(TCHAR* buff, UINT uiBuffLength, Real fVaule, UINT uiDigit)
{
    static TCHAR tmpBuff[10];
    appSprintf(tmpBuff, 10, _T("%s.%df"), _T("%"), uiDigit);
    appSprintf(buff, uiBuffLength, tmpBuff, fVaule);
}

#endif

void WriteStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->WriteAllText(sFileName, sContent);
}

void WriteStringFileRealArray(const CCString& sFileName, const TArray<Real>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }
    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i], iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;
        if (i != lst.GetCount() - 1)
        {
            file << _T(",");
        }
    }
    file.flush();
    file.close();
}

void WriteStringFileRealArray2(const CCString& sFileName, const TArray<TArray<Real>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }
    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j], iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            file << _T(" ");
            file << sReal;
            if (j != lst[i].GetCount() - 1)
            {
                file << _T(",");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

void WriteStringFileComplexArray(const CCString& sFileName, const TArray<CLGComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }
    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }

        file << _T(" ");
        file << sReal;
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}

void WriteStringFileComplexArray2(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }
    
    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            _gcvt_s(str, 50, lst[i][j].y, iDigital);
            CCString sImg = CCString(str);
            sImg = sImg.Replace(_T("e"), _T("*^"));
            CCString sMid = _T(" + ");
            if (sImg.Left(1) == _T("-"))
            {
                sImg = sImg.Right(sImg.GetLength() - 1);
                sMid = _T(" - ");
            }
            file << _T(" ");
            file << sReal;
            file << sMid;
            file << sImg;
            if (j == lst[i].GetCount() - 1)
            {
                file << _T(" I");
            }
            else
            {
                file << _T(" I,");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}

INT Measurement(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 0;
    params.FetchValueINT(_T("StartOmega"), iVaule);
    UINT iStartOmega = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("EndOmega"), iVaule);
    UINT iEndOmega = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("FermionMomentum"), iVaule);
    UBOOL bJF = 0 != iVaule;

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    //iVaule = 1;
    //params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    //UBOOL bCheckGaugeFixing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJobKS eJob = __STRING_TO_ENUM(EDistributionJobKS, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    Real fBeta = F(0.0);
    params.FetchValueReal(_T("GaugeBate"), fBeta);

    Real fOmega = F(0.1);
    params.FetchValueReal(_T("OmegaRange"), fOmega);
    fOmega = fOmega / iEndOmega;

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

#pragma endregion

    UINT uiMaxL = (_HC_Lx + 1) / 2 - 1;
    uiMaxL = uiMaxL * uiMaxL;
    TArray<TArray<CCString>> r_omega_idx;
    for (UINT i = 0; i < iEndOmega * iEndOmega * uiMaxL; ++i)
    {
        TArray<CCString> newlst;
        r_omega_idx.AddItem(newlst);
    }
    TArray<Real> lstR;
    TArray<TArray<CLGComplex>> lstPolyIn;
    TArray<TArray<CLGComplex>> lstPolyOut;

    CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureAngularMomentumKS* pFALight = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CMeasureAngularMomentumKS* pFAHeavy = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));

    CMeasureAction* pPE = dynamic_cast<CMeasureAction*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    //CActionFermionWilsonNf2* pAF = dynamic_cast<CActionFermionWilsonNf2*>(appGetLattice()->m_pActionList[1]);

    CActionGaugePlaquetteRotating* pAG = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    CFieldFermionKSSU3* pF1Heavy = NULL;
    CFieldFermionKSSU3* pF2Heavy = NULL;

    if (EDJKS_ChiralAndFermionMomentum == eJob
        || (EDJKS_AngularMomentum == eJob && bJF)
        || EDJKS_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
        pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
    }

    appSetLogDate(FALSE);


    for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
    {
        CCommonData::m_fOmega = fOmega * uiOmega;
        if (NULL != pAG)
        {
            pAG->SetOmega(CCommonData::m_fOmega);
        }
        appGeneral(_T("(* ==== Omega(%f) ========= *)\n"), fOmega * uiOmega);
        pPL->Reset();
        pCCLight->Reset();
        pCCHeavy->Reset();
        pFALight->Reset();
        pFAHeavy->Reset();
        pCCLight->SetFieldCount(iFieldCount);
        pCCHeavy->SetFieldCount(iFieldCount);
        pFALight->SetFieldCount(iFieldCount);
        pFAHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/O%d/%sR_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            else
            {
                sFileName.Format(_T("%sR_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            
            switch (eJob)
            {
                case EDJKS_Polyakov:
                {
                    pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EDJKS_Chiral:
                {
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1Light->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Light->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Light->FixBoundary();
                        pF1Light->CopyTo(pF2Light);
                        pF1Light->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Light->FixBoundary();

                        pCCLight->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        if (bZ4)
                        {
                            pF1Heavy->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Heavy->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Heavy->FixBoundary();
                        pF1Heavy->CopyTo(pF2Heavy);
                        pF1Heavy->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Heavy->FixBoundary();

                        pCCHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);
                    }

                }
                break;
                case EDJKS_AngularMomentum:
                {

                }
                break;
                case EDJKS_ChiralAndFermionMomentum:
                {
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1Light->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Light->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Light->FixBoundary();
                        pF1Light->CopyTo(pF2Light);
                        pF1Light->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Light->FixBoundary();

                        pCCLight->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        pFALight->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        if (bZ4)
                        {
                            pF1Heavy->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Heavy->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Heavy->FixBoundary();
                        pF1Heavy->CopyTo(pF2Heavy);
                        pF1Heavy->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Heavy->FixBoundary();

                        pCCHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);

                        pFAHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                }
                break;
                case EDJKS_PlaqutteEnergy:
                {
                    pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
            }

            if ((iEndN - uiN + 1) % uiNewLine == 0)
            {
                appSetLogDate(TRUE);
                appGeneral(_T("\n="));
                appSetLogDate(FALSE);
            }
            else
            {
                appGeneral(_T("="));
            }
            
        }
        appGeneral(_T("\n*)\n"));

#pragma endregion

        switch (eJob)
        {
            case EDJKS_Polyakov:
            {
                CCString sFileNameWrite1;
                CCString sFileNameWrite2;
                sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_R.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pPL->m_lstR.Num() == pPL->m_lstP.Num());
                
                if (uiOmega == iStartOmega)
                {
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        lstR.AddItem(_hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
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
            }
            break;
            case EDJKS_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCCLight->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pCCLight->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFileRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJKS_AngularMomentum:
            {

            }
            break;
            case EDJKS_ChiralAndFermionMomentum:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp);
                _CLG_EXPORT_CHIRAL(pFALight, OrbitalKS);
                _CLG_EXPORT_CHIRAL(pFALight, SpinKS);
                _CLG_EXPORT_CHIRAL(pFALight, PotentialKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, OrbitalKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, SpinKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, PotentialKS);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCCLight->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pCCLight->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFileRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJKS_PlaqutteEnergy:
            {
                CCString sFileName;
                sFileName.Format(_T("%s_plaqutte.csv"), sCSVSavePrefix.c_str());
                WriteStringFileRealArray(sFileName, pPE->m_lstData);
            }
            break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    switch (eJob)
    {
        case EDJKS_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyIn);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOut);
        }
        break;
        case EDJKS_Chiral:
        {
            //nothing to do
        }
        break;
        case EDJKS_AngularMomentum:
        {
            //nothing to do
        }
        break;
        case EDJKS_ChiralAndFermionMomentum:
        {
            //nothing to do
        }
        break;
        default:
            break;
    }

    appGeneral(_T("\n(*"));
    appSetLogDate(TRUE);

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));
    if (NULL != pF1Light)
    {
        pF1Light->Return();
        pF2Light->Return();
        pF1Heavy->Return();
        pF2Heavy->Return();
    }

    appQuitCLG();

    return 0;
}


