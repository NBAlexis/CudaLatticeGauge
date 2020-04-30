//=============================================================================
// FILENAME : PolykovDistrib.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [07/09/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

__DEFINE_ENUM(EDistributionJob,
    EDJ_Polyakov,
    EDJ_Chiral,
    EDJ_AngularMomentum,
    EDJ_ChiralAndFermionMomentum,
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

CCString ExportComplexArray(const TArray<CLGComplex> lst)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    CCString sSaveString = _T("");
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        TCHAR str[50];
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
        sSaveString = sSaveString + _T(" ") + sReal + sMid + sImg 
            + ((i == lst.GetCount() - 1) ? _T(" I") : _T(" I,"));
    }

    return sSaveString;
}

CCString ExportComplexArray2(const TArray<TArray<CLGComplex>> lst)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    CCString sSaveString = _T("");
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            TCHAR str[50];
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
            sSaveString = sSaveString + _T(" ") + sReal + sMid + sImg
            + ((j == lst[i].GetCount() - 1) ? _T(" I") : _T(" I,"));
        }
        sSaveString = sSaveString + _T("\n");
    }

    return sSaveString;
}

CCString ExportRealArray(const TArray<Real> lst)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    CCString sSaveString = _T("");
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        TCHAR str[50];
        _gcvt_s(str, 50, lst[i], iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        sSaveString = sSaveString + _T(" ") + sReal
        + ((i == lst.GetCount() - 1) ? _T("") : _T(","));
    }

    return sSaveString;
}

CCString ExportRealArray2(const TArray<TArray<Real>> lst)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    CCString sSaveString = _T("");
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            TCHAR str[50];
            _gcvt_s(str, 50, lst[i][j], iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            sSaveString = sSaveString + _T(" ") + sReal
                + ((j == lst[i].GetCount() - 1) ? _T("") : _T(","));
        }
        sSaveString = sSaveString + _T("\n");
    }
    return sSaveString;
}

void WriteStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->WriteAllText(sFileName, sContent);
}

void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}

INT MeasurePolyakovDist(CParameters& params)
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

    iVaule = 1;
    params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    UBOOL bCheckGaugeFixing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJob eJob = __STRING_TO_ENUM(EDistributionJob, sValue);

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
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureAMomentumStochastic* pJF = dynamic_cast<CMeasureAMomentumStochastic*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));

    appSetLogDate(FALSE);

    CFieldFermionWilsonSquareSU3* pF1 = NULL;
    CFieldFermionWilsonSquareSU3* pF2 = NULL;

    if (EDJ_ChiralAndFermionMomentum == eJob
     || (EDJ_AngularMomentum == eJob && bJF)
     || EDJ_Chiral == eJob)
    {
        pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    }

    for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
    {
        CCommonData::m_fOmega = fOmega * uiOmega;
        appGeneral(_T("(* ==== Omega(%f) ========= *)\n"), fOmega * uiOmega);
        pPL->Reset();
        pCC->Reset();
        pJG->Reset();
        pJF->Reset();
        pCC->SetFieldCount(iFieldCount);
        pJF->SetFieldCount(iFieldCount);

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/O%d/%sRotate_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            else
            {
                sFileName.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            
            switch (eJob)
            {
                case EDJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EDJ_Chiral:
                {
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
                        pF1->InverseD(appGetLattice()->m_pGaugeField);
                        pF1->FixBoundary();

                        pCC->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2,
                            pF1,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                }
                break;
                case EDJ_AngularMomentum:
                {
                    if (bCheckGaugeFixing && NULL != appGetLattice()->m_pGaugeFixing)
                    {
                        Real fError = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
                        if (appAbs(fError) > F(0.000000000001))
                        {
                            appGeneral(_T("Bad Gauge Fixing\n"));
                        }
                    }
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    if (bJF)
                    {
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
                            pF1->InverseD(appGetLattice()->m_pGaugeField);

                            pJF->OnConfigurationAcceptedZ4(
                                appGetLattice()->m_pGaugeField, 
                                NULL, 
                                pF2, 
                                pF1, 
                                0 == i, 
                                iFieldCount == i + 1);
                        }
                    }
                }
                break;
                case EDJ_ChiralAndFermionMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    if (NULL != pJF && NULL != pCC)
                    {
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
                            pF1->InverseD(appGetLattice()->m_pGaugeField);


                            pJF->OnConfigurationAcceptedZ4(
                                appGetLattice()->m_pGaugeField,
                                NULL,
                                pF2,
                                pF1,
                                0 == i,
                                iFieldCount == i + 1);

                            pCC->OnConfigurationAcceptedZ4(
                                appGetLattice()->m_pGaugeField,
                                NULL,
                                pF2,
                                pF1,
                                0 == i,
                                iFieldCount == i + 1);
                        }
                    }
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
            case EDJ_Polyakov:
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
                    WriteStringFile(sFileNameWrite1, ExportRealArray(lstR));
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
                WriteStringFile(sFileNameWrite2, ExportComplexArray2(polyakovOmgR));
            }
            break;
            case EDJ_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCC, Chiral);
                _CLG_EXPORT_CHIRAL(pCC, Gamma1);
                _CLG_EXPORT_CHIRAL(pCC, Gamma2);
                _CLG_EXPORT_CHIRAL(pCC, Gamma3);
                _CLG_EXPORT_CHIRAL(pCC, Gamma4);
                _CLG_EXPORT_CHIRAL(pCC, Gamma5);
                _CLG_EXPORT_CHIRAL(pCC, Gamma45);
                _CLG_EXPORT_CHIRAL(pCC, GammaX);
                _CLG_EXPORT_CHIRAL(pCC, GammaY);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pCC->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFile(sRadiousFile, ExportRealArray(lstRadius));
                }
            }
            break;
            case EDJ_AngularMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                //_CLG_EXPORT_ANGULAR(pJG, JGChenApprox);

                if (bJF && NULL != pJF)
                {
                    _CLG_EXPORT_ANGULAR(pJF, JL);
                    _CLG_EXPORT_ANGULAR(pJF, JS);
                    //_CLG_EXPORT_ANGULAR(pJF, JLPure);
                    //_CLG_EXPORT_ANGULAR(pJF, JLJM);
                    _CLG_EXPORT_ANGULAR(pJF, JPot);
                }

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pJG->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_angularR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFile(sRadiousFile, ExportRealArray(lstRadius));
                }
            }
            break;
            case EDJ_ChiralAndFermionMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                //_CLG_EXPORT_ANGULAR(pJG, JGChenApprox);

                if (NULL != pJF)
                {
                    _CLG_EXPORT_ANGULAR(pJF, JL);
                    _CLG_EXPORT_ANGULAR(pJF, JS);
                    //_CLG_EXPORT_ANGULAR(pJF, JLPure);
                    //_CLG_EXPORT_ANGULAR(pJF, JLJM);
                    _CLG_EXPORT_ANGULAR(pJF, JPot);
                }

                if (NULL != pCC)
                {
                    _CLG_EXPORT_CHIRAL(pCC, Chiral);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma1);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma2);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma3);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma4);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma5);
                    _CLG_EXPORT_CHIRAL(pCC, Gamma45);
                    _CLG_EXPORT_CHIRAL(pCC, GammaX);
                    _CLG_EXPORT_CHIRAL(pCC, GammaY);
                }

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pJG->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_angularR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFile(sRadiousFile, ExportRealArray(lstRadius));
                }
            }
            break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    switch (eJob)
    {
        case EDJ_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFile(sFileNameWrite1, ExportComplexArray2(lstPolyIn));
            WriteStringFile(sFileNameWrite2, ExportComplexArray2(lstPolyOut));
        }
        break;
        case EDJ_Chiral:
        {
            //nothing to do
        }
        break;
        case EDJ_AngularMomentum:
        {
            //nothing to do
        }
        break;
        case EDJ_ChiralAndFermionMomentum:
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
    if (NULL != pF1)
    {
        pF1->Return();
        pF2->Return();
    }

    appQuitCLG();

    return 0;
}


