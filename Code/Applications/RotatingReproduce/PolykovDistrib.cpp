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

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJob eJob = __STRING_TO_ENUM(EDistributionJob, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

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
            sFileName.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
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

                        pCC->OnConfigurationAcceptedZ4(appGetLattice()->m_pGaugeField, NULL, pF2, pF1, 0 == i, iFieldCount == i + 1);
                        pJF->OnConfigurationAcceptedZ4(appGetLattice()->m_pGaugeField, NULL, pF2, pF1, 0 == i, iFieldCount == i + 1);
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
                _CLG_EXPORT_CHIRAL(pCC, Pion);
                _CLG_EXPORT_CHIRAL(pCC, Rhon);
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
                _CLG_EXPORT_ANGULAR(pJG, JGChenApprox);

                if (bJF && NULL != pJF)
                {
                    _CLG_EXPORT_ANGULAR(pJF, JL);
                    _CLG_EXPORT_ANGULAR(pJF, JS);
                    _CLG_EXPORT_ANGULAR(pJF, JLPure);
                    _CLG_EXPORT_ANGULAR(pJF, JLJM);
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
#pragma region Chiral

                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pCC->m_lstR.Num() == pCC->m_lstChiral.Num());

                if (uiOmega == iStartOmega)
                {
                    appGeneral(_T("cr={"));

                    for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), _hostsqrt(static_cast<Real>(pCC->m_lstR[i])), (i == pCC->m_lstR.Num() - 1) ? _T("") : _T(", "));
                        if (pCC->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("Transpose[c%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pCC->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }

                    appGeneral(_T("};\n"));
                }
                else
                {
                    for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
                    {
                        if (pCC->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("Transpose[c%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pCC->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }
                }

                appGeneral(_T("c%d={\n"), uiOmega);

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), pCC->m_lstChiral[j * pCC->m_lstR.Num() + i], (i == pCC->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }

                appGeneral(_T("\n};\n"));

#pragma endregion

#pragma region JF

#pragma region JFL

                appGeneral(_T("jfl%d={\n"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 0; i < pJF->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), pJF->m_lstJL[j * pJF->m_lstR.Num() + i], (i == pJF->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }
                appGeneral(_T("\n};\n"));

                //================== JL Total ==================
                appGeneral(_T("jfl$in%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("%2.12f%s"),
                        pJF->m_lstJLInner[j],
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

                appGeneral(_T("jfl$out%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("%2.12f%s"),
                        pJF->m_lstJLAll[j],
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

#pragma endregion

#pragma region JFS

                appGeneral(_T("jfs%d={\n"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 0; i < pJF->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), pJF->m_lstJS[j * pJF->m_lstR.Num() + i], (i == pJF->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }
                appGeneral(_T("\n};\n"));

#pragma endregion

#pragma endregion
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
#pragma region Chiral

            appGeneral(_T("\nc0all={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("c%d%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("};\n"));

            appGeneral(_T("\nc0allmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[c%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\nc0allchi={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[c%d*c%d] - Mean[c%d]*Mean[c%d]%s"), uiOmega, uiOmega, uiOmega, uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n\n"));

#pragma region r times omega

            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("\ncrw%d=Join["), i);
                    for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                    {
                        appGeneral(_T("%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("];\n\n"));
                }
            }

            appGeneral(_T("\ncrwlist={"));
            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[crw%d]]}"), (0 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                }
            }
            appGeneral(_T("\n}"));

            appGeneral(_T("\ncrwchilist={"));
            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[crw%d*crw%d] - Mean[crw%d]*Mean[crw%d]]}"),
                        (0 == i ? _T("\n") : _T(",")),
                        _hostsqrt(i),
                        i, i, i, i);
                }
            }
            appGeneral(_T("\n}"));

#pragma endregion

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[c0allmean][[%d]]%s"), i + 1, (i == (pCC->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[c0allchi][[%d]]%s"), i + 1, (i == (pCC->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

#pragma endregion

#pragma region JF

#pragma region JFL

            appGeneral(_T("\njflall={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("jfl%d%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("};\n"));

            appGeneral(_T("\njflallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jfl%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

#pragma endregion

#pragma region JFS

            appGeneral(_T("\njfsall={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("jfs%d%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("};\n"));

            appGeneral(_T("\njfsallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jfs%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

#pragma endregion

#pragma endregion

#pragma region jf in out

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pJF->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[jflallmean][[%d]]%s"), i + 1, (i == (pJF->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

#pragma region in and out jl

            appGeneral(_T("\njflinallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jfl$in%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\njfloutallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jfl$out%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\nListLinePlot[{jflinallmean, jfloutallmean}]\n"));

#pragma endregion

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pJF->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[jfsallmean][[%d]]%s"), i + 1, (i == (pJF->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

#pragma endregion

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


