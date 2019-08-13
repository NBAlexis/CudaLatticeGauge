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

INT MeasurePolyakovDist(CParameters& params)
{
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
    params.FetchValueINT(_T("DoSmearing"), iVaule);
    UBOOL bSmearing = 0 != iVaule;

    iVaule = 1;
    params.FetchValueINT(_T("FermionMomentum"), iVaule);
    UBOOL bJF = 0 != iVaule;

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJob eJob = __STRING_TO_ENUM(EDistributionJob, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

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

    UINT uiMaxL = (_HC_Lx + 1) / 2 - 1;
    uiMaxL = uiMaxL * uiMaxL;
    TArray<TArray<CCString>> r_omega_idx;
    for (UINT i = 0; i < iEndOmega * iEndOmega * uiMaxL; ++i)
    {
        TArray<CCString> newlst;
        r_omega_idx.AddItem(newlst);
    }

    CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureAMomentumStochastic* pJF = dynamic_cast<CMeasureAMomentumStochastic*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));

    CFieldGaugeSU3D* pStaple = NULL;
    if (bSmearing)
    {
        pStaple = dynamic_cast<CFieldGaugeSU3D*>(appGetLattice()->m_pGaugeField->GetCopy());
    }

    appSetLogDate(FALSE);

    CFieldFermionWilsonSquareSU3* pF1 = NULL;
    CFieldFermionWilsonSquareSU3* pF2 = NULL;

    if (EDJ_ChiralAndFermionMomentum == eJob)
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

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            sFileName.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);

            if (bSmearing)
            {
                appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField, pStaple);
            }

            switch (eJob)
            {
                case EDJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EDJ_Chiral:
                {
                    pCC->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EDJ_AngularMomentum:
                {
                    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    if (bJF)
                    {
                        pJF->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    }
                }
                break;
                case EDJ_ChiralAndFermionMomentum:
                {
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        pF1->InitialField(EFIT_RandomZ4);
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

        switch (eJob)
        {
            case EDJ_Polyakov:
            {
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pPL->m_lstR.Num() == pPL->m_lstP.Num());

                if (uiOmega == iStartOmega)
                {
                    appGeneral(_T("pr={"));

                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), _hostsqrt(static_cast<Real>(pPL->m_lstR[i])), (i == pPL->m_lstR.Num() - 1) ? _T("") : _T(", "));
                        if (pPL->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("Transpose[p%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pPL->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }

                    appGeneral(_T("};\n"));
                }
                else
                {
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        if (pPL->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("Transpose[p%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pPL->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }
                }

                appGeneral(_T("p%d={\n"), uiOmega);

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        CLGComplex cV = pPL->m_lstP[j * pPL->m_lstR.Num() + i];
                        appGeneral(_T("%2.12f %s %2.12f I%s"), cV.x, cV.y < F(0.0) ? _T("") : _T("+"), cV.y, (i == pPL->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }

                appGeneral(_T("\n};\n"));

                appGeneral(_T("p$in%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    CLGComplex cV = pPL->m_lstLoopInner[j];
                    appGeneral(_T("%2.12f %s %2.12f I%s"), 
                        cV.x, cV.y < F(0.0) ? _T("") : _T("+"), cV.y, 
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

                appGeneral(_T("p$out%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    CLGComplex cV = pPL->m_lstLoop[j];
                    appGeneral(_T("%2.12f %s %2.12f I%s"),
                        cV.x, cV.y < F(0.0) ? _T("") : _T("+"), cV.y,
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

            }
            break;
            case EDJ_Chiral:
            {
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pCC->m_lstR.Num() == pCC->m_lstC.Num());

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
                        appGeneral(_T("%2.12f%s"), pCC->m_lstC[j * pCC->m_lstR.Num() + i], (i == pCC->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }

                appGeneral(_T("\n};\n"));
            }
            break;
            case EDJ_AngularMomentum:
            {
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pJG->m_lstR.Num() == pJG->m_lstJG.Num());

                if (uiOmega == iStartOmega)
                {
                    appGeneral(_T("jgr={"));

                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), _hostsqrt(static_cast<Real>(pJG->m_lstR[i])), (i == pJG->m_lstR.Num() - 1) ? _T("") : _T(", "));
                        if (pJG->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pJG->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }

                    appGeneral(_T("};\n"));
                }
                else
                {
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        if (pJG->m_lstR[i] < uiMaxL)
                        {
                            CCString tobeAdd;
                            tobeAdd.Format(_T("%d][[%d]]"), uiOmega, i + 1);
                            r_omega_idx[uiOmega * uiOmega * pJG->m_lstR[i]].AddItem(tobeAdd);
                        }
                    }
                }

#pragma region JG

                appGeneral(_T("jg%d={\n"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), pJG->m_lstJG[j * pJG->m_lstR.Num() + i], (i == pJG->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }
                appGeneral(_T("\n};\n"));

                //================== L / R ==================
                //================== L / Omega not calculated ==================
                appGeneral(_T("\njgi%d={\n"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("{"));
                    for (INT i = 1; i < pJG->m_lstR.Num(); ++i)
                    {
                        appGeneral(_T("%2.12f%s"), pJG->m_lstJG[j * pJG->m_lstR.Num() + i] / _hostsqrt(pJG->m_lstR[i]), (i == pJG->m_lstR.Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                }
                appGeneral(_T("\n};\n"));

                //================== JG Total ==================
                appGeneral(_T("jg$in%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("%2.12f%s"),
                        pJG->m_lstJGInner[j],
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

                appGeneral(_T("jg$out%d={"), uiOmega);
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    appGeneral(_T("%2.12f%s"),
                        pJG->m_lstJGAll[j],
                        (j == (iEndN - iStartN)) ? _T("") : _T(", "));
                }
                appGeneral(_T("};\n"));

#pragma endregion

                if (bJF)
                {
#pragma region JFL

                    appGeneral(_T("jfl%d={\n"), uiOmega);
                    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                    {
                        appGeneral(_T("{"));
                        for (INT i = 0; i < pJF->m_lstR.Num(); ++i)
                        {
                            appGeneral(_T("%2.12f%s"), pJF->m_lstJL[j * pJG->m_lstR.Num() + i], (i == pJF->m_lstR.Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                    }
                    appGeneral(_T("\n};\n"));

                    appGeneral(_T("\njfli%d={\n"), uiOmega);
                    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                    {
                        appGeneral(_T("{"));
                        for (INT i = 1; i < pJF->m_lstR.Num(); ++i)
                        {
                            appGeneral(_T("%2.12f%s"), pJF->m_lstJL[j * pJF->m_lstR.Num() + i] / _hostsqrt(pJF->m_lstR[i]), (i == pJF->m_lstR.Num() - 1) ? _T("") : _T(", "));
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

                    appGeneral(_T("\njfsi%d={\n"), uiOmega);
                    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                    {
                        appGeneral(_T("{"));
                        for (INT i = 1; i < pJF->m_lstR.Num(); ++i)
                        {
                            appGeneral(_T("%2.12f%s"), pJF->m_lstJS[j * pJF->m_lstR.Num() + i] / _hostsqrt(pJF->m_lstR[i]), (i == pJF->m_lstR.Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("}%s\n"), (j == (iEndN - iStartN)) ? _T("") : _T(","));
                    }
                    appGeneral(_T("\n};\n"));

#pragma endregion

                }

            }
            break;
            case EDJ_ChiralAndFermionMomentum:
            {
#pragma region Chiral

                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pCC->m_lstR.Num() == pCC->m_lstC.Num());

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
                        appGeneral(_T("%2.12f%s"), pCC->m_lstC[j * pCC->m_lstR.Num() + i], (i == pCC->m_lstR.Num() - 1) ? _T("") : _T(", "));
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
        }

        appGeneral(_T("\n"));
    }

    switch (eJob)
    {
        case EDJ_Polyakov:
        {
            appGeneral(_T("\np0all={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("p%d%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("};\n"));

            appGeneral(_T("\np0allmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p%d]]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\np0allchi={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p%d*p%d] - Mean[p%d]*Mean[p%d]]%s"), uiOmega, uiOmega, uiOmega, uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n\n"));

#pragma region in and out

            appGeneral(_T("\npinallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p$in%d]]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\npinallarg={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Arg[Mean[p$in%d]]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\npinallchi={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p$in%d*p$in%d] - Mean[p$in%d]*Mean[p$in%d]]%s"), uiOmega, uiOmega, uiOmega, uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n\n"));

            appGeneral(_T("\npoutallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p$out%d]]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\npoutallarg={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Arg[Mean[p$out%d]]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\npoutallchi={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Abs[Mean[p$out%d*p$out%d] - Mean[p$out%d]*Mean[p$out%d]]%s"), uiOmega, uiOmega, uiOmega, uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n\n"));

#pragma endregion

            #pragma region r times omega

            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("\nprw%d=Join["), i);
                    for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                    {
                        appGeneral(_T("%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("];\n\n"));
                }
            }

            appGeneral(_T("\nprwlist={"));
            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[prw%d]]}"), (0 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                }
            }
            appGeneral(_T("\n}"));

            appGeneral(_T("\nprwchilist={"));
            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[prw%d*prw%d] - Mean[prw%d]*Mean[prw%d]]}"), 
                        (0 == i ? _T("\n") : _T(",")), 
                        _hostsqrt(i), 
                        i, i, i, i);
                }
            }
            appGeneral(_T("\n}"));

            #pragma endregion

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[p0allmean][[%d]]%s"), i + 1, (i == (pPL->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[p0allchi][[%d]]%s"), i + 1, (i == (pPL->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

            appGeneral(_T("\nListLinePlot[{pinallmean, poutallmean}]\n"));
            appGeneral(_T("\nListLinePlot[{pinallarg, poutallarg}, PlotRange -> All]\n"));
            appGeneral(_T("\nListLinePlot[{pinallchi, poutallchi}, PlotRange -> All]\n"));
            appGeneral(_T("\nListPlot[prwlist]\n\n"));
            
        }
        break;
        case EDJ_Chiral:
        {
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
        }
        break;
        case EDJ_AngularMomentum:
        {
#pragma region JG

            appGeneral(_T("\njgall={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("jg%d%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("};\n"));

            appGeneral(_T("\njgallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jg%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

#pragma endregion
            if (bJF)
            {
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
            }

#pragma region r times omega

            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
#pragma region JG
                    appGeneral(_T("\njgrw%d=Join["), i);
                    for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                    {
                        appGeneral(_T("Transpose[jg%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("];\n\n"));

                    appGeneral(_T("\njgirw%d=Join["), i);
                    for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                    {
                        appGeneral(_T("Transpose[jgi%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                    }
                    appGeneral(_T("];\n\n"));

#pragma endregion

                    if (bJF)
                    {
#pragma region JFL

                        appGeneral(_T("\njflrw%d=Join["), i);
                        for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                        {
                            appGeneral(_T("Transpose[jfl%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("];\n\n"));

                        appGeneral(_T("\njflirw%d=Join["), i);
                        for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                        {
                            appGeneral(_T("Transpose[jfli%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("];\n\n"));

#pragma endregion

#pragma region JFS

                        appGeneral(_T("\njfsrw%d=Join["), i);
                        for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                        {
                            appGeneral(_T("Transpose[jfs%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("];\n\n"));

                        appGeneral(_T("\njfsirw%d=Join["), i);
                        for (INT j = 0; j < r_omega_idx[i].Num(); ++j)
                        {
                            appGeneral(_T("Transpose[jfsi%s%s"), r_omega_idx[i][j].c_str(), (j == r_omega_idx[i].Num() - 1) ? _T("") : _T(", "));
                        }
                        appGeneral(_T("];\n\n"));

#pragma endregion
                    }
                }
            }

#pragma region JG

            appGeneral(_T("\njgrwlist={"));
            for (INT i = 0; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[jgrw%d]]}"), (0 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                }
            }
            appGeneral(_T("\n}"));

            appGeneral(_T("\njgirwlist={"));
            for (INT i = 1; i < r_omega_idx.Num(); ++i)
            {
                if (r_omega_idx[i].Num() > 0)
                {
                    appGeneral(_T("%s{%f, Abs[Mean[jgirw%d]]}"), (1 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                }
            }
            appGeneral(_T("\n}"));

#pragma endregion
            if (bJF)
            {
#pragma region JFL

                appGeneral(_T("\njflrwlist={"));
                for (INT i = 0; i < r_omega_idx.Num(); ++i)
                {
                    if (r_omega_idx[i].Num() > 0)
                    {
                        appGeneral(_T("%s{%f, Abs[Mean[jflrw%d]]}"), (0 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                    }
                }
                appGeneral(_T("\n}"));

                appGeneral(_T("\njflirwlist={"));
                for (INT i = 1; i < r_omega_idx.Num(); ++i)
                {
                    if (r_omega_idx[i].Num() > 0)
                    {
                        appGeneral(_T("%s{%f, Abs[Mean[jflirw%d]]}"), (1 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                    }
                }
                appGeneral(_T("\n}"));

#pragma endregion

#pragma region JFS

                appGeneral(_T("\njfsrwlist={"));
                for (INT i = 0; i < r_omega_idx.Num(); ++i)
                {
                    if (r_omega_idx[i].Num() > 0)
                    {
                        appGeneral(_T("%s{%f, Abs[Mean[jfsrw%d]]}"), (0 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                    }
                }
                appGeneral(_T("\n}"));

                appGeneral(_T("\njfsirwlist={"));
                for (INT i = 1; i < r_omega_idx.Num(); ++i)
                {
                    if (r_omega_idx[i].Num() > 0)
                    {
                        appGeneral(_T("%s{%f, Abs[Mean[jfsirw%d]]}"), (1 == i ? _T("\n") : _T(",")), _hostsqrt(i), i);
                    }
                }
                appGeneral(_T("\n}"));

#pragma endregion
            }

#pragma endregion


            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[jgallmean][[%d]]%s"), i + 1, (i == (pJG->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

#pragma region in and out

            appGeneral(_T("\njginallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jg$in%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\njgoutallmean={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jg$out%d]%s"), uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n"));

            appGeneral(_T("\nListLinePlot[{jginallmean, jgoutallmean}]\n"));

#pragma endregion

            if (bJF)
            {
                appGeneral(_T("\nListLinePlot[{"));
                for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                {
                    appGeneral(_T("Transpose[jflallmean][[%d]]%s"), i + 1, (i == (pJG->m_lstR.Num() - 1)) ? _T("") : _T(", "));
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
                for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                {
                    appGeneral(_T("Transpose[jfsallmean][[%d]]%s"), i + 1, (i == (pJG->m_lstR.Num() - 1)) ? _T("") : _T(", "));
                }
                appGeneral(_T("}, PlotRange -> All]\n\n"));
            }
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
    }

    appGeneral(_T("\n(*"));
    appSetLogDate(TRUE);

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));
    if (EDJ_ChiralAndFermionMomentum == eJob)
    {
        pF1->Return();
        pF2->Return();
    }
    if (bSmearing)
    {
        delete pStaple;
    }
    appQuitCLG();

    return 0;
}
