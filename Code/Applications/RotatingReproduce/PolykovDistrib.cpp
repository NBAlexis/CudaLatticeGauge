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

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJob eJob = __STRING_TO_ENUM(EDistributionJob, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    Real fBeta = F(0.0);
    params.FetchValueReal(_T("GaugeBate"), fBeta);


    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
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

    for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
    {

        appGeneral(_T("(* ==== Omega(%d) ========= *)\n"), uiOmega);
        pPL->Reset();
        pCC->Reset();
        pJG->Reset();
        pJF->Reset();

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
                }
                break;
            }

            if ((iEndN - uiN + 1) % uiNewLine == 0)
            {
                appGeneral(_T("\n="));
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
                    }

                    appGeneral(_T("};\n"));
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
                    }

                    appGeneral(_T("};\n"));
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
                    }

                    appGeneral(_T("};\n"));
                }

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

            appGeneral(_T("\njgallchi={"));
            for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
            {
                appGeneral(_T("Mean[jg%d*jg%d] - Mean[jg%d]*Mean[jg%d]%s"), uiOmega, uiOmega, uiOmega, uiOmega, uiOmega == iEndOmega ? _T("") : _T(", "));
            }
            appGeneral(_T("}\n\n"));

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[jgallmean][[%d]]%s"), i + 1, (i == (pJG->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));

            appGeneral(_T("\nListLinePlot[{"));
            for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
            {
                appGeneral(_T("Transpose[jgallchi][[%d]]%s"), i + 1, (i == (pJG->m_lstR.Num() - 1)) ? _T("") : _T(", "));
            }
            appGeneral(_T("}, PlotRange -> All]\n\n"));
        }
        break;
    }

    appSetLogDate(TRUE);

    appGeneral(_T("\n=====================================\n========= finished! ==========\n"));
    if (bSmearing)
    {
        delete pStaple;
    }
    appQuitCLG();

    return 0;
}
