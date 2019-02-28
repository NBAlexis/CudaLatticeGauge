//=============================================================================
// FILENAME : SampleProducer.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/26/2019 nbale]
//=============================================================================

#include "SampleProducer.h"

int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("Applications/SampleProducer.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/Applications/SampleProducer.yaml"), params);
#endif

    CCString sJobType = _T("ESJ_BuildSample");
    params.FetchStringValue(_T("SampleJob"), sJobType);
    ESampleJob eJob = __STRING_TO_ENUM(ESampleJob, sJobType);

    switch (eJob)
    {
    case ESJ_BuildGauge:
    {
        CParameters buildGaugeParams = params.GetParameter(_T("BuildGauge"));
        appSetupLog(params);
        if (!appInitialCLG(buildGaugeParams))
        {
            return 1;
        }

        //we calculate staple energy from beta
        CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
        if (NULL == pAction)
        {
            return 1;
        }

        CCString sFileName;
        for (UINT i = 0; i < 100; ++i)
        {
            Real fBeta = i * F(0.1);
            sFileName.Format(_T("Gauge_%d.con"), i);
            pAction->SetBeta(fBeta);
            //Equilibration
            appGetLattice()->m_pUpdator->Update(10, FALSE);
            appGetLattice()->m_pGaugeField->SaveToFile(sFileName);
        }
    }
    break;
    case ESJ_BuildSample:
    {
        TArray<INT> sampleParam;
        TArray<Real> kappaRange;
        params.FetchValueArrayINT(_T("SampleParameter"), sampleParam);
        params.FetchValueArrayReal(_T("KappaRange"), kappaRange);
        if (4 != sampleParam.Num()
            || sampleParam[1] < 0
            || sampleParam[1] > 100
            || sampleParam[2] < 0
            || sampleParam[2] > 99
            || sampleParam[3] <= sampleParam[2])
        {
            sampleParam.RemoveAll();
            sampleParam.AddItem(1);
            sampleParam.AddItem(10);
            sampleParam.AddItem(0);
            sampleParam.AddItem(100);
        }

        if (2 != kappaRange.Num()
            || kappaRange[0] < F(0.01)
            || kappaRange[1] > F(0.21)
            || kappaRange[1] < kappaRange[0])
        {
            kappaRange.RemoveAll();
            kappaRange.AddItem(F(0.1));
            kappaRange.AddItem(F(0.2));
        }

        appSetupLog(params);
        if (!appInitialCLG(params))
        {
            return 1;
        }

        CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
        if (NULL == pAction)
        {
            return 1;
        }
        CFieldFermionWilsonSquareSU3* pFermionField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2));
        if (NULL == pFermionField)
        {
            return 1;
        }
        CFieldFermionWilsonSquareSU3* pResult = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());
        CFieldFermionWilsonSquareSU3* pResultCheck = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());
        CCString sFileName;
        CCString sSampleFileName;
        //4 * (1 + linkNum * 18 + siteNum * 48)
        //siteNum = 8192
        //linkNum = 32768
        //linkNum * 18 = 589824
        //siteNum * 48 = 393216
        //(1 + linkNum * 18 + siteNum * 48) = 983041
        //3932164
        //BYTE byData[3932164]; this is stack overflow
        BYTE * byData = (BYTE*)malloc(3932168);
        FLOAT fBeta = 0.0f;

        for (INT i = sampleParam[2]; i < sampleParam[3]; ++i)
        {
            //set beta
            fBeta = static_cast<FLOAT>(i * F(0.1));
#if _CLG_DEBUG
            sFileName.Format(_T("../Release/SampleProducer/Gauges/Gauge_%d.con"), i);
#else
            sFileName.Format(_T("SampleProducer/Gauges/Gauge_%d.con"), i);
#endif
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pAction->SetBeta(static_cast<Real>(fBeta));

            for (INT j = 0; j < sampleParam[1]; ++j)
            {
                //set kappa (0.1 - 0.2) only support single float
                FLOAT fKappa = static_cast<FLOAT>(kappaRange[0] + (kappaRange[1] - kappaRange[0]) * GetRandomReal());
                pFermionField->SetKai(static_cast<Real>(fKappa));

                appGetLattice()->m_pUpdator->Update(3, FALSE);

                //find the solution
                pFermionField->CopyTo(pResult);
                pResult->InverseDDdagger(appGetLattice()->m_pGaugeField);

                //check result
                pResult->CopyTo(pResultCheck);
                pResultCheck->DDdagger(appGetLattice()->m_pGaugeField);
                pResultCheck->AxpyMinus(pFermionField);
                Real fError = _cuCabsf(pResultCheck->Dot(pResultCheck));
                Real fLengthOfB = pFermionField->Dot(pFermionField).x;
                Real fLengthOfX = pResult->Dot(pResult).x;
                if (!isnan(fError) 
                    && fError < F(0.0000001) 
                    && !isnan(fLengthOfB)
                    && fLengthOfB < F(1000000.0) 
                    && !isnan(fLengthOfX)
                    && fLengthOfX < F(1000000.0))
                {
                    //save the sample
                    memcpy(byData, &fBeta, sizeof(FLOAT));
                    memcpy(byData + sizeof(FLOAT), &fKappa, sizeof(FLOAT));
                    UINT uiSize = 0;
                    BYTE* gaugeData = appGetLattice()->m_pGaugeField->CopyDataOut(uiSize);
                    memcpy(byData + sizeof(FLOAT) * 2, gaugeData, 589824 * sizeof(FLOAT));
                    free(gaugeData);
                    BYTE* fieldBData = pFermionField->CopyDataOut(uiSize);
                    memcpy(byData + sizeof(FLOAT) * (589824 + 2), fieldBData, 196608 * sizeof(FLOAT));
                    free(fieldBData);
                    BYTE* fieldXData = pResult->CopyDataOut(uiSize);
                    memcpy(byData + sizeof(FLOAT) * (196608 + 589824 + 2), fieldXData, 196608 * sizeof(FLOAT));
                    free(fieldXData);

                    sSampleFileName.Format(_T("Sample_%d.con"), sampleParam[0] + i * sampleParam[1] + j);
                    appGetFileSystem()->WriteAllBytes(sSampleFileName.c_str(), byData, 3932168);

                    appGeneral(_T("Saved sample kappa=%f, id=%d, lengthB=%f, lengthX=%f\n"), fKappa, sampleParam[0] + i * sampleParam[1] + j, fLengthOfB, fLengthOfX);
                }
                else
                {
                    //give up
                    appGeneral(_T("Solver give er = %1.12f > 0.0000001, give up! (kappa=%f)\n"), fError, fKappa);
                    j -= 1;

                    if (isnan(appGetLattice()->m_pGaugeField->Dot(appGetLattice()->m_pGaugeField).x))
                    {
                        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
                        pAction->SetBeta(fBeta);
                    }
                }
            }
        }

        free(byData);
    }
    break;
    case ESJ_CheckSample:
    {
        //the range to check
        TArray<INT> checkRange;
        params.FetchValueArrayINT(_T("CheckRange"), checkRange);
        if (2 != checkRange.Num()
         || checkRange[1] < checkRange[0])
        {
            checkRange[0] = 1;
            checkRange[1] = 1000;
        }

        appSetupLog(params);
        if (!appInitialCLG(params))
        {
            return 1;
        }

        CFieldFermionWilsonSquareSU3* pFermionField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2));
        if (NULL == pFermionField)
        {
            return 1;
        }
        CFieldFermionWilsonSquareSU3* pResult = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());
        CFieldFermionWilsonSquareSU3* pResultCheck = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());

        CCString sFileName;
        TArray<CCString> problemSamples;
        for (INT i = checkRange[0]; i <= checkRange[1]; ++i)
        {
#if _CLG_DEBUG
            sFileName.Format(_T("../Release/SampleProducer/Samples/Sample_%d.con"), i);
#else
            sFileName.Format(_T("SampleProducer/Samples/Sample_%d.con"), i);
#endif

            UINT uiSize = 0;
            BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
            FLOAT kappa[1];
            FLOAT beta[1];
            memcpy(beta, byData, sizeof(FLOAT));
            memcpy(kappa, byData + sizeof(FLOAT), sizeof(FLOAT));
            appGetLattice()->m_pGaugeField->InitialWithByte(byData + sizeof(FLOAT) * 2);
            pResult->InitialWithByte(byData + sizeof(FLOAT) * (589824 + 2));
            pResultCheck->InitialWithByte(byData + sizeof(FLOAT) * (196608 + 589824 + 2));
            pResultCheck->SetKai(kappa[0]);
            CCString sErr;
            Real fLengthB = pResult->Dot(pResult).x;
            if (isnan(fLengthB) || fLengthB > F(1000000.0))
            {
                CCString sErr1;
                sErr1.Format(_T("LengthB invalid %f"), fLengthB);
                sErr += sErr1;
            }
            Real fLengthX = pResultCheck->Dot(pResultCheck).x;
            if (isnan(fLengthX) || fLengthX > F(1000000.0))
            {
                CCString sErr1;
                sErr1.Format(_T("LengthX invalid %f"), fLengthX);
                sErr += sErr1;
            }
            pResultCheck->DDdagger(appGetLattice()->m_pGaugeField);
            pResultCheck->AxpyMinus(pResult);
            Real fErr = pResultCheck->Dot(pResultCheck).x;
            if (isnan(fErr) || fErr > F(0.00000010001))
            {
                CCString sErr1;
                sErr1.Format(_T("Err invalid %f(>1e-7)"), fErr);
                sErr += sErr1;
            }
            if (sErr.GetLength() > 0)
            {
                CCString sFile;
                sFile.Format(_T("Sample_%d.con is invalid(beta=%f,kappa=%f): %s\n"), i, beta[0], kappa[0], sErr.c_str());
                problemSamples.AddItem(sFile);
            }
            //appGeneral(_T("%d=%f,%f\n"), i, fLengthB, fLengthX);

            appGeneral(_T("="));
            free(byData);
        }
        appGeneral(_T("\n"));
        if (problemSamples.Num() < 1)
        {
            appGeneral(_T("All checked, no problem!\n"));
        }
        else
        {
            for (INT i = 0; i < problemSamples.Num(); ++i)
            {
                appGeneral(problemSamples[i]);
            }
        }
    }
    break;
    case ESJ_FillSample:
    {
        TArray<INT> fillList;
        TArray<INT> gaugeList;
        TArray<Real> kappaRange;
        params.FetchValueArrayINT(_T("SampleList"), fillList);
        params.FetchValueArrayINT(_T("SampleListGauge"), gaugeList);

        params.FetchValueArrayReal(_T("KappaRange"), kappaRange);
        if (fillList.Num() != gaugeList.Num())
        {
            appCrucial(_T("length of SampleList != SampleListGauge"));
            return 1;
        }
        if (2 != kappaRange.Num()
            || kappaRange[0] < F(0.01)
            || kappaRange[1] > F(0.21)
            || kappaRange[1] < kappaRange[0])
        {
            kappaRange.RemoveAll();
            kappaRange.AddItem(F(0.1));
            kappaRange.AddItem(F(0.2));
        }

        appSetupLog(params);
        if (!appInitialCLG(params))
        {
            return 1;
        }

        CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
        if (NULL == pAction)
        {
            return 1;
        }
        CFieldFermionWilsonSquareSU3* pFermionField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2));
        if (NULL == pFermionField)
        {
            return 1;
        }
        CFieldFermionWilsonSquareSU3* pResult = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());
        CFieldFermionWilsonSquareSU3* pResultCheck = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermionField->GetCopy());

        //4 * (1 + linkNum * 18 + siteNum * 48)
        //siteNum = 8192
        //linkNum = 32768
        //linkNum * 18 = 589824
        //siteNum * 48 = 393216
        //(1 + linkNum * 18 + siteNum * 48) = 983041
        //3932164
        //BYTE byData[3932164]; this is stack overflow
        BYTE * byData = (BYTE*)malloc(3932164);
        Real fBeta = F(0.0);
        CCString sFileName;
        CCString sFileNameGauge;
        for (INT i = 0; i < fillList.Num(); ++i)
        {
            sFileName.Format(_T("Sample_%d.con"), fillList[i]);
            FLOAT fKappa = static_cast<FLOAT>(kappaRange[0] + (kappaRange[1] - kappaRange[0]) * GetRandomReal());
            pFermionField->SetKai(static_cast<Real>(fKappa));

            INT iGauge = gaugeList[i];
            if (iGauge < 0 || iGauge > 99)
            {
                iGauge = 50;
            }
            fBeta = iGauge * F(0.1);
#if _CLG_DEBUG
            sFileNameGauge.Format(_T("../Release/SampleProducer/Gauges/Gauge_%d.con"), iGauge);
#else
            sFileNameGauge.Format(_T("SampleProducer/Gauges/Gauge_%d.con"), iGauge);
#endif
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileNameGauge, EFFT_CLGBin);
            pAction->SetBeta(fBeta);


            appGetLattice()->m_pUpdator->Update(3, FALSE);

            //find the solution
            pFermionField->CopyTo(pResult);
            pResult->InverseDDdagger(appGetLattice()->m_pGaugeField);

            //check result
            pResult->CopyTo(pResultCheck);
            pResultCheck->DDdagger(appGetLattice()->m_pGaugeField);
            pResultCheck->AxpyMinus(pFermionField);
            Real fError = _cuCabsf(pResultCheck->Dot(pResultCheck));
            Real fLengthOfB = pFermionField->Dot(pFermionField).x;
            Real fLengthOfX = pResult->Dot(pResult).x;

            if (!isnan(fError)
                && fError < F(0.0000001)
                && !isnan(fLengthOfB)
                && fLengthOfB < F(1000000.0)
                && !isnan(fLengthOfX)
                && fLengthOfX < F(1000000.0))
            {
                //save the sample
                memcpy(byData, &fKappa, sizeof(FLOAT));
                UINT uiSize = 0;
                BYTE* gaugeData = appGetLattice()->m_pGaugeField->CopyDataOut(uiSize);
                memcpy(byData + sizeof(FLOAT), gaugeData, 589824 * sizeof(FLOAT));
                free(gaugeData);
                BYTE* fieldBData = pFermionField->CopyDataOut(uiSize);
                memcpy(byData + sizeof(FLOAT) * (589824 + 1), fieldBData, 196608 * sizeof(FLOAT));
                free(fieldBData);
                BYTE* fieldXData = pResult->CopyDataOut(uiSize);
                memcpy(byData + sizeof(FLOAT) * (196608 + 589824 + 1), fieldXData, 196608 * sizeof(FLOAT));
                free(fieldXData);

                appGetFileSystem()->WriteAllBytes(sFileName.c_str(), byData, 3932164);

                appGeneral(_T("Saved sample kappa=%f, file=%s, lengthB=%f, lengthX=%f\n"), fKappa, sFileName.c_str(), fLengthOfB, fLengthOfX);
            }
            else
            {
                //give up
                appGeneral(_T("Solver give er = %1.12f > 0.0000001, give up! (kappa=%f)\n"), fError, fKappa);
                i -= 1;
            }
        }

        free(byData);

    }
    break;
    }

    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
