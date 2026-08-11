//=============================================================================
// FILENAME : TestFermionMatrix.cpp
// 
// DESCRIPTION:
//
//     Test the operations on fields
//
// REVISION:
//  [mm/dd/yy]
//  [08/26/2023 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestFermionDMatrixKS(CParameters& param)
{
    //TArray<UINT> siteIndexes;
    //param.FetchValueArrayUINT(_T("Sites"), siteIndexes);
    CFieldFermionKSSU3* pField = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetFieldById(2));
    UINT uiSiteCount = pField->GetSiteCount();
    CLGComplex* toprint = (CLGComplex*)malloc(sizeof(CLGComplex) * uiSiteCount * uiSiteCount);
    CFieldFermionKSSU3* pFieldCopy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pField->m_byFieldId, _T(__FILE__), __LINE__));

    //CFieldGaugeSU3* pSU3Gauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);

    for (UINT i = 0; i < uiSiteCount; ++i)
    {
        SFermionBosonSource source;
        source.m_eSourceType = EFS_Point;
        source.m_byColorIndex = 0;
        source.m_sSourcePoint = __hostSiteIndexToInt4(i);
        pField->InitialAsSource(source);
        pField->D(_FIELDS);
        //pField->InverseD(appGetLattice()->m_pGaugeField);
        //pField->ApplyGamma(appGetLattice()->m_pGaugeField, GAMMA53);

        //pField->CopyTo(pFieldCopy);
        //CMeasureAngularMomentumKS::ApplyOrbitalMatrix(pField->m_pDeviceData, pFieldCopy->m_pDeviceData, pSU3Gauge->m_pDeviceData, pField->m_byFieldId);


        UINT size;
        BYTE* data = pField->CopyDataOut(size);
        Real* datar = (Real*)data;
        for (UINT j = 0; j < uiSiteCount; ++j)
        {
            toprint[i * uiSiteCount + j] = _make_cuComplex(datar[6 * j], datar[6 * j + 1]);
        }
    }

    pFieldCopy->Return();



    appPushLogDate(FALSE);
    appGeneral(_T("\nres=\n{\n"));
    for (UINT i = 0; i < uiSiteCount; ++i)
    {
        appGeneral(_T("\n{"));
        for (UINT uiSite = 0; uiSite < uiSiteCount - 1; ++uiSite)
        {
            if (appAbs(toprint[i * uiSiteCount + uiSite].x) < F(0.00001) && appAbs(toprint[i * uiSiteCount + uiSite].y) < F(0.00001))
            {
                appGeneral(_T("0,"));
            }
            else if (appAbs(toprint[i * uiSiteCount + uiSite].x) < F(0.00001))
            {
                appGeneral(_T("%fI,"), toprint[i * uiSiteCount + uiSite].y);
            }
            else if (appAbs(toprint[i * uiSiteCount + uiSite].y) < F(0.00001))
            {
                appGeneral(_T("%f,"), toprint[i * uiSiteCount + uiSite].x);
            }
            else
            {
                appGeneral(_T("%f%s%fI,"),
                    toprint[i * uiSiteCount + uiSite].x,
                    toprint[i * uiSiteCount + uiSite].y > F(0.0) ? _T("+") : _T("-"),
                    appAbs(toprint[i * uiSiteCount + uiSite].y)
                );
            }
        }

        if (i + 1 == uiSiteCount)
        {
            if (appAbs(toprint[(i + 1) * uiSiteCount - 1].x) < F(0.00001) && appAbs(toprint[(i + 1) * uiSiteCount - 1].y) < F(0.00001))
            {
                appGeneral(_T("0}\n"));
            }
            else if (appAbs(toprint[(i + 1) * uiSiteCount - 1].x) < F(0.00001))
            {
                appGeneral(_T("%fI}\n"), toprint[(i + 1) * uiSiteCount - 1].y);
            }
            else if (appAbs(toprint[(i + 1) * uiSiteCount - 1].y) < F(0.00001))
            {
                appGeneral(_T("%f}\n"), toprint[(i + 1) * uiSiteCount - 1].x);
            }
            else
            {
                appGeneral(_T("%f%s%fI}\n"),
                    toprint[(i + 1) * uiSiteCount - 1].x,
                    toprint[(i + 1) * uiSiteCount - 1].y > F(0.0) ? _T("+") : _T("-"),
                    appAbs(toprint[(i + 1) * uiSiteCount - 1].y)
                );
            }
            appGeneral(_T("};\n"));
        }
        else
        {
            if (appAbs(toprint[(i + 1) * uiSiteCount - 1].x) < F(0.00001) && appAbs(toprint[(i + 1) * uiSiteCount - 1].y) < F(0.00001))
            {
                appGeneral(_T("0},\n"));
            }
            else if (appAbs(toprint[(i + 1) * uiSiteCount - 1].x) < F(0.00001))
            {
                appGeneral(_T("%fI},\n"), toprint[(i + 1) * uiSiteCount - 1].y);
            }
            else if (appAbs(toprint[(i + 1) * uiSiteCount - 1].y) < F(0.00001))
            {
                appGeneral(_T("%f},\n"), toprint[(i + 1) * uiSiteCount - 1].x);
            }
            else
            {
                appGeneral(_T("%f%s%fI},\n"),
                    toprint[(i + 1) * uiSiteCount - 1].x,
                    toprint[(i + 1) * uiSiteCount - 1].y > F(0.0) ? _T("+") : _T("-"),
                    appAbs(toprint[(i + 1) * uiSiteCount - 1].y)
                );
            }
        }
    }
    appPopLogDate();

    //print it


    return 0;
}

//__REGIST_TEST(TestFermionDMatrixKS, Misc, TestDMatrixKS);


UINT TestFermionMatrixKS(CParameters& param)
{
    CFieldFermionKSSU3* pFermion = static_cast<CFieldFermionKSSU3*>(appGetLattice()->GetFieldById(2));
    //CFieldGaugeSU3* pGauge = static_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);
    CCString sFile = _T("KS_");
    param.FetchStringValue(_T("SaveHead"), sFile);
    CCString sFile2 = _T(".csv");

    for (INT i = 0; i < static_cast<INT>(EMD_Max); ++i)
    {
        CCString enumName = __ENUM_TO_STRING(EMeasureDiagnal, static_cast<EMeasureDiagnal>(i));
        appGeneral(_T("Working on %s ... \n"), enumName.c_str());
        ExportDiagnalStaggeredSU3(sFile + enumName + sFile2, static_cast<EMeasureDiagnal>(i), _FIELDS, pFermion);
    }

    return 0;
}

__REGIST_TEST(TestFermionMatrixKS, Tools, TestFermionMatrixKS, ExportFermionMatrix);

//UINT TestSmearingGauge(CParameters& param)
//{
//    //CIndexData::DebugStapleTable(1);
//
//    CFieldGauge* pGaugeCopy1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(appGetLattice()->m_pGaugeField[0]));
//    CFieldGauge* pGaugeCopy2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(appGetLattice()->m_pGaugeField[0]));
//    CGaugeSmearing* pGaugeSmearing = appGetGaugeSmearing(pGaugeCopy1->m_byFieldId);
//
//    pGaugeCopy1->ScalarMultply(F(0.25));
//    INT allpathes[24][3] = {
//        {2, 1, -2},
//        {-2, 1, 2},
//        {3, 1, -3},
//        {-3, 1, 3},
//        {4, 1, -4},
//        {-4, 1, 4},
//
//        {1, 2, -1},
//        {-1, 2, 1},
//        {3, 2, -3},
//        {-3, 2, 3},
//        {4, 2, -4},
//        {-4, 2, 4},
//
//        {2, 3, -2},
//        {-2, 3, 2},
//        {1, 3, -1},
//        {-1, 3, 1},
//        {4, 3, -4},
//        {-4, 3, 4},
//
//        {2, 4, -2},
//        {-2, 4, 2},
//        {3, 4, -3},
//        {-3, 4, 3},
//        {1, 4, -1},
//        {-1, 4, 1}
//    };
//    INT* pathbuffer = NULL;
//    appSimpleMalloc((void**)&pathbuffer, sizeof(INT) * 7);
//    //INT* start = allpathes[0];
//    //appSimpleCopyHD(pathbuffer, start, sizeof(INT) * 3);
//    //pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 3, 0, F(0.0625));
//
//    for (BYTE mu = 0; mu < 4; ++mu)
//    {
//        for (BYTE p = 0; p < 6; ++p)
//        {
//            BYTE pid = mu * 6 + p;
//            INT* start = allpathes[pid];
//            appSimpleCopyHD(pathbuffer, start, sizeof(INT) * 3);
//            appGeneral(_T("path: %d %d %d\n"), start[0], start[1], start[2]);
//            pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 3, mu, F(0.0625));
//        }
//    }
//
//
//    for (BYTE dir = 0; dir < 4; ++dir)
//    {
//        for (INT mu = 0; mu < 4; ++mu)
//        {
//            if (mu == dir)
//            {
//                continue;
//            }
//
//            for (INT nu = 0; nu < 4; ++nu)
//            {
//                if (nu == dir || nu == mu)
//                {
//                    continue;
//                }
//
//                INT forwardforward[5] = { mu + 1, nu + 1, dir + 1, -nu - 1, -mu - 1 };
//                INT forwardbackward[5] = { mu + 1, -nu - 1, dir + 1, nu + 1, -mu - 1 };
//                INT backwardforward[5] = { -mu - 1, nu + 1, dir + 1, -nu - 1, mu + 1 };
//                INT backwardbackward[5] = { -mu - 1, -nu - 1, dir + 1, nu + 1, mu + 1 };
//
//                appSimpleCopyHD(pathbuffer, forwardforward, sizeof(INT) * 5);
//                pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 5, dir, F(0.015625));
//                appSimpleCopyHD(pathbuffer, forwardbackward, sizeof(INT) * 5);
//                pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 5, dir, F(0.015625));
//                appSimpleCopyHD(pathbuffer, backwardforward, sizeof(INT) * 5);
//                pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 5, dir, F(0.015625));
//                appSimpleCopyHD(pathbuffer, backwardbackward, sizeof(INT) * 5);
//                pGaugeCopy2->AddLinkTo(pGaugeCopy1, pathbuffer, 5, dir, F(0.015625));
//            }
//        }
//    }
//
//    appSimpleFree(pathbuffer);
//    pGaugeSmearing->GaugeSmearing(pGaugeCopy2, NULL);
//
//    pGaugeCopy1->AxpyMinus(pGaugeCopy2);
//    DOUBLE res = pGaugeCopy1->Dot(pGaugeCopy1).x;
//
//    appGeneral(_T("\n\n========= res = %f ===============\n\n"), res);
//    return 0;
//}
//
//UINT TestSmearing(CParameters& param)
//{
//    CFieldGauge* pGauge[1] = { dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(appGetLattice()->m_pGaugeField[0]))};
//    CGaugeSmearing* pGaugeSmearing = appGetGaugeSmearing(pGauge[0]->m_byFieldId);
//    CFieldFermionKS* pFermion = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetFieldById(2));
//
//    CFieldFermionKS* pFermionCopy1 = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledCopy(pFermion));
//    CFieldFermionKS* pFermionCopy2 = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledCopy(pFermion));
//
//    //fermion copy use link
//    pFermionCopy1->D0(1, 0, pGauge, NULL);
//    pFermionCopy1->ScalarMultply(F(0.25));
//
//    SCHAR allpathes[24][3] = {
//    {2, 1, -2},
//    {-2, 1, 2},
//    {3, 1, -3},
//    {-3, 1, 3},
//    {4, 1, -4},
//    {-4, 1, 4},
//
//    {1, 2, -1},
//    {-1, 2, 1},
//    {3, 2, -3},
//    {-3, 2, 3},
//    {4, 2, -4},
//    {-4, 2, 4},
//
//    {2, 3, -2},
//    {-2, 3, 2},
//    {1, 3, -1},
//    {-1, 3, 1},
//    {4, 3, -4},
//    {-4, 3, 4},
//
//    {2, 4, -2},
//    {-2, 4, 2},
//    {3, 4, -3},
//    {-3, 4, 3},
//    {1, 4, -1},
//    {-1, 4, 1}
//    };
//    SCHAR* pathbuffer = NULL;
//    appSimpleMalloc((void**)&pathbuffer, sizeof(SCHAR) * 7);
//    //TArray<CField*> toreturn = f->CalculateRationalFields(m_pEffectiveGauge);
//    for (BYTE mu = 0; mu < 4; ++mu)
//    {
//        for (BYTE p = 0; p < 6; ++p)
//        {
//            BYTE pid = mu * 6 + p;
//            SCHAR* start = allpathes[pid];
//            appSimpleCopyHD(pathbuffer, start, sizeof(SCHAR) * 3);
//            pFermionCopy2->OneLinkS(pGauge[0]->GetData(), pGauge[0]->m_byFieldId, pFermionCopy1->GetData(), F(0.0625), pathbuffer, 3, mu, FALSE, EOCT_None, F(1.0), _onec);
//        }
//    }
//
//    for (BYTE dir = 0; dir < 4; ++dir)
//    {
//        for (SCHAR mu = 0; mu < 4; ++mu)
//        {
//            if (mu == dir)
//            {
//                continue;
//            }
//
//            for (SCHAR nu = 0; nu < 4; ++nu)
//            {
//                if (nu == dir || nu == mu)
//                {
//                    continue;
//                }
//
//                SCHAR forwardforward[5] = { mu + 1, nu + 1, dir + 1, -nu - 1, -mu - 1 };
//                SCHAR forwardbackward[5] = { mu + 1, -nu - 1, dir + 1, nu + 1, -mu - 1 };
//                SCHAR backwardforward[5] = { -mu - 1, nu + 1, dir + 1, -nu - 1, mu + 1 };
//                SCHAR backwardbackward[5] = { -mu - 1, -nu - 1, dir + 1, nu + 1, mu + 1 };
//
//                appSimpleCopyHD(pathbuffer, forwardforward, sizeof(SCHAR) * 5);
//                pFermionCopy2->OneLinkS(pGauge[0]->GetData(), pGauge[0]->m_byFieldId, pFermionCopy1->GetData(), F(0.015625), pathbuffer, 5, dir, FALSE, EOCT_None, F(1.0), _onec);
//                appSimpleCopyHD(pathbuffer, forwardbackward, sizeof(SCHAR) * 5);
//                pFermionCopy2->OneLinkS(pGauge[0]->GetData(), pGauge[0]->m_byFieldId, pFermionCopy1->GetData(), F(0.015625), pathbuffer, 5, dir, FALSE, EOCT_None, F(1.0), _onec);
//                appSimpleCopyHD(pathbuffer, backwardforward, sizeof(SCHAR) * 5);
//                pFermionCopy2->OneLinkS(pGauge[0]->GetData(), pGauge[0]->m_byFieldId, pFermionCopy1->GetData(), F(0.015625), pathbuffer, 5, dir, FALSE, EOCT_None, F(1.0), _onec);
//                appSimpleCopyHD(pathbuffer, backwardbackward, sizeof(SCHAR) * 5);
//                pFermionCopy2->OneLinkS(pGauge[0]->GetData(), pGauge[0]->m_byFieldId, pFermionCopy1->GetData(), F(0.015625), pathbuffer, 5, dir, FALSE, EOCT_None, F(1.0), _onec);
//            }
//        }
//    }
//
//    appSimpleFree(pathbuffer);
//
//    pGaugeSmearing->GaugeSmearing(pGauge[0], NULL);
//    pFermionCopy2->D0(1, 0, pGauge, NULL);
//
//    pFermionCopy1->AxpyMinus(pFermionCopy2);
//    DOUBLE res = pFermionCopy1->Dot(pFermionCopy1).x;
//
//    appGeneral(_T("\n\n========= res = %f ===============\n\n"), res);
//
//    return 0;
//}
//
//___REGIST_TEST(TestSmearing, Tools, TestCompareGaugeSmearing, Smearing, _TEST_NOCHECK);
//___REGIST_TEST(TestSmearingGauge, Tools, TestCompareGaugeSmearing, SmearingGauge, _TEST_NOCHECK);


//=============================================================================
// END OF FILE
//=============================================================================
