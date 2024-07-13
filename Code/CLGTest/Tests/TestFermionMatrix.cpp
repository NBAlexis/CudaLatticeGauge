//=============================================================================
// FILENAME : TestFermionMatrix.cpp
// 
// DESCRIPTION:
//
//     Test the operations on fields
//
// REVISION:
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
    CFieldFermionKSSU3* pFieldCopy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pField->m_byFieldId));

    //CFieldGaugeSU3* pSU3Gauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);

    for (UINT i = 0; i < uiSiteCount; ++i)
    {
        SFermionBosonSource source;
        source.m_eSourceType = EFS_Point;
        source.m_byColorIndex = 0;
        source.m_sSourcePoint = __hostSiteIndexToInt4(i);
        pField->InitialAsSource(source);
        CCommonData::m_fOmega = F(0.2);
        pField->D(_FIELDS);
        //pField->InverseD(appGetLattice()->m_pGaugeField);
        //pField->ApplyGammaKS(appGetLattice()->m_pGaugeField, GAMMA53);

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

//=============================================================================
// END OF FILE
//=============================================================================
