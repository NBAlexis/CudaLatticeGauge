//=============================================================================
// FILENAME : GaugeFixing.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/21/2020 nbale]
//=============================================================================

#include "StaggeredSpectrum.h"

INT StaggeredGaugeFixing(CParameters& params)
{
    appSetupLog(params);
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    INT iVaule = 1;
    params.FetchValueINT(_T("IndexStart"), iVaule);
    const UINT iIndexStart = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("IndexEnd"), iVaule);
    const UINT iIndexEnd = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OnlyCheck"), iVaule);
    const UBOOL bOnlyCheck = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("CheckAndFix"), iVaule);
    const UBOOL bAlsoFix = 0 != iVaule;

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());


    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);
    appGeneral(_T("load prefix: %s\n"), sLoadPrefix.c_str());

    if (bOnlyCheck)
    {
        appGeneral(_T("====== Start: %d to %d ======\n"), iIndexStart, iIndexEnd);
        appSetLogDate(FALSE);
        for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
        {
            CCString sSaveFile;
            sSaveFile.Format(_T("%sMatching_%d.con"), sSavePrefix.c_str(), uiIndex);
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sSaveFile, EFFT_CLGBin);

#if !_CLG_DOUBLEFLOAT
            const DOUBLE fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
            if (fRes >= 0.0 && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy)
#else
            const Real fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
            if (fRes >= F(0.0) && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy)
#endif
            {
                if (0 == (uiIndex % 20))
                {
                    appGeneral(_T("%1.2f,\n"), fRes);
                }
                else
                {
                    appGeneral(_T("%1.2f,"), fRes);
                }

                continue;
            }

            if (bAlsoFix)
            {
                if (fRes >= F(0.0) && fRes < F(0.01))
                {
                    appGeneral(_T("\nNot good enough : %d \n"), uiIndex);
                    appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                    appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
                }
                else
                {
                    appGeneral(_T("\nBad : %d \n"), uiIndex);
                    CCString sLoadFile;
                    sLoadFile.Format(_T("%sMatching_%d.con"), sLoadPrefix.c_str(), uiIndex);
                    appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, EFFT_CLGBin);
                    appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                    appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
                }
            }
            else
            {
                appGeneral(_T("\nBad : %d, is %2.20f\n"), uiIndex, fRes);
            }
        }
        appSetLogDate(TRUE);
    }
    else
    {
        for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
        {
            CCString sLoadFile;
            CCString sSaveFile;
            sLoadFile.Format(_T("%sMatching_%d.con"), sLoadPrefix.c_str(), uiIndex);
            sSaveFile.Format(_T("%sMatching_%d.con"), sSavePrefix.c_str(), uiIndex);
            appGeneral(_T("Fixing : %d \n"), uiIndex);
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, EFFT_CLGBin);
            appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
            appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
        }
    }

    appGeneral(_T("\n=====================================\n========= finished! ==========\n"));
    appQuitCLG();

    return 0;
}
