//=============================================================================
// FILENAME : Measure.cpp
//
// DESCRIPTION:
// Measurement driver for the PSU(3) fundamental-lift action with a Z3 2-form
// boundary field B.
//
// For every configuration N in [StartN, EndN]:
//   - loads <prefix>_<N>.con (gauge)
//   - loads the dynamic tensor2 companion <prefix>_<N>_t<fieldId>.con (B)
//     via CUpdator::LoadTensor2Configuration
//   - measures:
//       * Polyakov loop (CMeasurePolyakovXY, measure id 1)
//       * the PSU(3) action energy (needs the boundary field B)
//       * the Z3 monopole count of B (CFieldTensor2Z3::CountMonopoles)
//
// REVISION:
//  [08/11/26]
//=============================================================================
#include "PSU3WithBoundary.h"

INT Measure(CParameters& params)
{
    appSetupLog(params);

#pragma region read parameters
    INT iVaule = 0;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    CCString sSavePrefix = _T("PSU3");
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);

    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadFileType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }
#pragma endregion

    if (!appInitialCLG(params))
    {
        appCrucial(_T("PSU3WithBoundary Measure: Initial Failed!\n"));
        return 1;
    }

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));

    if (NULL == pGauge)
    {
        appCrucial(_T("PSU3WithBoundary Measure: gauge field (id 1) missing\n"));
        appQuitCLG();
        return 1;
    }

    UINT uiMeasured = 0;
    for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
    {
        CCString sFileName;
        sFileName.Format(_T("%s_%d.con"), sSavePrefix.c_str(), uiN);
        if (!CFileSystem::IsFileExist(sFileName))
        {
            appGeneral(_T("PSU3WithBoundary Measure: %s not found, skip\n"), sFileName.c_str());
            continue;
        }

        // load gauge configuration
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);
        // load the dynamic tensor2 companion (B field)
        appGetLattice()->m_pUpdator->LoadTensor2Configuration(sSavePrefix, uiN, eLoadType);
        appSynchronize();

        CCString sOut;
        sOut.Format(_T("PSU3WithBoundary Measure N=%u:"), uiN);

        // Polyakov loop
        if (NULL != pPL)
        {
            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            appSynchronize();
            if (pPL->m_lstLoop.Num() > 0)
            {
                CCString sItem;
                sItem.Format(_T(" |Polyakov|=%f"), cuCabs(pPL->m_lstLoop[0]));
                sOut = sOut + sItem;
            }
        }

        // PSU(3) action energy (needs B)
        if (NULL != pAction && NULL != pB)
        {
            const CFieldGauge* g[1] = { pGauge };
            const CFieldTensor2* t[1] = { pB };
            const DOUBLE fEnergy = pAction->Energy(FALSE, 1, 0, 1, g, NULL, t, NULL);
            CCString sItem;
            sItem.Format(_T(" |S|=%f"), fEnergy);
            sOut = sOut + sItem;
        }

        // Z3 monopole count of B
        if (NULL != pB)
        {
            const UINT uiMonopoles = pB->CountMonopoles();
            CCString sItem;
            sItem.Format(_T(" monopole=%u"), uiMonopoles);
            sOut = sOut + sItem;
        }

        appGeneral(_T("%s\n"), sOut.c_str());
        ++uiMeasured;
    }

    appGeneral(_T("PSU3WithBoundary Measure: %u configurations measured\n"), uiMeasured);
    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
