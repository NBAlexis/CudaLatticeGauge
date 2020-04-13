//=============================================================================
// FILENAME : ConfigurationCompresser.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/10/2020 nbale]
//=============================================================================

#include "ConfigurationCompresser.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("ConfigurationCompresser.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/ConfigurationCompresser.yaml"), params);
#endif

    appSetupLog(params);

    INT iVaule = 0;
    params.FetchValueINT(_T("StartIndex"), iVaule);
    const UINT iStartIndex = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("EndIndex"), iVaule);
    const UINT iEndIndex = static_cast<UINT>(iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("StartOmega"), iVaule);
    //const UINT iStartOmega = static_cast<UINT>(iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("EndOmega"), iVaule);
    //const UINT iEndOmega = static_cast<UINT>(iVaule);
    
    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());
    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);
    appGeneral(_T("load prefix: %s\n"), sLoadPrefix.c_str());

    CCString sJobType = _T("EJT_Matching");
    params.FetchStringValue(_T("JobType"), sJobType);
    EJobType eJob = __STRING_TO_ENUM(EJobType, sJobType);

    //=========================================================
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    if (EJT_Matching == eJob)
    {
        appSetLogDate(FALSE);
        for (UINT iIndex = iStartIndex; iIndex <= iEndIndex; ++iIndex)
        {
            CCString sFileNameLoad;
            CCString sFileNameSave;

            sFileNameLoad.Format(_T("%sMatching_%d.con"), sLoadPrefix, iIndex);
            sFileNameSave.Format(_T("%sMatching_%d.cco"), sSavePrefix, iIndex);

            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileNameLoad, EFFT_CLGBin);
            appGetLattice()->m_pGaugeField->SaveToCompressedFile(sFileNameSave);
            appGeneral(_T("=%s"), 0 == (iIndex % 50) ? _T("\n") : _T(""));
        }
        appSetLogDate(TRUE);
    }

    appGeneral(_T("\n================= Done ===============\n"));
    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
