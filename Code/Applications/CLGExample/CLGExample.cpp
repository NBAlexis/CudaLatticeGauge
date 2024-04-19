//=============================================================================
// FILENAME : CLGExample.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/18/2024 nbale]
//=============================================================================

#include "CLGLib.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("CLGExample.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/CLGExample.yaml"), params);
#endif
    appSetupLog(params);
    appInitialCLG(params);

    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("CLGExample"));
    appGetLattice()->m_pUpdator->UpdateUntileAccept(10, FALSE);

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
