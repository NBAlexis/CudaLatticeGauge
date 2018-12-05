//=============================================================================
// FILENAME : CLGTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGTest.h"

int main(int argc, char * argv[])
{
    appInitialTracer(GENERAL);

    CCommonData::InitialWithDefault();
    
    CLatticeData* lattice = new CLatticeData();
    CFieldGaugeSU3* field = new CFieldGaugeSU3(lattice);

    appGeneral(_T("\nPress any key to continue\n"));
    CIN.get();

    appSafeDelete(field);
    appSafeDelete(lattice);

    appGeneral(_T("\nPress any key to quit\n"));
    CIN.get();

    

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
