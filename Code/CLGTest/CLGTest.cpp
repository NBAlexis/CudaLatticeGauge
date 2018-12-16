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
    CTimer timer;
    timer.Start();

    appInitialCLG(_T("TestSuit.yaml"));

    timer.Stop();
    appGeneral(_T("\n =========== CLG Lib Initialed, cost: %f(ms)\n\n"), timer.Elapsed());
    timer.Reset();
    

    //CCommonData::InitialWithDefault();
    //
    //CLatticeData* lattice = CLatticeData::GetInstance();
    //CFieldGaugeSU3* field = new CFieldGaugeSU3(lattice);

    //CParameters params;
    //CYAMLParser::ParseFile(_T("TestSuit.yaml"), params);

    //params.Dump();

    //params.GetParameter(_T("Quark_2")).SetStringVaule(_T("__abc__"), _T("__123__"));

    //params.Dump();
    //THashMap<CCString, CCString> testhash;
    //testhash[CCString(_T("123"))] = _T("456");
    //testhash.SetAt(_T("1231"), _T("4567"));
    //TArray<CCString> keys = testhash.GetAllKeys();
    //appGeneral(_T("t1:%d %s->%s\n"), testhash.Exist(keys[0]), keys[0], testhash[keys[0]]);
    //appGeneral(_T("t2:%d %s->%s"), testhash.Exist(keys[1]), keys[1], testhash[keys[1]]);
    //THashMap<CCString, INT> testhash = appGetEnumTable(_T("    Ex = 0,\n        Ey,\n        Ec = 0x77, "));
    //TArray<CCString> keys = testhash.GetAllKeys();
    //appGeneral(_T("t1:%d %s->\n"), testhash.Exist(keys[0]), keys[0]);
    //appGeneral(_T("t1:%d %s->%d\n"), testhash.Exist(keys[0]), keys[0], testhash[keys[0]]);
    //appGeneral(_T("t2:%d %s->\n"), testhash.Exist(keys[1]), keys[1]);
    //appGeneral(_T("t2:%d %s->%d"), testhash.Exist(keys[1]), keys[1], testhash[keys[1]]);

    //CCString a1 = _T("Ex");
    //CCString a2 = _T("Ex");
    //BYTE dwa1 = *((BYTE*)(&a1));
    //BYTE dwa2 = *((BYTE*)(&a2));
    //appGeneral(_T("test hash map key %d %d\n"), dwa1, dwa2);
    //appGeneral(_T("test hash map key %d %d"), TMapHashKey<const CCString&>(a1), TMapHashKey<const CCString&>(a2));

    //appGeneral(_T("\ntest enum %d\n"), appStringToEnumEabc(_T("Ex")));
    //
    //appGeneral(_T("\ntest enum %s\n"), __ENUM_TO_STRING(EFieldType, EFT_GaugeSU3).c_str());

    //appGeneral(_T("\ntest enum %d\n"), __STRING_TO_ENUM(EFieldType, _T("EFT_Max")));

    //appGeneral(_T("\nPress any key to continue\n"));
    //CIN.get();

    //appSafeDelete(field);
    //CLatticeData::Release();

    //appGetLattice()->m_pGaugeField->DebugPrintMe();

    UINT iErrors = 0;
    timer.Start();
    iErrors += TestRandom();

    
    //test update
    appGetLattice()->m_pUpdator->Update(10);

    timer.Stop();
    appGeneral(_T("\n =========== CLG Test Finished, errors: %d, cost: %f(ms)\n"), iErrors, timer.Elapsed());

    appGeneral(_T("\nPress any key to quit\n"));
    CIN.get();

    appQuitCLG();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
