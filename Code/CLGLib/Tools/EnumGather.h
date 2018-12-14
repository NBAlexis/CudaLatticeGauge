//=============================================================================
// FILENAME : EnumGather.h
// 
// DESCRIPTION:
// This is the macro tools to gather the enum for YAML file
// We assume this is only use for initialize, so we do not care the speed
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _ENUMGATHER_H_
#define _ENUMGATHER_H_

__BEGIN_NAMESPACE

/**
* Parse the string
* 
* "Ex = 0,
*  Ey,
*  Ec = 0x77,"
*
* To Hash map :
* "Ex"->0
* "Ey"->1
* "Ec"->0x77
*/
inline THashMap<CCString, INT> appGetEnumTable(const CCString &inGatheredEnum)
{
    TArray <INT> inSeps;
    inSeps.AddItem(_T(','));
    TArray <CCString> sEnumArray = appGetStringList(inGatheredEnum, inSeps, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
    INT iCurrentNum = -1;
    THashMap<CCString, INT> sOutTable;
    for (INT i = 0; i < sEnumArray.Num(); ++i)
    {
        TArray <INT> inSeps2;
        inSeps2.AddItem(_T('='));
        TArray <CCString> sEnumArray2 = appGetStringList(sEnumArray[i], inSeps2, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
        assert(1 == sEnumArray2.Num() || 2 == sEnumArray2.Num());
        if (2 == sEnumArray2.Num())
            iCurrentNum = appStrToINT(sEnumArray2[1]);
        else
            ++iCurrentNum;
        sOutTable.SetAt(sEnumArray2[0], iCurrentNum);
    }

    return sOutTable;
}

//This is an example
enum Eabc
{
    Ex = 0,
    Ey,
    Ec = 0x77,
};

inline CCString appEnumToStringEabc(Eabc v)
{
    TArray <INT> inSeps;
    inSeps.AddItem(_T(','));
    TArray <CCString> sEnumArray = appGetStringList(_T("    Ex = 0,\n        Ey,\n        Ec = 0x77, "), inSeps, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
    INT iCurrentNum = -1;
    for (INT i = 0; i < sEnumArray.Num(); ++i)
    {
        TArray <INT> inSeps2;
        inSeps2.AddItem(_T('='));
        TArray <CCString> sEnumArray2 = appGetStringList(sEnumArray[i], inSeps2, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
        assert(1 == sEnumArray2.Num() || 2 == sEnumArray2.Num());
        if (2 == sEnumArray2.Num())
        {
            iCurrentNum = appStrToINT(sEnumArray2[1]);
        }
        else
        {
            ++iCurrentNum;
        }
        if (iCurrentNum == (INT)v)
        {
            return sEnumArray2[0];
        }
    }
    assert(FALSE);
    return "";
}

inline Eabc appStringToEnumEabc(const CCString& s)
{
    THashMap<CCString, INT> theTable = appGetEnumTable(_T("    Ex = 0,\n        Ey,\n        Ec = 0x77, "));
    assert(theTable.Exist(s));
    return (Eabc)theTable[s];
}

__END_NAMESPACE

#define __DEFINE_ENUM_FUNCTION(enumname, ...) inline CCString appEnumToString##enumname(enumname v) \
{ \
    TArray <INT> inSeps; \
    inSeps.AddItem(_T(',')); \
    TArray <CCString> sEnumArray = appGetStringList(_T(#__VA_ARGS__), inSeps, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety); \
    INT iCurrentNum = -1; \
    for (INT i = 0; i < sEnumArray.Num(); ++i) \
    { \
        TArray <INT> inSeps2; \
        inSeps2.AddItem(_T('=')); \
        TArray <CCString> sEnumArray2 = appGetStringList(sEnumArray[i], inSeps2, EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety); \
        assert(1 == sEnumArray2.Num() || 2 == sEnumArray2.Num()); \
        if (2 == sEnumArray2.Num()) \
        { \
            iCurrentNum = appStrToINT(sEnumArray2[1]); \
        } \
        else \
        { \
            ++iCurrentNum; \
        } \
        if (iCurrentNum == (INT)v) \
        { \
            return sEnumArray2[0]; \
        } \
    } \
    assert(FALSE); \
    return ""; \
} \
inline enumname appStringToEnum##enumname(const CCString& s) \
{ \
    THashMap<CCString, INT> theTable = appGetEnumTable(_T(#__VA_ARGS__)); \
    assert(theTable.Exist(s)); \
    return (enumname)theTable[s]; \
}




#define __DEFINE_ENUM(enumname, ...) \
    enum enumname \
    { \
        __VA_ARGS__ \
    }; \
    __DEFINE_ENUM_FUNCTION(enumname, __VA_ARGS__);



#define __ENUM_TO_STRING(enumname, enumvalue) (appEnumToString##enumname(enumvalue))

#define __STRING_TO_ENUM(enumname, stringvalue) (appStringToEnum##enumname(stringvalue))


#endif //#ifndef _ENUMGATHER_H_


//=============================================================================
// END OF FILE
//=============================================================================
