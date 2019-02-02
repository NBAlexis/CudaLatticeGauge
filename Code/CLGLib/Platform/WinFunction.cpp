//=============================================================================
// FILENAME : WinFunction.cpp
// 
// DESCRIPTION:
// This is the system functions for MS-VC platform
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CLGAPI INT appUnicodeToAnsi(ANSICHAR* mbstr, const UNICHAR* wcstr, INT bufsize)
{
    if (bufsize == 0 && mbstr != NULL)
        return 0;

    INT result = ::WideCharToMultiByte(CP_ACP, 0, wcstr, -1,
        mbstr, bufsize, NULL, NULL);
    assert(NULL == mbstr || result <= static_cast<INT>(bufsize));
    if (result > 0)
        mbstr[result - 1] = 0;
    return result;
}

CLGAPI INT appAnsiToUnicode(UNICHAR* wcstr, const ANSICHAR* mbstr, INT bufsize)
{
    INT count = bufsize / 2;
    if (count == 0 && wcstr != NULL)
        return 0;

    INT result = ::MultiByteToWideChar(CP_ACP, 0, mbstr, -1,
        wcstr, count);
    assert(wcstr == NULL || result <= static_cast<INT>(count));
    if (result > 0)
        wcstr[result - 1] = 0;
    else
    {
        //show_last_error();
    }
    return result * 2;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
