//=============================================================================
// FILENAME : CUUID.h
// 
// DESCRIPTION:
//  Use random number
//
// REVISION:
//  [mm/dd/yy]
//  [4/26/2025 nbale]
//=============================================================================
#pragma once

#include <array>
#include <random>
//#include <sstream>
//#include <string>
//#include <iomanip>

#ifndef _CUUID_H_
#define _CUUID_H_

__BEGIN_NAMESPACE

inline CCString GetUUID()
{
    TArray<BYTE> data;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<UINT> dis(0, 255);

    for (UINT i = 0; i < 16; ++i)
    {
        data.AddItem(static_cast<BYTE>(dis(gen)));
    }

    data[6] = (data[6] & 0x0F) | 0x40;
    data[8] = (data[8] & 0x3F) | 0x80;
    const TCHAR* alphabetaTable = "0123456789ABCDEF";
    CCString ret = _T("");
    TCHAR strs[3];
    strs[2] = 0;
    for (BYTE i = 0; i < 16; ++i)
    {
        if (i == 4 || i == 6 || i == 8 || i == 10)
        {
            ret = ret + _T("-");
        }
        strs[0] = alphabetaTable[data[i] / 16];
        strs[1] = alphabetaTable[data[i] % 16];
        ret = ret + CCString(strs);
    }

    return ret;
}

__END_NAMESPACE

#endif//#ifndef _CUUID_H_

//=============================================================================
// END OF FILE
//=============================================================================