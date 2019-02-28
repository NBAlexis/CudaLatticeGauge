//=============================================================================
// FILENAME : CFile.h
// 
// DESCRIPTION:
// 
// REVISION:
//  [01/31/2019 nbale]
//=============================================================================
#pragma once

#ifndef _CFILE_H_
#define _CFILE_H_

///TODO: Use boost 

__BEGIN_NAMESPACE

enum { DETAILED_STATS = 0 };
enum { kMaxPathNameLength = 1024 };
enum { kInvalidReturn = -1 };

class CLGAPI CFileSystem
{
public:
    void Init() {}

    /**
    * Need to free the pointer
    */
    BYTE* ReadAllBytes(const TCHAR* sFilename, UINT& size);
    UBOOL WriteAllBytes(const TCHAR* sFilename, BYTE* data, UINT uiSize);
    CCString ReadAllText(const TCHAR* sFilename);
    UBOOL WriteAllText(const TCHAR* sFilename, const CCString& data);

    //UBOOL MakeDir(const CCString& dirPath);

    UBOOL SetDefaultDirectory(const TCHAR* Filename);
    class CCString GetDefaultDirectory();
};

__END_NAMESPACE

#endif //#ifndef _CFILE_H_

//=============================================================================
// END OF FILE
//=============================================================================

