//=============================================================================
// FILENAME : CLGLibMananger.h
// 
// DESCRIPTION:
// This is the class for global start-up, control, shut-down
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CLGLIBMANAGER_H_
#define _CLGLIBMANAGER_H_

__BEGIN_NAMESPACE

class CLGAPI CCLGLibManager
{
public:
    CCLGLibManager()
    {
        appInitCriticalSection();
    }
    ~CCLGLibManager()
    {
        appUninitCriticalSection();
    }
};

extern CLGAPI CCLGLibManager GCLGManager;

__END_NAMESPACE

#endif //#ifndef _CLGLIBMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================