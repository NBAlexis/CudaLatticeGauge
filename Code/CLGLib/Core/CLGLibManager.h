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
        : m_pCudaHelper(NULL)
        , m_pLatticeData(NULL)
    {
        //appInitCriticalSection();
        appInitialTracer(GENERAL);
    }
    ~CCLGLibManager()
    {
        //appUninitCriticalSection();
    }

    void InitialWithParameter(class CParameters& params);

    /**
    * Free all data
    */
    void Quit();

    class CCudaHelper* m_pCudaHelper;
    class CLatticeData* m_pLatticeData;

};

extern CLGAPI CCLGLibManager GCLGManager;

extern void CLGAPI appInitialCLG(const TCHAR* paramFileName);

extern void CLGAPI appInitialCLG(class CParameters& params);

extern void CLGAPI appQuitCLG();

inline class CCudaHelper* appGetCudaHelper()
{
    return GCLGManager.m_pCudaHelper;
}

inline class CLatticeData* appGetLattice()
{
    return GCLGManager.m_pLatticeData;
}

__END_NAMESPACE

#endif //#ifndef _CLGLIBMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================