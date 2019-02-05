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

struct CLGAPI SCLGLibManangerInitialCache
{
    ERandom eR;
    UINT constIntegers[kContentLength];
    Real constFloats[kContentLength];
};

class CLGAPI CCLGLibManager
{
public:
    CCLGLibManager()
        : m_pCudaHelper(NULL)
        , m_pLatticeData(NULL)
        , m_pFileSystem(NULL)
    {
        //appInitCriticalSection();
        appInitialTracer(GENERAL);
    }
    ~CCLGLibManager()
    {
        //appUninitCriticalSection();
    }

    UBOOL InitialWithParameter(class CParameters& params);

    /**
    * Free all data
    */
    void Quit();

    class CCudaHelper* m_pCudaHelper;
    class CLatticeData* m_pLatticeData;
    class CFileSystem* m_pFileSystem;

    void SetupLog(class CParameters& params);

protected:

    SCLGLibManangerInitialCache m_InitialCache;

    //Required
    void InitialLatticeAndConstant(class CParameters& params);
    void InitialRandom(class CParameters& params);

    //Optional
    void CreateGaugeField(class CParameters& params);
    void CreateFermionFields(class CParameters& params);
    void CreateIndexAndBoundary(class CParameters& params);
    void CreateActionList(class CParameters& params);
    void CreateUpdator(class CParameters& params);
    void CreateMeasurement(class CParameters& params);
    void CreateSolver(class CParameters& params);

};

extern CLGAPI CCLGLibManager GCLGManager;

inline void CLGAPI appSetupLog(class CParameters& params)
{
    GCLGManager.SetupLog(params);
}

extern UBOOL CLGAPI appInitialCLG(const TCHAR* paramFileName);
extern UBOOL CLGAPI appInitialCLG(class CParameters& params);

extern void CLGAPI appQuitCLG();

inline class CCudaHelper* appGetCudaHelper()
{
    return GCLGManager.m_pCudaHelper;
}

inline class CLatticeData* appGetLattice()
{
    return GCLGManager.m_pLatticeData;
}

inline class CFileSystem* appGetFileSystem()
{
    return GCLGManager.m_pFileSystem;
}

inline class CSLASolver* appGetFermionSolver()
{
    return appGetLattice()->m_pFermionSolver;
}

__END_NAMESPACE

#endif //#ifndef _CLGLIBMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================