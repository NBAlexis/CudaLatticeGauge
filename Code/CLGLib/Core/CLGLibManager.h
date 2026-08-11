//=============================================================================
// FILENAME : CLGLibMananger.h
// 
// DESCRIPTION:
// This is the class for global start-up, control, shut-down
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CLGLIBMANAGER_H_
#define _CLGLIBMANAGER_H_

#define __Divisible(a, b) ( (b * (a/b)) == a )

__BEGIN_NAMESPACE

extern UINT CLGAPI appVersion();

struct CLGAPI SCLGLibManangerInitialCache
{
    ERandom eR;
    UINT constIntegers[kContentLength];
    INT constSignedIntegers[kContentLength];
    Real constFloats[kContentLength];
};

class CLGAPI CCLGLibManager
{
public:
    CCLGLibManager()
        : m_pCudaHelper(NULL)
        , m_pLatticeData(NULL)
        , m_pFileSystem(NULL)
        , m_pBuffer(NULL)
        , m_pFieldPool(NULL)
        , m_pComm(NULL)
        , m_pHaloManager(NULL)
        , m_ullLayoutGeneration(0)
        , m_iDeviceId(0)
        , m_InitialCache()
        , m_byLoadingFieldId(0)
    {
    }

    ~CCLGLibManager()
    {

    }

    UBOOL InitialWithParameter(class CParameters& params);

    /**
    * Free all data
    */
    void Quit();

    class CCudaHelper* m_pCudaHelper;
    class CLatticeData* m_pLatticeData;
    class CFileSystem* m_pFileSystem;
    class CCudaBuffer* m_pBuffer;
    class CFieldPool* m_pFieldPool;
    //Multi-GPU: allocated in every build, but a lone rank with a [1,1,1,1] grid
    //when _CLG_MULTI_GPU is 0, so the single-GPU path is unchanged.
    class CLGComm* m_pComm;
    class CHaloManager* m_pHaloManager;

    //Improve-1 (multi-GPU-improve1.md 3.1): bumped by InitialLatticeAndConstant
    //every time the lattice/process-grid/halo-width constants are (re)baked.
    //CHaloBufferHandle Binds snapshot it; a stale generation means the handle's
    //halo can never be reused across a layout rebuild.
    ULONGLONG m_ullLayoutGeneration;
    TArray<class CRegisteredBufferCache*> m_lstBufferCaches;

    void RegisterCache(class CRegisteredBufferCache* cacher)
    {
        m_lstBufferCaches.AddItem(cacher);
    }
    
    void SetupLog(class CParameters& params);

    INT m_iDeviceId;

protected:

    SCLGLibManangerInitialCache m_InitialCache;

    //==================================
    //Create

    //Required
    void InitialLatticeAndConstant(class CParameters& params);
    void InitialRandom(class CParameters& params);

    //Optional
    class CField* CreateGaugeFields(class CParameters& params) const;
    class CField* CreateFermionFields(class CParameters& params) const;
    class CField* CreateBosonFields(class CParameters& params) const;
    class CField* CreateTensor2Fields(class CParameters& params) const;

    //other gauge fields is to be removed...
    //void CreateOtherGaugeFields(class CParameters& params) const;

    void CreateBoundaryFields(class CParameters& params, const CCString& sDefaultName) const;
    void CreateIndexAndBoundary(class CParameters& params) const;
    void CreateActionList(class CParameters& params);
    void CreateUpdator(class CParameters& params) const;
    void CreateMeasurement(class CParameters& params);
    void CreateSolver(class CParameters& params) const;
    void CreateMultiShiftSolver(class CParameters& params) const;
    void CreateGaugeSmearing(class CParameters& params) const;
    void CreateGaugeFixing(class CParameters& params) const;
    void CreateGaugeStapleCache(class CParameters& params) const;

    //==================================
    //Cache

    //Requared
    void InitialIndexBuffer() const;

    void InitialFieldBuffer() const;

    BYTE m_byLoadingFieldId = 0;
};

extern CLGAPI CCLGLibManager GCLGManager;

inline void CLGAPI appSetupLog(class CParameters& params)
{
    GCLGManager.SetupLog(params);
}

extern UBOOL CLGAPI appInitialCLG(const TCHAR* paramFileName);
extern UBOOL CLGAPI appInitialCLG(class CParameters& params);

extern void CLGAPI appQuitCLG();
extern void CLGAPI appFailQuitCLG();

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

inline class CSLASolver* appGetFermionSolver(BYTE byFieldId)
{
    return appGetLattice()->m_pFermionSolver[byFieldId];
}

inline class CMultiShiftSolver* appGetMultiShiftSolver(BYTE byFieldId)
{
    return appGetLattice()->m_pFermionMultiShiftSolver[byFieldId];
}

inline class CGaugeSmearing* appGetGaugeSmearing(BYTE byFieldId)
{
    return appGetLattice()->m_pGaugeSmearing[byFieldId];
}

inline class CStapleCache* appGetStapleCache(BYTE byFieldId)
{
    return appGetLattice()->m_pStapleCaches[byFieldId];
}

inline class CFieldPool* appGetFieldPool()
{
    return GCLGManager.m_pFieldPool;
}

inline CCudaBuffer* GetBuffer()
{
    return GCLGManager.m_pBuffer;
}

inline class CLGComm* appGetComm()
{
    return GCLGManager.m_pComm;
}

inline class CHaloManager* appGetHaloManager()
{
    return GCLGManager.m_pHaloManager;
}

//Improve-1: current layout generation (see CCLGLibManager::m_ullLayoutGeneration).
inline ULONGLONG appGetLayoutGeneration()
{
    return GCLGManager.m_ullLayoutGeneration;
}

inline INT appGetDeviceId()
{
    return GCLGManager.m_iDeviceId;
}

/**
* find all factors of input number
*/
inline TArray<UINT> _getFactors(UINT length)
{
    TArray<UINT> ret;
    ret.AddItem(1);
    for (UINT i = 2; i < (length / 2); ++i)
    {
        if (__Divisible(length, i))
        {
            ret.AddItem(i);
        }
    }
    if (length > 1)
    {
        ret.AddItem(length);
    }
    return ret;
}

/**
* find the max block size for thread decompose
*/
inline TArray<UINT> _getDecompose(const TArray<UINT>& contraints, const TArray<UINT>& latticeLength)
{
    UINT uiBlockSize = 1;
    TArray<UINT> ret;

    //number of blocks
    ret.AddItem(latticeLength[0]);
    ret.AddItem(latticeLength[1]);
    ret.AddItem(latticeLength[2]);
    //block size
    ret.AddItem(1);
    ret.AddItem(1);
    ret.AddItem(1);

    TArray<UINT> factorsOfX = _getFactors(latticeLength[0]);
    TArray<UINT> factorsOfY = _getFactors(latticeLength[1]);
    TArray<UINT> factorsOfZ = _getFactors(latticeLength[2]);
    for (INT i = 0; i < factorsOfX.Num(); ++i)
    {
        for (INT j = 0; j < factorsOfY.Num(); ++j)
        {
            for (INT k = 0; k < factorsOfZ.Num(); ++k)
            {
                if (factorsOfX[i] <= (UINT)contraints[1]
                    && factorsOfY[j] <= (UINT)contraints[2]
                    && factorsOfZ[k] <= (UINT)contraints[3])
                {
                    const UINT uiThreadPerBlcok = factorsOfX[i] * factorsOfY[j] * factorsOfZ[k];
                    if (uiThreadPerBlcok <= (UINT)contraints[0]
                        && uiThreadPerBlcok > uiBlockSize)
                    {
                        uiBlockSize = uiThreadPerBlcok;

                        //number of blocks
                        ret[0] = latticeLength[0] / factorsOfX[i];
                        ret[1] = latticeLength[1] / factorsOfY[j];
                        ret[2] = latticeLength[2] / factorsOfZ[k];

                        //block size
                        ret[3] = factorsOfX[i];
                        ret[4] = factorsOfY[j];
                        ret[5] = factorsOfZ[k];
                    }
                }
            }
        }
    }

    return ret;
}


inline class CCString GetCLGVersion()
{
    CCString sRet;
    sRet.Format(_T("%d.%d"), __GVERSION, __GVERSION_S);
    return sRet;
}

__END_NAMESPACE

#endif //#ifndef _CLGLIBMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================