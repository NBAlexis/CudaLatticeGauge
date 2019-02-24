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

#define __Divisible(a, b) ( (b * (a/b)) == a )

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
        , m_pBuffer(NULL)
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
    void CreateGaugeSmearing(class CParameters& params);

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

inline class CGaugeSmearing* appGetGaugeSmearing()
{
    return appGetLattice()->m_pGaugeSmearing;
}

inline CCudaBuffer* GetBuffer()
{
    return GCLGManager.m_pBuffer;
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
                    UINT uiThreadPerBlcok = factorsOfX[i] * factorsOfY[j] * factorsOfZ[k];
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

__host__ __inline__ static UINT _hostGetSiteIndex(const SSmallInt4& coord)
{
    return static_cast<UINT>(coord.x * _HC_MultX + coord.y * _HC_MultY + coord.z * _HC_MultZ + coord.w);
}

__END_NAMESPACE

#endif //#ifndef _CLGLIBMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================