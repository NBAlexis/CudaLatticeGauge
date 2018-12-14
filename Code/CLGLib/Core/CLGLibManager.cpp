//=============================================================================
// FILENAME : CLGLibMananger.cpp
// 
// DESCRIPTION:
// This is the class for global start-up, control, shut-down
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

#define __CheckTag(tagname,...) iTag = 0; if (params.FetchValueINT(tagname, iTag) && (0 != iTag)) {__VA_ARGS__;}

#define __FetchIntWithDefaultSub(paramname, tagname, defaultv) if (!paramname.FetchValueINT(tagname, iVaules)) { iVaules = defaultv; }

#define __FetchIntWithDefault(tagname, defaultv) __FetchIntWithDefaultSub(params, tagname, defaultv)

#define __FetchStringWithDefaultSub(paramname, tagname, defaultv) if (!paramname.FetchStringValue(tagname, sValues)) { sValues = defaultv; }

#define __FetchStringWithDefault(tagname, defaultv) __FetchStringWithDefaultSub(params, tagname, defaultv)

#define __Divisible(a, b) ( (b * (a/b)) == a )

__BEGIN_NAMESPACE

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
        for (INT j = 0; j < factorsOfX.Num(); ++j)
        {
            for (INT k = 0; k < factorsOfX.Num(); ++k)
            {
                if (factorsOfX[i] <= (UINT)contraints[1]
                 && factorsOfY[j] <= (UINT)contraints[2]
                 && factorsOfZ[k] <= (UINT)contraints[3])
                {
                    UINT uiThreadPerBlcok = factorsOfX[i] * factorsOfX[j] * factorsOfX[k];
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

CLGAPI CCLGLibManager GCLGManager;

void CCLGLibManager::InitialWithParameter(CParameters &params)
{
    //GClassGather.TraceAllClass();
    checkCudaErrors(cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1 << 30));

#pragma region Verbose

    //Setup outputs
    CCString verboselevel;
    EVerboseLevel eVerbLevel = CRUCIAL;
    CCString sVerbFile = _T("stdout");
    UBOOL fetchVerbLevel = params.FetchStringValue(_T("VerboseLevel"), verboselevel);
    UBOOL fetchVerbFile = params.FetchStringValue(_T("VerboseOutput"), sVerbFile);
    if (fetchVerbLevel || fetchVerbFile) //do NOT put fetch string in if, it will enter if when the first is TRUE
    {
        eVerbLevel = __STRING_TO_ENUM(EVerboseLevel, verboselevel);
        appSetTracer(eVerbLevel, sVerbFile);
    }

    //check whether to log parameter file
    INT iTag = 0;
    appGeneral(_T("============================== Parameter =============================\n\n"));
    __CheckTag(_T("ShowParameterContent"), params.Dump());
    appGeneral(_T("============================== GPU =============================\n\n"));
    __CheckTag(_T("ShowDeviceInformation"), CCudaHelper::DeviceQuery());

    appGeneral(_T("============================== Log Start =============================\n\n"));

#pragma endregion

    UINT constIntegers[kContentLength];
    FLOAT constFloats[kContentLength];

    INT iVaules = 0;

#pragma region Lattice Size and Threads

    __FetchIntWithDefault(_T("Dim"), 4);
    assert(iVaules > 1 && iVaules < 5);
    constIntegers[ECI_Dim] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Dir"), 4);
    assert(iVaules > 1);
    constIntegers[ECI_Dir] = static_cast<UINT>(iVaules);

    TArray<INT> intValues;
    if (!params.FetchValueArrayINT(_T("LatticeLength"), intValues))
    {
        appCrucial(_T("LatticeLength not found, will use 8x8x8x8"));
        intValues.RemoveAll();
        intValues.AddItem(8);
        intValues.AddItem(8);
        intValues.AddItem(8);
        intValues.AddItem(8);
    }
    else if (intValues[0] < 1
          || intValues[1] < 1
          || intValues[2] < 1
          || intValues[3] < 1)
    {
        appCrucial(_T("Lattice length is invalid, will use 8x8x8x8"));
        intValues.RemoveAll();
        intValues.AddItem(8);
        intValues.AddItem(8);
        intValues.AddItem(8);
        intValues.AddItem(8);
    }

    constIntegers[ECI_Lx] = static_cast<UINT>(intValues[0]);
    constIntegers[ECI_Ly] = static_cast<UINT>(intValues[1]);
    constIntegers[ECI_Lz] = static_cast<UINT>(intValues[2]);
    constIntegers[ECI_Lt] = static_cast<UINT>(intValues[3]);
    constIntegers[ECI_Volumn] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2] * intValues[3]);
    constIntegers[ECI_MultX] = static_cast<UINT>(intValues[1] * intValues[2] * intValues[3]);
    constIntegers[ECI_MultY] = static_cast<UINT>(intValues[2] * intValues[3]);
    constIntegers[ECI_MultZ] = static_cast<UINT>(intValues[3]);

    UBOOL bAutoDecompose = TRUE;
    __FetchIntWithDefault(_T("ThreadAutoDecompose"), 0);

    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();

    if (0 == iVaules)
    {
        bAutoDecompose = FALSE;
        if (!params.FetchValueArrayINT(_T("ThreadDecompose"), intValues))
        {
            appCrucial(_T("ThreadAutoDecompose = 0 but not asign ThreadDecompose, will use auto decompose"));
            bAutoDecompose = TRUE;
        }
        else if (3 != intValues.Num())
        {
            appCrucial(_T("ThreadAutoDecompose = 0 but not ThreadDecompose is invalid, will use auto decompose"));
            bAutoDecompose = TRUE;
        }
        else if (intValues[0] < 1 
              || intValues[1] < 1 
              || intValues[2] < 1
              || deviceConstraints[1] < (UINT)intValues[0]
              || deviceConstraints[2] < (UINT)intValues[1]
              || deviceConstraints[3] < (UINT)intValues[2]
              || deviceConstraints[0] < (UINT)(intValues[0] * intValues[1] * intValues[2])
              || !__Divisible(constIntegers[ECI_Lx], (UINT)intValues[0])
              || !__Divisible(constIntegers[ECI_Ly], (UINT)intValues[1])
              || !__Divisible(constIntegers[ECI_Lz], (UINT)intValues[2])
            )
        {
            appCrucial(_T("ThreadAutoDecompose = 0 but not ThreadDecompose is invalid (should >= 1, should be divisible by lattice length, should < max thread constraints), will use auto decompose"));
            bAutoDecompose = TRUE;
        }
        else
        {
            //use the decompose in param
            constIntegers[ECI_DecompLx] = static_cast<UINT>(intValues[0]);
            constIntegers[ECI_DecompLy] = static_cast<UINT>(intValues[1]);
            constIntegers[ECI_DecompLz] = static_cast<UINT>(intValues[2]);
            constIntegers[ECI_DecompX] = constIntegers[ECI_Lx] / constIntegers[ECI_DecompLx];
            constIntegers[ECI_DecompY] = constIntegers[ECI_Ly] / constIntegers[ECI_DecompLy];
            constIntegers[ECI_DecompZ] = constIntegers[ECI_Lz] / constIntegers[ECI_DecompLz];
        }
    }

    if (bAutoDecompose)
    {
        TArray <UINT> latticeSize;
        latticeSize.AddItem(constIntegers[ECI_Lx]);
        latticeSize.AddItem(constIntegers[ECI_Ly]);
        latticeSize.AddItem(constIntegers[ECI_Lz]);
        TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeSize);

        constIntegers[ECI_DecompX] = decomp[0];
        constIntegers[ECI_DecompY] = decomp[1];
        constIntegers[ECI_DecompZ] = decomp[2];
        constIntegers[ECI_DecompLx] = decomp[3];
        constIntegers[ECI_DecompLy] = decomp[4];
        constIntegers[ECI_DecompLz] = decomp[5];
    }
    appGeneral(_T("\n will run on lattice (%d,%d,%d,%d) with (%d x %d x %d) blocks and (%d x %d x %d) threads per block\n")
        , constIntegers[ECI_Lx]
        , constIntegers[ECI_Ly]
        , constIntegers[ECI_Lz]
        , constIntegers[ECI_Lt]
        , constIntegers[ECI_DecompX]
        , constIntegers[ECI_DecompY]
        , constIntegers[ECI_DecompZ]
        , constIntegers[ECI_DecompLx]
        , constIntegers[ECI_DecompLy]
        , constIntegers[ECI_DecompLz]
    );

#pragma endregion

#pragma region Fill constant table

    __FetchIntWithDefault(_T("RandomSeed"), 1234567);
    constIntegers[ECI_RandomSeed] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("ExponentialPrecision"), 4);
    constIntegers[ECI_ExponentPrecision] = static_cast<UINT>(iVaules);

    const CParameters subparam_updator = params.GetParameter(_T("Updator"));
    if (!subparam_updator.FetchValueINT(_T("IntegratorStep"), iVaules))
    { 
        iVaules = 50; 
    }
    constIntegers[ECI_IntegratorStepCount] = static_cast<UINT>(iVaules);

    CCString sValues;

    __FetchStringWithDefault(_T("RandomType"), _T("ER_Schrage"));
    ERandom eR = __STRING_TO_ENUM(ERandom, sValues);

    if (eR == ER_Schrage)
    {
        constIntegers[ECI_UsingSchrageRandom] = 1;
    }
    else 
    {
        constIntegers[ECI_UsingSchrageRandom] = 0;
    }

#pragma endregion

    m_pCudaHelper = new CCudaHelper();
    memcpy(m_pCudaHelper->m_ConstIntegers, constIntegers, sizeof(UINT) * kContentLength);
    memcpy(m_pCudaHelper->m_ConstFloats, constFloats, sizeof(FLOAT) * kContentLength);
    m_pCudaHelper->CopyConstants();

#pragma region Create Random

    m_pLatticeData = new CLatticeData();
    //appGeneral("are we here?");
    if (eR != ER_Schrage)
    {
        //create another random
    }
    else 
    {
        m_pLatticeData->m_pRandomSchrage = new CRandomSchrage(constIntegers[ECI_RandomSeed]);
        checkCudaErrors(cudaMalloc((void**)&(m_pLatticeData->m_pDeviceRandomSchrage), sizeof(CRandomSchrage)));
        checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceRandomSchrage, m_pLatticeData->m_pRandomSchrage, sizeof(CRandomSchrage), cudaMemcpyHostToDevice));
        appGeneral(_T("Create the Schrage random with seed:%d\n", constIntegers[ECI_RandomSeed]));
    }
    m_pCudaHelper->CopyRandomPointer(m_pLatticeData->m_pDeviceRandom, m_pLatticeData->m_pDeviceRandomSchrage);

#pragma endregion

#pragma region Create Fields

#pragma region Gauge

    CParameters gauge = params.GetParameter(_T("Gauge"));

    CCString sGaugeClassName;
    __FetchStringWithDefaultSub(gauge, _T("FieldName"), _T("CFieldGaugeSU3"));
    sGaugeClassName = sValues;
    __FetchStringWithDefaultSub(gauge, _T("FieldInitialType"), _T("EFIT_Random"));
    EFieldInitialType eGaugeInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    CBase* pGaugeField = appCreate(sGaugeClassName);
    CFieldGauge* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldGauge*>(pGaugeField)) : NULL;
    if (NULL == pGauge)
    {
        appCrucial(_T("Unable to create the gauge field! with name %s!"), sGaugeClassName);
    }

    if (EFIT_ReadFromFile != eGaugeInitial)
    {
        pGauge->m_pOwner = m_pLatticeData;
        pGauge->InitialField(eGaugeInitial);
    }

    m_pLatticeData->m_pGaugeField = pGauge;
    checkCudaErrors(cudaMalloc((void**)&(m_pLatticeData->m_pDeviceGaugeField), pGauge->GetClass()->GetSize()));
    checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceGaugeField, m_pLatticeData->m_pGaugeField, pGauge->GetClass()->GetSize(), cudaMemcpyHostToDevice));
    
    appGeneral(_T("Create the gauge %s with initial: %s\n"), sGaugeClassName.c_str(), sValues.c_str());

#pragma endregion

#pragma endregion

#pragma region Create Index and Boundary

    __FetchStringWithDefault(_T("LatticeBoundary"), _T("CBoundaryConditionTorusSquare"));

    CBoundaryCondition* pBoundary = dynamic_cast<CBoundaryCondition*>(appCreate(sValues));
    if (NULL == pBoundary)
    {
        appCrucial(_T("Cannot create boundary with name: %s"), sValues.c_str());
        exit(EXIT_FAILURE);
    }
    pBoundary->m_pOwner = m_pLatticeData;

    appGeneral(_T("Create the boundary %s\n"), sValues.c_str());

    __FetchStringWithDefault(_T("LatticeIndex"), _T("CIndexSquare"));

    CIndex* pIndex = dynamic_cast<CIndex*>(appCreate(sValues));
    if (NULL == pIndex)
    {
        appCrucial(_T("Cannot create index with name: %s"), sValues.c_str());
        exit(EXIT_FAILURE);
    }
    pIndex->m_pOwner = m_pLatticeData;
    pIndex->m_pBoundaryCondition = pBoundary;
    checkCudaErrors(cudaMalloc((void**)&(pIndex->m_pDeviceBoundaryCondition), pBoundary->GetClass()->GetSize()));
    checkCudaErrors(cudaMemcpy(pIndex->m_pDeviceBoundaryCondition, pIndex->m_pBoundaryCondition, pBoundary->GetClass()->GetSize(), cudaMemcpyHostToDevice));

    m_pLatticeData->m_pIndex = pIndex;

    checkCudaErrors(cudaMalloc((void**)&(m_pLatticeData->m_pDeviceIndex), pIndex->GetClass()->GetSize()));
    checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceIndex, m_pLatticeData->m_pIndex, pIndex->GetClass()->GetSize(), cudaMemcpyHostToDevice));

    appGeneral(_T("Create the index %s\n"), sValues.c_str());


#pragma endregion

#pragma region Create Updator

#pragma endregion

    checkCudaErrors(cudaDeviceSynchronize());
}

void CCLGLibManager::Quit()
{
    appSafeDelete(m_pLatticeData);
    appSafeDelete(m_pCudaHelper);
}

void CLGAPI appInitialCLG(const TCHAR* paramFileName) 
{ 
    CParameters params;
    CYAMLParser::ParseFile(paramFileName, params);
    GCLGManager.InitialWithParameter(params);
}

void CLGAPI appInitialCLG(CParameters& params)
{
    GCLGManager.InitialWithParameter(params);
}

void CLGAPI appQuitCLG() 
{ 
    GCLGManager.Quit(); 
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================