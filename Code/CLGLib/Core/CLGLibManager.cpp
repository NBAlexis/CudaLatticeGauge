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
    Real constFloats[kContentLength];

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
            constIntegers[ECI_ThreadCount] = intValues[0] * intValues[1] * intValues[2];
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
        constIntegers[ECI_ThreadCount] = decomp[3] * decomp[4] * decomp[5];
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

    __FetchIntWithDefault(_T("ActionListLength"), 1);
    constIntegers[ECI_ActionListLength] = static_cast<UINT>(iVaules);

    const CParameters subparam_updator = params.GetParameter(_T("Updator"));
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

    m_pCudaHelper = new CCudaHelper();
    memcpy(m_pCudaHelper->m_ConstIntegers, constIntegers, sizeof(UINT) * kContentLength);
    memcpy(m_pCudaHelper->m_ConstFloats, constFloats, sizeof(Real) * kContentLength);
    m_pCudaHelper->CopyConstants();
    m_pCudaHelper->AllocateTemeraryBuffers(_HC_ThreadCount);

#pragma endregion

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
        appGeneral(_T("Create the Schrage random with seed:%d\n"), constIntegers[ECI_RandomSeed]);
    }
    m_pCudaHelper->CopyRandomPointer(m_pLatticeData->m_pDeviceRandom, m_pLatticeData->m_pDeviceRandomSchrage);

#pragma endregion

#pragma region Create Gamma matrix set

    m_pCudaHelper->CreateGammaMatrix();

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
        appCrucial(_T("Unable to create the gauge field! with name %s!"), sGaugeClassName.c_str());
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

    __FetchStringWithDefault(_T("LatticeBoundary"), _T("EBC_TorusSquare"));

    UINT sizeBufferHost[1];
    sizeBufferHost[0] = 0;
    UINT* sizeBuffer;
    checkCudaErrors(cudaMalloc((void**)&sizeBuffer, sizeof(UINT)));
    checkCudaErrors(cudaMemcpy(sizeBuffer, sizeBufferHost, sizeof(UINT), cudaMemcpyHostToDevice));

    EBoundaryCondition eBC = __STRING_TO_ENUM(EBoundaryCondition, sValues);
    deviceBoundaryCondition ** devicePtrBC;
    checkCudaErrors(cudaMalloc((void**)&devicePtrBC, sizeof(deviceBoundaryCondition *)));

    _cCreateBC((void**)devicePtrBC, sizeBuffer, eBC);
    checkCudaErrors(cudaMemcpy(sizeBufferHost, sizeBuffer, sizeof(UINT), cudaMemcpyDeviceToHost));
    UINT iBoundaryClassSize = sizeBufferHost[0];

    if (0 == iBoundaryClassSize)
    {
        appCrucial(_T("Create Boundary Condition failed! %s"), sValues.c_str());
        exit(EXIT_FAILURE);
    }

    __FetchStringWithDefault(_T("LatticeIndex"), _T("EIndexType_Square"));
    EIndexType eIT = __STRING_TO_ENUM(EIndexType, sValues);
    
    CIndex ** devicePtrIndex;
    checkCudaErrors(cudaMalloc((void**)&devicePtrIndex, sizeof(CIndex *)));
    sizeBufferHost[0] = 0;
    checkCudaErrors(cudaMemcpy(sizeBuffer, sizeBufferHost, sizeof(UINT), cudaMemcpyHostToDevice));
    _cCreateIndex((void**)devicePtrIndex, devicePtrBC, sizeBuffer, eIT);
    checkCudaErrors(cudaMemcpy(sizeBufferHost, sizeBuffer, sizeof(UINT), cudaMemcpyDeviceToHost));
    UINT indexClassSize = sizeBufferHost[0];

    if (0 == indexClassSize)
    {
        appCrucial(_T("Create Index Failed!!!: %s"), sValues.c_str());
        exit(EXIT_FAILURE);
    }

    //Now, we need to copy the content of ptr to lattice
    CIndex* ppIndexHost[1];
    checkCudaErrors(cudaMemcpy(ppIndexHost, devicePtrIndex, sizeof(CIndex**), cudaMemcpyDeviceToHost));
    m_pLatticeData->m_pDeviceIndex = ppIndexHost[0];
    m_pCudaHelper->SetDeviceIndex(devicePtrIndex);

    checkCudaErrors(cudaFree(sizeBuffer));
    checkCudaErrors(cudaFree(devicePtrBC));
    checkCudaErrors(cudaFree(devicePtrIndex));

    appGeneral(_T("Create the index %s\n"), sValues.c_str());

#pragma endregion

#pragma region Craete Actions

    CCString sActionNameList;
    TArray<CAction*> actions;
    for (UINT i = 0; i < constIntegers[ECI_ActionListLength]; ++i)
    {
        CCString sActionParamName;
        sActionParamName.Format(_T("Action%d"), i + 1);
        const CParameters subparam_action = params.GetParameter(sActionParamName);
        CCString sActionName;
        CAction* pAction = NULL;
        if (subparam_action.FetchStringValue(_T("ActionName"), sActionName))
        {
            pAction = dynamic_cast<CAction*>(appCreate(sActionName));
            if (NULL != pAction)
            {
                pAction->Initial(m_pLatticeData, subparam_action);
                actions.AddItem(pAction);

                sActionNameList += (CCString(_T(" ")) + pAction->GetClass()->GetName() + _T(" "));
            }
            else
            {
                //We have already set the constant ECI_ActionListLength
                //So, NULL is not allowed!
                appCrucial(_T("Create Action Failed: %s\n"), sActionName.c_str());
                exit(EXIT_FAILURE);
            }
        }
    }    
    m_pLatticeData->m_pActionList = actions;

    appGeneral(_T("Create the action list, with %d actions: %s\n"), actions.Num(), sActionNameList.c_str());

#pragma endregion

#pragma region Create Updator

    __FetchStringWithDefaultSub(subparam_updator, _T("UpdatorType"), _T("CHMC"));
    CUpdator* updator = dynamic_cast<CUpdator*>(appCreate(sValues));
    CCString sUpdatorInfo = sValues;
    if (NULL != updator && EUT_HMC == updator->GetUpdatorType())
    {
        CHMC* pHMC = dynamic_cast<CHMC*>(updator);
        __FetchStringWithDefaultSub(subparam_updator, _T("IntegratorType"), _T("CIntegratorLeapFrog"));
        CIntegrator * integrator = dynamic_cast<CIntegrator *>(appCreate(sValues));

        if (NULL == pHMC || NULL == integrator)
        {
            appCrucial(_T("HMC need a integrator!, but s = %s"), sValues.c_str());
            exit(EXIT_FAILURE);
        }

        sUpdatorInfo += (" Integrator:" + sValues);
        integrator->Initial(pHMC, m_pLatticeData, subparam_updator);
        pHMC->Initial(m_pLatticeData, subparam_updator);
        pHMC->m_pIntegrator = integrator;
        m_pLatticeData->m_pUpdator = pHMC;
    }
    else
    {
        appCrucial(_T("Failed to create Updator! s = %s"), sValues);
        exit(EXIT_FAILURE);
    }
    
    appGeneral(_T("Create Updator %s\n"), sUpdatorInfo);

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