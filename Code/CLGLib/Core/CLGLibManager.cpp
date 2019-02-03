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

void CCLGLibManager::SetupLog(CParameters &params)
{
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
}

void CCLGLibManager::InitialLatticeAndConstant(CParameters& params)
{
    INT iVaules = 0;
    CCString sValues;

#pragma region Lattice Size and Threads

    __FetchIntWithDefault(_T("Dim"), 4);
    assert(iVaules > 1 && iVaules < 5);
    m_InitialCache.constIntegers[ECI_Dim] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Dir"), 4);
    assert(iVaules > 1);
    m_InitialCache.constIntegers[ECI_Dir] = static_cast<UINT>(iVaules);

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

    m_InitialCache.constIntegers[ECI_Lx] = static_cast<UINT>(intValues[0]);
    m_InitialCache.constIntegers[ECI_Ly] = static_cast<UINT>(intValues[1]);
    m_InitialCache.constIntegers[ECI_Lz] = static_cast<UINT>(intValues[2]);
    m_InitialCache.constIntegers[ECI_Lt] = static_cast<UINT>(intValues[3]);
    m_InitialCache.constIntegers[ECI_Volumn] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultX] = static_cast<UINT>(intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultY] = static_cast<UINT>(intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultZ] = static_cast<UINT>(intValues[3]);

    m_InitialCache.constIntegers[ECI_PlaqutteCount] = m_InitialCache.constIntegers[ECI_Volumn] * m_InitialCache.constIntegers[ECI_Dir] * (m_InitialCache.constIntegers[ECI_Dir] - 1) / 2;
    m_InitialCache.constIntegers[ECI_LinkCount] = m_InitialCache.constIntegers[ECI_Volumn] * m_InitialCache.constIntegers[ECI_Dir];

    m_InitialCache.constFloats[ECF_InverseSqrtLink16] = 1 / _sqrt(16 * m_InitialCache.constIntegers[ECI_LinkCount]);
    appGeneral(_T("ECF_InverseSqrtLink16:%f"), m_InitialCache.constFloats[ECF_InverseSqrtLink16]);
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
            || !__Divisible(m_InitialCache.constIntegers[ECI_Lx], (UINT)intValues[0])
            || !__Divisible(m_InitialCache.constIntegers[ECI_Ly], (UINT)intValues[1])
            || !__Divisible(m_InitialCache.constIntegers[ECI_Lz], (UINT)intValues[2])
            )
        {
            appCrucial(_T("ThreadAutoDecompose = 0 but not ThreadDecompose is invalid (should >= 1, should be divisible by lattice length, should < max thread constraints), will use auto decompose"));
            bAutoDecompose = TRUE;
        }
        else
        {
            //use the decompose in param
            m_InitialCache.constIntegers[ECI_DecompLx] = static_cast<UINT>(intValues[0]);
            m_InitialCache.constIntegers[ECI_DecompLy] = static_cast<UINT>(intValues[1]);
            m_InitialCache.constIntegers[ECI_DecompLz] = static_cast<UINT>(intValues[2]);
            m_InitialCache.constIntegers[ECI_DecompX] = m_InitialCache.constIntegers[ECI_Lx] / m_InitialCache.constIntegers[ECI_DecompLx];
            m_InitialCache.constIntegers[ECI_DecompY] = m_InitialCache.constIntegers[ECI_Ly] / m_InitialCache.constIntegers[ECI_DecompLy];
            m_InitialCache.constIntegers[ECI_DecompZ] = m_InitialCache.constIntegers[ECI_Lz] / m_InitialCache.constIntegers[ECI_DecompLz];
            m_InitialCache.constIntegers[ECI_ThreadCount] = m_InitialCache.constIntegers[ECI_Lx] * m_InitialCache.constIntegers[ECI_Ly] * m_InitialCache.constIntegers[ECI_Lz];
            m_InitialCache.constIntegers[ECI_ThreadCountPerBlock] = intValues[0] * intValues[1] * intValues[2];
        }
    }

    if (bAutoDecompose)
    {
        TArray <UINT> latticeSize;
        latticeSize.AddItem(m_InitialCache.constIntegers[ECI_Lx]);
        latticeSize.AddItem(m_InitialCache.constIntegers[ECI_Ly]);
        latticeSize.AddItem(m_InitialCache.constIntegers[ECI_Lz]);
        TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeSize);

        m_InitialCache.constIntegers[ECI_DecompX] = decomp[0];
        m_InitialCache.constIntegers[ECI_DecompY] = decomp[1];
        m_InitialCache.constIntegers[ECI_DecompZ] = decomp[2];
        m_InitialCache.constIntegers[ECI_DecompLx] = decomp[3];
        m_InitialCache.constIntegers[ECI_DecompLy] = decomp[4];
        m_InitialCache.constIntegers[ECI_DecompLz] = decomp[5];
        m_InitialCache.constIntegers[ECI_ThreadCountPerBlock] = decomp[3] * decomp[4] * decomp[5];
        m_InitialCache.constIntegers[ECI_ThreadCount] = decomp[0] * decomp[1] * decomp[2] * decomp[3] * decomp[4] * decomp[5];
    }
    appGeneral(_T("\n will run on lattice (%d,%d,%d,%d) with (%d x %d x %d) blocks and (%d x %d x %d) threads per block\n")
        , m_InitialCache.constIntegers[ECI_Lx]
        , m_InitialCache.constIntegers[ECI_Ly]
        , m_InitialCache.constIntegers[ECI_Lz]
        , m_InitialCache.constIntegers[ECI_Lt]
        , m_InitialCache.constIntegers[ECI_DecompX]
        , m_InitialCache.constIntegers[ECI_DecompY]
        , m_InitialCache.constIntegers[ECI_DecompZ]
        , m_InitialCache.constIntegers[ECI_DecompLx]
        , m_InitialCache.constIntegers[ECI_DecompLy]
        , m_InitialCache.constIntegers[ECI_DecompLz]
    );

#pragma endregion

#pragma region Fill constant table

    __FetchIntWithDefault(_T("RandomSeed"), 1234567);
    m_InitialCache.constIntegers[ECI_RandomSeed] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("ExponentialPrecision"), 8);
    m_InitialCache.constIntegers[ECI_ExponentPrecision] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("ActionListLength"), 0);
    m_InitialCache.constIntegers[ECI_ActionListLength] = static_cast<UINT>(iVaules);

    __FetchStringWithDefault(_T("RandomType"), _T("ER_Schrage"));
    m_InitialCache.eR = __STRING_TO_ENUM(ERandom, sValues);

    __FetchIntWithDefault(_T("MeasureListLength"), 0);
    m_InitialCache.constIntegers[ECI_MeasureListLength] = static_cast<UINT>(iVaules);

    m_InitialCache.constIntegers[ECI_SUN] = 1;
    if (params.Exist(_T("Gauge")))
    {
        CParameters gauge = params.GetParameter(_T("Gauge"));
        if (gauge.Exist(_T("FieldName")))
        {
            CCString sFieldName;
            gauge.FetchStringValue(_T("FieldName"), sFieldName);
            if (sFieldName == _T("CFieldGaugeSU3"))
            {
                m_InitialCache.constIntegers[ECI_SUN] = 3;
            }
        }
    }

    memcpy(m_pCudaHelper->m_ConstIntegers, m_InitialCache.constIntegers, sizeof(UINT) * kContentLength);
    memcpy(m_pCudaHelper->m_ConstFloats, m_InitialCache.constFloats, sizeof(Real) * kContentLength);
    m_pCudaHelper->CopyConstants();
    m_pCudaHelper->AllocateTemeraryBuffers(_HC_ThreadCount);

#pragma endregion

    m_pCudaHelper->CreateGammaMatrix();
}

void CCLGLibManager::InitialRandom(CParameters &)
{
    //INT iVaules = 0;
    //CCString sValues;

    m_pLatticeData->m_pRandom = new CRandom(m_InitialCache.constIntegers[ECI_RandomSeed], m_InitialCache.eR);
    checkCudaErrors(cudaMalloc((void**)&(m_pLatticeData->m_pDeviceRandom), sizeof(CRandom)));
    checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceRandom, m_pLatticeData->m_pRandom, sizeof(CRandom), cudaMemcpyHostToDevice));
    appGeneral(_T("Create the %s random with seed:%d\n"), __ENUM_TO_STRING(ERandom, m_InitialCache.eR).c_str(), m_InitialCache.constIntegers[ECI_RandomSeed]);

    m_pCudaHelper->CopyRandomPointer(m_pLatticeData->m_pDeviceRandom);
}

void CCLGLibManager::CreateGaugeField(class CParameters& params)
{
    INT iVaules = 0;
    CCString sValues;

    CCString sGaugeClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldGaugeSU3"));
    sGaugeClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    EFieldInitialType eGaugeInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    CBase* pGaugeField = appCreate(sGaugeClassName);
    CFieldGauge* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldGauge*>(pGaugeField)) : NULL;
    if (NULL == pGauge)
    {
        appCrucial(_T("Unable to create the gauge field! with name %s!"), sGaugeClassName.c_str());
    }
    pGauge->m_byFieldId = 1;
    if (EFIT_ReadFromFile != eGaugeInitial)
    {
        pGauge->m_pOwner = m_pLatticeData;
        pGauge->InitialField(eGaugeInitial);
    }
    else
    {
        CCString sFileType, sFileName;
        if (!params.FetchStringValue(_T("GaugeFileType"), sFileType)
         || !params.FetchStringValue(_T("GaugeFileName"), sFileName))
        {
            appCrucial(_T("Gauge initial type is EFIT_ReadFromFile, but cannot find GaugeFileType or GaugeFileName!\n"));
            exit(EXIT_FAILURE);
        }
        EFieldFileType eFileType = __STRING_TO_ENUM(EFieldFileType, sFileType);
        pGauge->m_pOwner = m_pLatticeData;
        pGauge->InitialFieldWithFile(sFileName, eFileType);
    }

    m_pLatticeData->m_pGaugeField = pGauge;
    m_pLatticeData->m_pFieldMap.SetAt(1, pGauge);
    if (NULL != m_pLatticeData->m_pDeviceIndex)
    {
        pGauge->CachePlaqutteIndexes();
    }
    //checkCudaErrors(cudaMalloc((void**)&(m_pLatticeData->m_pDeviceGaugeField), pGauge->GetClass()->GetSize()));
    //checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceGaugeField, m_pLatticeData->m_pGaugeField, pGauge->GetClass()->GetSize(), cudaMemcpyHostToDevice));

    appGeneral(_T("Create the gauge %s with initial: %s\n"), sGaugeClassName.c_str(), sValues.c_str());
}

void CCLGLibManager::CreateOtherFields(class CParameters& params)
{

}

void CCLGLibManager::CreateIndexAndBoundary(class CParameters& params)
{
    INT iVaules = 0;
    CCString sValues;

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
}

void CCLGLibManager::CreateActionList(class CParameters& params)
{
    CCString sActionNameList;
    TArray<CAction*> actions;
    for (UINT i = 0; i < m_InitialCache.constIntegers[ECI_ActionListLength]; ++i)
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
                pAction->Initial(m_pLatticeData, subparam_action, i + 1);
                actions.AddItem(pAction);
                m_pLatticeData->m_pActionMap.SetAt(i + 1, pAction);
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
}

void CCLGLibManager::CreateUpdator(class CParameters& params)
{
    CCString sValues;

    __FetchStringWithDefault( _T("UpdatorType"), _T("CHMC"));
    CUpdator* updator = dynamic_cast<CUpdator*>(appCreate(sValues));
    CCString sUpdatorInfo = sValues;
    if (NULL != updator && EUT_HMC == updator->GetUpdatorType())
    {
        CHMC* pHMC = dynamic_cast<CHMC*>(updator);
        __FetchStringWithDefault(_T("IntegratorType"), _T("CIntegratorLeapFrog"));
        CIntegrator * integrator = dynamic_cast<CIntegrator *>(appCreate(sValues));

        if (NULL == pHMC || NULL == integrator)
        {
            appCrucial(_T("HMC need a integrator!, but s = %s"), sValues.c_str());
            exit(EXIT_FAILURE);
        }

        sUpdatorInfo += (" Integrator:" + sValues);
        integrator->Initial(pHMC, m_pLatticeData, params);
        pHMC->Initial(m_pLatticeData, params);
        pHMC->m_pIntegrator = integrator;
        m_pLatticeData->m_pUpdator = pHMC;
    }
    else
    {
        appCrucial(_T("Failed to create Updator! s = %s"), sValues.c_str());
        exit(EXIT_FAILURE);
    }

    appGeneral(_T("Create Updator %s\n"), sUpdatorInfo.c_str());
}

void CCLGLibManager::CreateMeasurement(class CParameters& params)
{
    CCString sMeasureNameList;
    CMeasurementManager* pMeasurements = new CMeasurementManager(m_pLatticeData);
    for (UINT i = 0; i < m_InitialCache.constIntegers[ECI_MeasureListLength]; ++i)
    {
        CCString sMeasureParamName;
        sMeasureParamName.Format(_T("Measure%d"), i + 1);

        const CParameters subparam_measure = params.GetParameter(sMeasureParamName);
        CCString sMeasureName;
        CMeasure* pMeasure = NULL;
        if (subparam_measure.FetchStringValue(_T("MeasureName"), sMeasureName))
        {
            pMeasure = dynamic_cast<CMeasure*>(appCreate(sMeasureName));
            if (NULL != pMeasure)
            {
                pMeasure->Initial(pMeasurements, m_pLatticeData, subparam_measure, i + 1);
                pMeasurements->m_lstAllMeasures.AddItem(pMeasure);
                pMeasurements->m_mapMeasures.SetAt(i + 1, pMeasure);
                sMeasureNameList += (CCString(_T(" ")) + pMeasure->GetClass()->GetName() + _T(" "));
            }
            else
            {
                //We have already set the constant ECI_ActionListLength
                //So, NULL is not allowed!
                appCrucial(_T("Create Measure Failed: %s\n"), sMeasureName.c_str());
                exit(EXIT_FAILURE);
            }
        }
    }

    m_pLatticeData->m_pMeasurements = pMeasurements;

    appGeneral(_T("Create the measure list, with %d measures: %s\n"), pMeasurements->m_lstAllMeasures.Num(), sMeasureNameList.c_str());
}

void CCLGLibManager::CreateSolver(class CParameters& params)
{

}

UBOOL CCLGLibManager::InitialWithParameter(CParameters &params)
{
    checkCudaErrors(cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1 << 30));

    m_pCudaHelper = new CCudaHelper();
    m_pLatticeData = new CLatticeData();
    m_pFileSystem = new CFileSystem();

    InitialLatticeAndConstant(params);
    InitialRandom(params);

    if (params.Exist(_T("LatticeIndex")))
    {
        CreateIndexAndBoundary(params);
    }

    if (params.Exist(_T("Gauge")))
    {
        CParameters gauge = params.GetParameter(_T("Gauge"));
        CreateGaugeField(gauge);
    }

    if (m_InitialCache.constIntegers[ECI_ActionListLength] > 0)
    {
        CreateActionList(params);
    }

    if (params.Exist(_T("Updator")))
    {
        CParameters updator = params.GetParameter(_T("Updator"));
        CreateUpdator(updator);
    }

    if (m_InitialCache.constIntegers[ECI_MeasureListLength] > 0)
    {
        CreateMeasurement(params);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    cudaError_t cudaEr = cudaGetLastError();
    if (cudaEr != cudaSuccess)
    {
        return FALSE;
    }
    return TRUE;
}

void CCLGLibManager::Quit()
{
    appSafeDelete(m_pLatticeData);
    appSafeDelete(m_pCudaHelper);
    appSafeDelete(m_pFileSystem);

    INT devCount;
    cudaGetDeviceCount(&devCount);
    for (INT i = 0; i < devCount; ++i)
    {
        cudaSetDevice(i);
        cudaDeviceReset();
    }
}

UBOOL CLGAPI appInitialCLG(const TCHAR* paramFileName)
{ 
    CParameters params;
    CYAMLParser::ParseFile(paramFileName, params);
    return GCLGManager.InitialWithParameter(params);
}

UBOOL CLGAPI appInitialCLG(CParameters& params)
{
    return GCLGManager.InitialWithParameter(params);
}

void CLGAPI appQuitCLG() 
{ 
    GCLGManager.Quit(); 
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================