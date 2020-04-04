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

__BEGIN_NAMESPACE

CLGAPI CCLGLibManager GCLGManager;

void CCLGLibManager::SetupLog(CParameters &params)
{
    //Setup outputs
    CCString verboselevel;
    EVerboseLevel eVerbLevel = CRUCIAL;
    CCString sVerbFile = _T("stdout");
    const UBOOL fetchVerbLevel = params.FetchStringValue(_T("VerboseLevel"), verboselevel);
    const UBOOL fetchVerbFile = params.FetchStringValue(_T("VerboseOutput"), sVerbFile);
    if (fetchVerbLevel || fetchVerbFile) //do NOT put fetch string in if, it will enter if when the first is TRUE
    {
        eVerbLevel = __STRING_TO_ENUM(EVerboseLevel, verboselevel);
        appSetTracer(eVerbLevel, sVerbFile);
    }

    //check whether to log parameter file
    INT iTag = 0;
    //appGeneral(_T("============================== Parameter =============================\n\n"));
    __CheckTag(_T("ShowParameterContent"), params.Dump());
    //appGeneral(_T("============================== GPU =============================\n\n"));
    __CheckTag(_T("ShowDeviceInformation"), CCudaHelper::DeviceQuery());

    appGeneral(_T("============================== Log Start =============================\n\n"));
}

#pragma region Creates

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

    __FetchIntWithDefault(_T("MaxThreadPerBlock"), 0);
    if (iVaules > 0)
    {
        CCommonData::m_uiMaxThreadPerBlock = static_cast<UINT>(iVaules);
    }

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
    CCommonData::m_sCenter.x = static_cast<SBYTE>(m_InitialCache.constIntegers[ECI_Lx] / 2);
    CCommonData::m_sCenter.y = static_cast<SBYTE>(m_InitialCache.constIntegers[ECI_Ly] / 2);
    CCommonData::m_sCenter.z = static_cast<SBYTE>(m_InitialCache.constIntegers[ECI_Lz] / 2);
    CCommonData::m_sCenter.w = static_cast<SBYTE>(m_InitialCache.constIntegers[ECI_Lt] / 2);
    m_InitialCache.constIntegers[ECI_Volume] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_Volume_xyz] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2]);
    m_InitialCache.constIntegers[ECI_MultX] = static_cast<UINT>(intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultY] = static_cast<UINT>(intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultZ] = static_cast<UINT>(intValues[3]);
    m_InitialCache.constIntegers[ECI_GridDimZT] = m_InitialCache.constIntegers[ECI_Lz] * m_InitialCache.constIntegers[ECI_Lt];
    TArray<UINT> latticeDim;
    latticeDim.AddItem(intValues[0] * intValues[1]); //xy
    latticeDim.AddItem(intValues[2]); //z
    latticeDim.AddItem(intValues[3]); //t

    m_InitialCache.constIntegers[ECI_PlaqutteCount] = m_InitialCache.constIntegers[ECI_Volume] * m_InitialCache.constIntegers[ECI_Dir] * (m_InitialCache.constIntegers[ECI_Dir] - 1) / 2;
    m_InitialCache.constIntegers[ECI_LinkCount] = m_InitialCache.constIntegers[ECI_Volume] * m_InitialCache.constIntegers[ECI_Dir];

    m_InitialCache.constFloats[ECF_InverseSqrtLink16] = F(1.0) / _sqrt(F(16.0) * m_InitialCache.constIntegers[ECI_LinkCount]);
    
    UBOOL bAutoDecompose = TRUE;
    __FetchIntWithDefault(_T("ThreadAutoDecompose"), 1);

    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();


    m_InitialCache.constIntegers[ECI_ThreadConstaint] = deviceConstraints[0];
    m_InitialCache.constIntegers[ECI_ThreadConstaintX] = deviceConstraints[1];
    m_InitialCache.constIntegers[ECI_ThreadConstaintY] = deviceConstraints[2];
    m_InitialCache.constIntegers[ECI_ThreadConstaintZ] = deviceConstraints[3];

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
            || !__Divisible(latticeDim[0], (UINT)intValues[0])
            || !__Divisible(latticeDim[1], (UINT)intValues[1])
            || !__Divisible(latticeDim[2], (UINT)intValues[2])
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
            m_InitialCache.constIntegers[ECI_DecompX] = latticeDim[0] / m_InitialCache.constIntegers[ECI_DecompLx];
            m_InitialCache.constIntegers[ECI_DecompY] = latticeDim[1] / m_InitialCache.constIntegers[ECI_DecompLy];
            m_InitialCache.constIntegers[ECI_DecompZ] = latticeDim[2] / m_InitialCache.constIntegers[ECI_DecompLz];
            m_InitialCache.constIntegers[ECI_ThreadCountPerBlock] = intValues[0] * intValues[1] * intValues[2];
        }
    }

    if (bAutoDecompose)
    {
        TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);

        m_InitialCache.constIntegers[ECI_DecompX] = decomp[0];
        m_InitialCache.constIntegers[ECI_DecompY] = decomp[1];
        m_InitialCache.constIntegers[ECI_DecompZ] = decomp[2];
        m_InitialCache.constIntegers[ECI_DecompLx] = decomp[3];
        m_InitialCache.constIntegers[ECI_DecompLy] = decomp[4];
        m_InitialCache.constIntegers[ECI_DecompLz] = decomp[5];
        m_InitialCache.constIntegers[ECI_ThreadCountPerBlock] = decomp[3] * decomp[4] * decomp[5];
    }
    appGeneral(_T("\n will run on lattice ( (%d,%d),%d,%d) with (xy %d x z %d x t %d) blocks and (xy %d x z %d x t %d) threads per block\n")
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
    CCString sRandomSeedType;
    if (params.FetchStringValue(_T("RandomSeedType"), sRandomSeedType))
    {
        const ERandomSeedType eRST = __STRING_TO_ENUM(ERandomSeedType, sRandomSeedType);
        if (ERST_Timestamp == eRST)
        {
            m_InitialCache.constIntegers[ECI_RandomSeed] = appGetTimeStamp();
        }
    }

    __FetchIntWithDefault(_T("ExponentialPrecision"), 0);
    m_InitialCache.constIntegers[ECI_ExponentPrecision] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("CacheStaple"), 1);
    CCommonData::m_bStoreStaple = (0 != iVaules);

    __FetchIntWithDefault(_T("StochasticGaussian"), 0);
    CCommonData::m_bStochasticGaussian = (0 != iVaules);
    
    __FetchIntWithDefault(_T("CacheSolution"), 1);
    CCommonData::m_bStoreLastSolution = (0 != iVaules);

    __FetchIntWithDefault(_T("ActionListLength"), 0);
    m_InitialCache.constIntegers[ECI_ActionListLength] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("FermionFieldCount"), 0);
    m_InitialCache.constIntegers[ECI_FermionFieldLength] = static_cast<UINT>(iVaules);

    __FetchStringWithDefault(_T("RandomType"), _T("ER_Schrage"));
    m_InitialCache.eR = __STRING_TO_ENUM(ERandom, sValues);

    __FetchIntWithDefault(_T("MeasureListLength"), 0);
    m_InitialCache.constIntegers[ECI_MeasureListLength] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("UseLogADefinition"), 0);
    m_InitialCache.constIntegers[ECI_UseLogADefinition] = static_cast<UINT>(iVaules);

    const UINT iThreadConstraint = m_InitialCache.constIntegers[ECI_ThreadConstaint];
    __FetchIntWithDefault(_T("SummationDecompose"), iThreadConstraint);
    m_InitialCache.constIntegers[ECI_SummationDecompose] = static_cast<UINT>(iVaules);
    appDetailed(_T("Summation decompose: %d\n"), m_InitialCache.constIntegers[ECI_SummationDecompose]);

    m_InitialCache.constIntegers[ECI_SUN] = 1;
    if (params.Exist(_T("Gauge")))
    {
        const CParameters gauge = params.GetParameter(_T("Gauge"));
        if (gauge.Exist(_T("FieldName")))
        {
            CCString sFieldName;
            gauge.FetchStringValue(_T("FieldName"), sFieldName);
            if (sFieldName == _T("CFieldGaugeSU3")
             || sFieldName == _T("CFieldGaugeSU3D"))
            {
                m_InitialCache.constIntegers[ECI_SUN] = 3;
            }
        }
    }

    memcpy(m_pCudaHelper->m_ConstIntegers, m_InitialCache.constIntegers, sizeof(UINT) * kContentLength);
    memcpy(m_pCudaHelper->m_ConstFloats, m_InitialCache.constFloats, sizeof(Real) * kContentLength);
    m_pCudaHelper->CopyConstants();
    m_pCudaHelper->AllocateTemeraryBuffers(_HC_Volume);

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
    m_pLatticeData->m_uiRandomType = static_cast<UINT>(m_InitialCache.eR);
    m_pLatticeData->m_uiRandomSeed = m_InitialCache.constIntegers[ECI_RandomSeed];
}

void CCLGLibManager::CreateGaugeField(class CParameters& params) const
{
    //INT iVaules = 0;
    CCString sValues;

    CCString sGaugeClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldGaugeSU3"));
    sGaugeClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    const EFieldInitialType eGaugeInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    CBase* pGaugeField = appCreate(sGaugeClassName);
    CFieldGauge* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldGauge*>(pGaugeField)) : NULL;
    if (NULL == pGauge)
    {
        appCrucial(_T("Unable to create the gauge field! with name %s!"), sGaugeClassName.c_str());
        return;
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
            _FAIL_EXIT;
        }
        const EFieldFileType eFileType = __STRING_TO_ENUM(EFieldFileType, sFileType);
        pGauge->m_pOwner = m_pLatticeData;
        pGauge->InitialFieldWithFile(sFileName, eFileType);
    }
    pGauge->InitialOtherParameters(params);

    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SBYTE>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SBYTE>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SBYTE>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SBYTE>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(1, bc);
    }

    m_pLatticeData->m_pGaugeField = pGauge;
    m_pLatticeData->m_pFieldMap.SetAt(1, pGauge);

    appGeneral(_T("Create the gauge %s with initial: %s\n"), sGaugeClassName.c_str(), sValues.c_str());
}

void CCLGLibManager::CreateGaugeBoundaryField(class CParameters& params) const
{
    CCString sValues;

    CCString sGaugeClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldBoundaryGaugeSU3"));
    sGaugeClassName = sValues;

    CBase* pGaugeField = appCreate(sGaugeClassName);
    CFieldBoundary* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldBoundary*>(pGaugeField)) : NULL;

    if (NULL == pGauge)
    {
        appCrucial(_T("Unable to create the gauge field! with name %s!"), sGaugeClassName.c_str());
    }
    pGauge->InitialField(params);

    m_pLatticeData->m_pBoundaryFieldMap.SetAt(1, pGauge);
    appGeneral(_T("Create the boundary gauge %s with initial: %s\n"), sGaugeClassName.c_str(), sValues.c_str());
}

void CCLGLibManager::CreateFermionFields(class CParameters& params) const
{
    INT iVaules = 0;
    CCString sValues;

    CCString sFermionClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldFermionWilsonSquareSU3"));
    sFermionClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_RandomGaussian"));
    const EFieldInitialType eFieldInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    CBase* pFermionField = appCreate(sFermionClassName);
    CFieldFermion* pFermion = (NULL != pFermionField) ? (dynamic_cast<CFieldFermion*>(pFermionField)) : NULL;
    if (NULL == pFermion)
    {
        appCrucial(_T("Unable to create the fermion field! with name %s!"), sFermionClassName.c_str());
        return;
    }

    __FetchIntWithDefault(_T("FieldId"), -1);
    const BYTE byFieldId = static_cast<BYTE>(iVaules);
    assert(byFieldId < kMaxFieldCount && byFieldId > 1);
    if (byFieldId >= kMaxFieldCount || byFieldId <= 1)
    {
        appCrucial(_T("The field Id must > 1 and < %d\n"), kMaxFieldCount);
        _FAIL_EXIT;
    }
    if (m_pLatticeData->m_pFieldMap.Exist(byFieldId))
    {
        appCrucial(_T("Unable to create the fermion field! with wrong field ID %s %d!"), sFermionClassName.c_str(), byFieldId);
        return;
    }

    pFermion->m_byFieldId = byFieldId;
    pFermion->m_pOwner = m_pLatticeData;
    pFermion->InitialField(eFieldInitial);
    pFermion->InitialOtherParameters(params);
    m_pLatticeData->m_pFieldMap.SetAt(byFieldId, pFermion);
    m_pLatticeData->m_pOtherFields.AddItem(pFermion);
    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SBYTE>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SBYTE>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SBYTE>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SBYTE>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
    }

    __FetchIntWithDefault(_T("PoolNumber"), 0);
    if (iVaules > 0)
    {
        m_pLatticeData->CreateFieldPool(byFieldId, iVaules);
    }

    appGeneral(_T("Create the fermion field %s with id %d and initial: %s\n"), sFermionClassName.c_str(), byFieldId, sValues.c_str());
}

void CCLGLibManager::CreateFermionBoundaryField(class CParameters& params) const
{
    CCString sValues;
    INT iVaules;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldBoundaryWilsonSquareSU3"));
    const CCString sFieldClassName = sValues;

    CBase* pBCField = appCreate(sFieldClassName);
    CFieldBoundary* pBC = (NULL != pBCField) ? (dynamic_cast<CFieldBoundary*>(pBCField)) : NULL;

    if (NULL == pBC)
    {
        appCrucial(_T("Unable to create the boundary fermion field! with name %s!"), sFieldClassName.c_str());
    }
    pBC->InitialField(params);

    __FetchIntWithDefault(_T("FieldId"), -1);
    const BYTE byFieldId = static_cast<BYTE>(iVaules);
    m_pLatticeData->m_pBoundaryFieldMap.SetAt(byFieldId, pBC);
    appGeneral(_T("Create the boundary fermion field %s with initial: %s\n"), sFieldClassName.c_str(), sValues.c_str());
}

void CCLGLibManager::CreateIndexAndBoundary(class CParameters& params) const
{
    //INT iVaules = 0;
    CCString sValues;

    __FetchStringWithDefault(_T("LatticeBoundary"), _T("CBoundaryConditionTorusSquare"));

    CBoundaryCondition * pBc = dynamic_cast<CBoundaryCondition *>(appCreate(sValues));

    if (NULL == pBc)
    {
        appCrucial(_T("Create Boundary Condition failed! %s"), sValues.c_str());
        _FAIL_EXIT;
    }

    __FetchStringWithDefault(_T("LatticeIndex"), _T("CIndexSquare"));

    CIndex * pIndex = dynamic_cast<CIndex *>(appCreate(sValues));
    if (NULL == pIndex)
    {
        appCrucial(_T("Create Index failed! %s"), sValues.c_str());
        _FAIL_EXIT;
    }
    pIndex->SetBoundaryCondition(pBc);
    m_pLatticeData->m_pIndex = pIndex;

    appGeneral(_T("Create the index %s\n"), sValues.c_str());

    //Index Cache
    m_pLatticeData->m_pIndexCache = new CIndexData();
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
                BYTE byId = static_cast<BYTE>(i + 1);
                pAction->Initial(m_pLatticeData, subparam_action, byId);
                actions.AddItem(pAction);
                m_pLatticeData->m_pActionMap.SetAt(byId, pAction);
                sActionNameList += (CCString(_T(" ")) + pAction->GetClass()->GetName() + _T(" "));
            }
            else
            {
                //We have already set the constant ECI_ActionListLength
                //So, NULL is not allowed!
                appCrucial(_T("Create Action Failed: %s\n"), sActionName.c_str());
                _FAIL_EXIT;
            }
        }
    }
    m_pLatticeData->m_pActionList = actions;

    appGeneral(_T("Create the action list, with %d actions: %s\n"), actions.Num(), sActionNameList.c_str());
}

void CCLGLibManager::CreateUpdator(class CParameters& params) const
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
            appCrucial(_T("HMC need a integrator!, but s = %s\n"), sValues.c_str());
            _FAIL_EXIT;
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
        _FAIL_EXIT;
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
                BYTE byId = static_cast<BYTE>(i + 1);
                pMeasure->Initial(pMeasurements, m_pLatticeData, subparam_measure, byId);
                pMeasurements->m_lstAllMeasures.AddItem(pMeasure);
                pMeasurements->m_mapMeasures.SetAt(byId, pMeasure);
                sMeasureNameList += (CCString(_T(" ")) + pMeasure->GetClass()->GetName() + _T(" "));
            }
            else
            {
                //We have already set the constant ECI_ActionListLength
                //So, NULL is not allowed!
                appCrucial(_T("Create Measure Failed: %s\n"), sMeasureName.c_str());
                _FAIL_EXIT;
            }
        }
    }

    m_pLatticeData->m_pMeasurements = pMeasurements;

    appGeneral(_T("Create the measure list, with %d measures: %s\n"), pMeasurements->m_lstAllMeasures.Num(), sMeasureNameList.c_str());
}

void CCLGLibManager::CreateSolver(class CParameters& params) const
{
    CCString sSolverName = _T("CSLASolverBiCGStab");
    params.FetchStringValue(_T("SolverName"), sSolverName);
    INT byFieldId = 2;
    params.FetchValueINT(_T("SolverForFieldId"), byFieldId);
    CField * pField = m_pLatticeData->GetFieldById(static_cast<BYTE>(byFieldId));
    if (NULL == pField)
    {
        appCrucial(_T("Solver must be created for a specified field!\n"));
    }
    m_pLatticeData->CreateFermionSolver(sSolverName, params, pField, static_cast<BYTE>(byFieldId));
}

void CCLGLibManager::CreateGaugeSmearing(class CParameters& params) const
{
    CCString sSmearingName = _T("CGaugeSmearingAPEStout");
    params.FetchStringValue(_T("SmearingName"), sSmearingName);
    m_pLatticeData->m_pGaugeSmearing = dynamic_cast<CGaugeSmearing*>(appCreate(sSmearingName));
    if (NULL != m_pLatticeData->m_pGaugeSmearing)
    {
        m_pLatticeData->m_pGaugeSmearing->Initial(m_pLatticeData, params);
    }
}

void CCLGLibManager::CreateGaugeFixing(class CParameters& params) const
{
    CCString sSmearingName = _T("CGaugeFixingLandauCornell");
    params.FetchStringValue(_T("Name"), sSmearingName);
    m_pLatticeData->m_pGaugeFixing = dynamic_cast<CGaugeFixing*>(appCreate(sSmearingName));
    if (NULL != m_pLatticeData->m_pGaugeFixing)
    {
        m_pLatticeData->m_pGaugeFixing->Initial(m_pLatticeData, params);
    }
}

#pragma endregion

#pragma region Caches

void CCLGLibManager::InitialIndexBuffer() const
{
    if (NULL == m_pLatticeData->m_pIndexCache)
    {
        appGeneral(_T("No Index Cache"));
        return;
    }

    m_pLatticeData->m_pIndex->BakeAllIndexBuffer(m_pLatticeData->m_pIndexCache);
    if (NULL != m_pLatticeData->m_pGaugeField)
    {
        UBOOL bHasStaggeredFermion = FALSE;
        assert(1 == m_pLatticeData->m_pGaugeField->m_byFieldId);

        m_pLatticeData->m_pIndex->BakePlaquttes(m_pLatticeData->m_pIndexCache, 1);

        for (BYTE i = 2; i < kMaxFieldCount; ++i)
        {
            if (NULL != m_pLatticeData->GetFieldById(i))
            {
                m_pLatticeData->m_pIndex->BakeMoveIndex(m_pLatticeData->m_pIndexCache, i);
                if (EFT_FermionStaggeredSU3 == m_pLatticeData->GetFieldById(i)->GetFieldType())
                {
                    bHasStaggeredFermion = TRUE;
                }
            }
        }
        if (bHasStaggeredFermion)
        {
            m_pLatticeData->m_pIndex->BakeEtaMuTable(m_pLatticeData->m_pIndexCache);
        }
    }
    m_pLatticeData->m_pIndex->CalculateSiteCount(m_pLatticeData->m_pIndexCache);

    m_pCudaHelper->SetDeviceIndex(m_pLatticeData->m_pIndexCache);
}

#pragma endregion

UBOOL CCLGLibManager::InitialWithParameter(CParameters &params)
{
    m_pCudaHelper = new CCudaHelper();
    m_pLatticeData = new CLatticeData();
    m_pFileSystem = new CFileSystem();
    m_pBuffer = new CCudaBuffer();

    //Allocate Buffer
    Real fBufferSize = F(0.0);
    if (params.FetchValueReal(_T("AllocateBuffer"), fBufferSize))
    {
        if (fBufferSize > F(0.1) && fBufferSize < F(32.0))
        {
            m_pBuffer->Initial(static_cast<FLOAT>(fBufferSize));
        }
    }

    InitialLatticeAndConstant(params);
    InitialRandom(params);
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("LatticeIndex")))
    {
        CreateIndexAndBoundary(params);
    }
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("Gauge")))
    {
        CParameters gauge = params.GetParameter(_T("Gauge"));
        CreateGaugeField(gauge);
    }
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("GaugeBoundary")))
    {
        CParameters gaugeboundary = params.GetParameter(_T("GaugeBoundary"));
        CreateGaugeBoundaryField(gaugeboundary);
    }
    checkCudaErrors(cudaGetLastError());
    if (m_InitialCache.constIntegers[ECI_FermionFieldLength] > 0)
    {
        for (UINT i = 1; i <= m_InitialCache.constIntegers[ECI_FermionFieldLength]; ++i)
        {
            CCString sFermionSubParamName;
            sFermionSubParamName.Format(_T("FermionField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters fermionField = params.GetParameter(sFermionSubParamName);
                CreateFermionFields(fermionField);
            }
            checkCudaErrors(cudaGetLastError());
            sFermionSubParamName.Format(_T("BoundaryFermionField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters bcfermionField = params.GetParameter(sFermionSubParamName);
                CreateFermionBoundaryField(bcfermionField);
            }
            checkCudaErrors(cudaGetLastError());
        }
    }
    checkCudaErrors(cudaGetLastError());
    if (m_InitialCache.constIntegers[ECI_ActionListLength] > 0)
    {
        CreateActionList(params);
    }
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("Solver")))
    {
        CParameters solver = params.GetParameter(_T("Solver"));
        CreateSolver(solver);
    }
    for (INT i = 0; i < _kMaxFieldCount; ++i)
    {
        CCString sSolverName = _T("Solver") + appIntToString(i);
        if (params.Exist(sSolverName))
        {
            CParameters solver = params.GetParameter(sSolverName);
            CreateSolver(solver);
        }
    }

    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("GaugeSmearing")))
    {
        CParameters gaugesmearing = params.GetParameter(_T("GaugeSmearing"));
        CreateGaugeSmearing(gaugesmearing);
    }
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("GaugeFixing")))
    {
        CParameters gaugesmearing = params.GetParameter(_T("GaugeFixing"));
        CreateGaugeFixing(gaugesmearing);
    }
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("Updator")))
    {
        CParameters updator = params.GetParameter(_T("Updator"));
        CreateUpdator(updator);
    }
    checkCudaErrors(cudaGetLastError());
    if (m_InitialCache.constIntegers[ECI_MeasureListLength] > 0)
    {
        CreateMeasurement(params);
    }
    checkCudaErrors(cudaGetLastError());
    //=============================================
    // at last, fill the field pointers
    // and copy the index data to device
    InitialIndexBuffer();
    m_pCudaHelper->SetFieldPointers();
    m_pLatticeData->FixAllFieldBoundary();

    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    appGeneral(_T("\n =========== Initialized ! ==============\n"));
    return TRUE;
}

void CCLGLibManager::Quit()
{
    appSafeDelete(m_pLatticeData);
    appSafeDelete(m_pCudaHelper);
    appSafeDelete(m_pFileSystem);
    appSafeDelete(m_pBuffer);

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

void CLGAPI appFailQuitCLG()
{
    GCLGManager.Quit();
    _FAIL_EXIT;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================