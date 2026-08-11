//=============================================================================
// FILENAME : CLGLibMananger.cpp
// 
// DESCRIPTION:
// This is the class for global start-up, control, shut-down
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Update/CStapleCache.h"

#define __CheckTag(tagname,...) iTag = 0; if (params.FetchValueINT(tagname, iTag) && (0 != iTag)) {__VA_ARGS__;}

#define __FetchIntWithDefaultSub(paramname, tagname, defaultv) if (!paramname.FetchValueINT(tagname, iVaules)) { iVaules = defaultv; }

#define __FetchRealWithDefaultSub(paramname, tagname, defaultv) if (!paramname.FetchValueReal(tagname, fValues)) { fValues = defaultv; }

#define __FetchIntWithDefault(tagname, defaultv) __FetchIntWithDefaultSub(params, tagname, defaultv)

#define __FetchRealWithDefault(tagname, defaultv) __FetchRealWithDefaultSub(params, tagname, defaultv)

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
    //__CheckTag(_T("ShowDeviceInformation"), CCudaHelper::DeviceQuery());

    appGeneral(_T("============================== Log Start =============================\n\n"));
}

#pragma region Creates

void CCLGLibManager::InitialLatticeAndConstant(CParameters& params)
{
    //Improve-1: every (re)bake of the lattice/process-grid/halo constants is a
    //new layout generation; handles Bound earlier snapshot the old one and
    //their halos must not be reused (multi-GPU-improve1.md 3.1).
    ++m_ullLayoutGeneration;

    INT iVaules = 0;
    Real fValues = F(0.0);
    CCString sValues;

#pragma region Lattice Size and Threads

    __FetchIntWithDefault(_T("Dim"), 4);
    appAssert(iVaules > 1 && iVaules < 5);
    m_InitialCache.constIntegers[ECI_Dim] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Dir"), 4);
    appAssert(iVaules > 1);
    m_InitialCache.constIntegers[ECI_Dir] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("MaxThreadPerBlock"), _CLG_LAUNCH_MAX_THREAD);
    if (iVaules > 0)
    {
        CCommonData::m_uiMaxThreadPerBlock = iVaules;
    }

    TArray<INT> intValues;

#if _CLG_DEBUG
    if (!params.FetchValueArrayINT(_T("LatticeLengthDebug"), intValues))
    {
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
#else
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
#endif

    //==================================
    // multi-GPU: read and validate the process grid against the GLOBAL lattice.
    // Phase 0 only validates and records the topology; the constants below are
    // still the global lattice. Switching them to per-rank sub-lattice sizes is
    // Phase 1 (see Docs/MultiGPU-Plan.md).
    //==================================
    {
        UINT uiGlobalLattice[4] = {
            static_cast<UINT>(intValues[0]), static_cast<UINT>(intValues[1]),
            static_cast<UINT>(intValues[2]), static_cast<UINT>(intValues[3]) };
        UINT uiGpuGrid[4] = { 1, 1, 1, 1 };

        TArray<INT> gridValues;
        if (params.FetchValueArrayINT(_T("GpuGrid"), gridValues) && gridValues.Num() >= 4)
        {
            for (UINT i = 0; i < 4; ++i)
            {
                uiGpuGrid[i] = static_cast<UINT>(gridValues[i] < 1 ? 1 : gridValues[i]);
            }
        }

        //Halo width: widest stencil reach in one direction. Naik (HISQ) needs 2,
        //so be conservative here rather than guessing small (section 8.3).
        UINT uiHaloWidth = 2;
        INT iHaloWidth = 0;
        if (params.FetchValueINT(_T("HaloWidth"), iHaloWidth))
        {
            //Improve-1 (multi-GPU-improve1.md 3.8): an explicitly configured
            //width must be inside [1, kCacheIndexEdge] -- the boundary/index
            //bakes cover at most kCacheIndexEdge layers, so a wider halo could
            //never be filled. Out-of-range used to fall back to the default
            //silently, which hides a mis-configured run; fail instead. (An
            //absent key still keeps the default above.)
            if (iHaloWidth < 1 || iHaloWidth > CIndexData::kCacheIndexEdge)
            {
                appCrucial(_T("HaloWidth %d is out of range [1, %d]; the index bake cannot represent such a halo.\n"),
                    iHaloWidth, static_cast<INT>(CIndexData::kCacheIndexEdge));
                _FAIL_EXIT;
            }
            uiHaloWidth = static_cast<UINT>(iHaloWidth);
        }

        if (NULL != m_pComm && !m_pComm->SetGpuGrid(uiGpuGrid, uiGlobalLattice, uiHaloWidth))
        {
            //SetGpuGrid already reported which hard constraint failed. An invalid
            //decomposition would silently produce wrong physics, so stop here
            //rather than continue. (This function returns void.)
            appCrucial(_T("GpuGrid validation failed, cannot continue.\n"));
            _FAIL_EXIT;
        }

        //Phase 1: remember the GLOBAL lattice, then shrink intValues to this
        //rank's SUB-LATTICE. Everything below (ECI_Lx, volumes, strides, thread
        //decomposition, field allocation) is then automatically per-rank.
        //The global sizes stay available for position-dependent physics and for
        //gather/scatter -- see Docs/MultiGPU-Plan.md section 1.4-R1.
        for (UINT i = 0; i < 4; ++i)
        {
            m_InitialCache.constIntegers[ECI_GlobalLx + i] = uiGlobalLattice[i];
            m_InitialCache.constIntegers[ECI_GpuGridX + i] = uiGpuGrid[i];
            m_InitialCache.constIntegers[ECI_GlobalOffsetX + i] =
                (NULL == m_pComm) ? 0 : m_pComm->GlobalOffset()[i];
            intValues[i] = static_cast<INT>(uiGlobalLattice[i] / uiGpuGrid[i]);
        }

        //Improve-1 (multi-GPU-improve1.md 3.8): on a split direction the halo
        //is the ONLY source of cross-rank neighbours, so the sub-lattice must
        //be at least HaloWidth long there; a thinner slice would need the same
        //neighbour rank on both sides (or skip a rank entirely) and the
        //face-exchange bookkeeping breaks.
        for (UINT i = 0; i < 4; ++i)
        {
            if (uiGpuGrid[i] > 1 && uiGlobalLattice[i] / uiGpuGrid[i] < uiHaloWidth)
            {
                appCrucial(_T("Split direction %d: local length %d (= global %d / grid %d) is smaller than HaloWidth %d.\n"),
                    static_cast<INT>(i), static_cast<INT>(uiGlobalLattice[i] / uiGpuGrid[i]),
                    static_cast<INT>(uiGlobalLattice[i]), static_cast<INT>(uiGpuGrid[i]),
                    static_cast<INT>(uiHaloWidth));
                _FAIL_EXIT;
            }
        }

        m_InitialCache.constIntegers[ECI_HaloWidth] = uiHaloWidth;
    }

    m_InitialCache.constIntegers[ECI_Lx] = static_cast<UINT>(intValues[0]);
    m_InitialCache.constIntegers[ECI_Ly] = static_cast<UINT>(intValues[1]);
    m_InitialCache.constIntegers[ECI_Lz] = static_cast<UINT>(intValues[2]);
    m_InitialCache.constIntegers[ECI_Lt] = static_cast<UINT>(intValues[3]);
    m_InitialCache.constIntegers[ECI_Volume] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_Volume_xyz] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[2]);
    m_InitialCache.constIntegers[ECI_Volume_xyt] = static_cast<UINT>(intValues[0] * intValues[1] * intValues[3]);
    m_InitialCache.constIntegers[ECI_Volume_xzt] = static_cast<UINT>(intValues[0] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_Volume_yzt] = static_cast<UINT>(intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultX] = static_cast<UINT>(intValues[1] * intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultY] = static_cast<UINT>(intValues[2] * intValues[3]);
    m_InitialCache.constIntegers[ECI_MultZ] = static_cast<UINT>(intValues[3]);
    m_InitialCache.constIntegers[ECI_GridDimZT] = m_InitialCache.constIntegers[ECI_Lz] * m_InitialCache.constIntegers[ECI_Lt];
    TArray<UINT> latticeDim;
    latticeDim.AddItem(intValues[0] * intValues[1]); //xy
    latticeDim.AddItem(intValues[2]); //z
    latticeDim.AddItem(intValues[3]); //t

    if (params.FetchValueArrayINT(_T("Center"), intValues))
    {
        if (4 == intValues.Num())
        {
            m_InitialCache.constSignedIntegers[ECSI_CenterX] = intValues[0];
            m_InitialCache.constSignedIntegers[ECSI_CenterY] = intValues[1];
            m_InitialCache.constSignedIntegers[ECSI_CenterZ] = intValues[2];
            m_InitialCache.constSignedIntegers[ECSI_CenterT] = intValues[3];
            SSmallInt4 sCenter(
                static_cast<SCHAR>(intValues[0]),
                static_cast<SCHAR>(intValues[1]),
                static_cast<SCHAR>(intValues[2]),
                static_cast<SCHAR>(intValues[3]));
            m_InitialCache.constIntegers[ECI_Center] = sCenter.m_uiData;
        }
    }
    else
    {
        //Default rotation center = centre of the GLOBAL lattice. Under a
        //multi-GPU decomposition ECI_Lx.. are local lengths, so using them
        //here would give every rank its own local centre (wrong force
        //coefficients for position-dependent actions, e.g. Acceleration g*t).
        m_InitialCache.constSignedIntegers[ECSI_CenterX] = static_cast<INT>(m_InitialCache.constIntegers[ECI_GlobalLx] / 2);;
        m_InitialCache.constSignedIntegers[ECSI_CenterY] = static_cast<INT>(m_InitialCache.constIntegers[ECI_GlobalLy] / 2);;
        m_InitialCache.constSignedIntegers[ECSI_CenterZ] = static_cast<INT>(m_InitialCache.constIntegers[ECI_GlobalLz] / 2);;
        m_InitialCache.constSignedIntegers[ECSI_CenterT] = static_cast<INT>(m_InitialCache.constIntegers[ECI_GlobalLt] / 2);;
        SSmallInt4 sCenter(
            static_cast<SCHAR>(m_InitialCache.constSignedIntegers[ECSI_CenterX]),
            static_cast<SCHAR>(m_InitialCache.constSignedIntegers[ECSI_CenterY]),
            static_cast<SCHAR>(m_InitialCache.constSignedIntegers[ECSI_CenterZ]),
            static_cast<SCHAR>(m_InitialCache.constSignedIntegers[ECSI_CenterT]));
        m_InitialCache.constIntegers[ECI_Center] = sCenter.m_uiData;
    }

    m_InitialCache.constIntegers[ECI_PlaqutteCount] = m_InitialCache.constIntegers[ECI_Volume] * m_InitialCache.constIntegers[ECI_Dir] * (m_InitialCache.constIntegers[ECI_Dir] - 1) / 2;
    m_InitialCache.constIntegers[ECI_LinkCount] = m_InitialCache.constIntegers[ECI_Volume] * m_InitialCache.constIntegers[ECI_Dir];

    __FetchRealWithDefault(_T("GaugeMomentumFactor"), F(1.0));
    m_InitialCache.constFloats[ECF_GaugeMomentumFactor] = fValues;
    
    UBOOL bAutoDecompose = TRUE;
    __FetchIntWithDefault(_T("ThreadAutoDecompose"), 1);

    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(m_iDeviceId);

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

    const UINT constrainx = min(deviceConstraints[0], deviceConstraints[1]);
    const UINT threadsneed = m_InitialCache.constIntegers[ECI_Volume];
    const UINT cdir = m_InitialCache.constIntegers[ECI_Dir];
    UINT uib = threadsneed > constrainx ? appCeil(threadsneed, constrainx) : 1;
    UINT uit = threadsneed > constrainx ? appCeil(threadsneed, uib) : threadsneed;
    UINT uitd = cdir * (constrainx / cdir);
    UINT uibd = (threadsneed * cdir + uitd - 1) / uitd;
    m_InitialCache.constIntegers[ECI_DecompAllBlock] = uib;
    m_InitialCache.constIntegers[ECI_DecompAllThread] = uit;
    m_InitialCache.constIntegers[ECI_DecompAllBlockDir] = uibd;
    m_InitialCache.constIntegers[ECI_DecompAllThreadDir] = uitd;

    const UINT threadsneedhalf = m_InitialCache.constIntegers[ECI_Volume] / 2;
    m_InitialCache.constIntegers[ECI_VolumeHalf] = threadsneedhalf;
    uib = threadsneedhalf > constrainx ? appCeil(threadsneedhalf, constrainx) : 1;
    uit = threadsneedhalf > constrainx ? appCeil(threadsneedhalf, uib) : threadsneedhalf;
    uitd = cdir * (constrainx / cdir);
    uibd = (threadsneedhalf * cdir + uitd - 1) / uitd;
    m_InitialCache.constIntegers[ECI_DecompAllBlockHalf] = uib;
    m_InitialCache.constIntegers[ECI_DecompAllThreadHalf] = uit;
    m_InitialCache.constIntegers[ECI_DecompAllBlockDirHalf] = uibd;
    m_InitialCache.constIntegers[ECI_DecompAllThreadDirHalf] = uitd;

    TArray<UINT> latticeDim2;
    latticeDim2.AddItem(m_InitialCache.constIntegers[ECI_Lx]);
    latticeDim2.AddItem(m_InitialCache.constIntegers[ECI_Ly]);
    latticeDim2.AddItem(m_InitialCache.constIntegers[ECI_Lz]);
    TArray <UINT> decomp2 = _getDecompose(deviceConstraints, latticeDim2);
    m_InitialCache.constIntegers[ECI_DecompX3D] = decomp2[0];
    m_InitialCache.constIntegers[ECI_DecompY3D] = decomp2[1];
    m_InitialCache.constIntegers[ECI_DecompZ3D] = decomp2[2];
    m_InitialCache.constIntegers[ECI_DecompLx3D] = decomp2[3];
    m_InitialCache.constIntegers[ECI_DecompLy3D] = decomp2[4];
    m_InitialCache.constIntegers[ECI_DecompLz3D] = decomp2[5];

    TArray<UINT> latticeDim3;
    latticeDim3.AddItem(m_InitialCache.constIntegers[ECI_Lx]);
    latticeDim3.AddItem(m_InitialCache.constIntegers[ECI_Ly]);
    latticeDim3.AddItem(m_InitialCache.constIntegers[ECI_Lt]);
    TArray <UINT> decomp3 = _getDecompose(deviceConstraints, latticeDim3);
    m_InitialCache.constIntegers[ECI_DecompX3DXYT] = decomp3[0];
    m_InitialCache.constIntegers[ECI_DecompY3DXYT] = decomp3[1];
    m_InitialCache.constIntegers[ECI_DecompZ3DXYT] = decomp3[2];
    m_InitialCache.constIntegers[ECI_DecompLx3DXYT] = decomp3[3];
    m_InitialCache.constIntegers[ECI_DecompLy3DXYT] = decomp3[4];
    m_InitialCache.constIntegers[ECI_DecompLz3DXYT] = decomp3[5];

    TArray<UINT> latticeDim4;
    latticeDim4.AddItem(m_InitialCache.constIntegers[ECI_Lx]);
    latticeDim4.AddItem(m_InitialCache.constIntegers[ECI_Lz]);
    latticeDim4.AddItem(m_InitialCache.constIntegers[ECI_Lt]);
    TArray <UINT> decomp4 = _getDecompose(deviceConstraints, latticeDim4);
    m_InitialCache.constIntegers[ECI_DecompX3DXZT] = decomp4[0];
    m_InitialCache.constIntegers[ECI_DecompY3DXZT] = decomp4[1];
    m_InitialCache.constIntegers[ECI_DecompZ3DXZT] = decomp4[2];
    m_InitialCache.constIntegers[ECI_DecompLx3DXZT] = decomp4[3];
    m_InitialCache.constIntegers[ECI_DecompLy3DXZT] = decomp4[4];
    m_InitialCache.constIntegers[ECI_DecompLz3DXZT] = decomp4[5];

    TArray<UINT> latticeDim5;
    latticeDim5.AddItem(m_InitialCache.constIntegers[ECI_Ly]);
    latticeDim5.AddItem(m_InitialCache.constIntegers[ECI_Lz]);
    latticeDim5.AddItem(m_InitialCache.constIntegers[ECI_Lt]);
    TArray <UINT> decomp5 = _getDecompose(deviceConstraints, latticeDim5);
    m_InitialCache.constIntegers[ECI_DecompX3DYZT] = decomp5[0];
    m_InitialCache.constIntegers[ECI_DecompY3DYZT] = decomp5[1];
    m_InitialCache.constIntegers[ECI_DecompZ3DYZT] = decomp5[2];
    m_InitialCache.constIntegers[ECI_DecompLx3DYZT] = decomp5[3];
    m_InitialCache.constIntegers[ECI_DecompLy3DYZT] = decomp5[4];
    m_InitialCache.constIntegers[ECI_DecompLz3DYZT] = decomp5[5];

#pragma endregion

#pragma region Fill constant table

    __FetchIntWithDefault(_T("RandomSeed"), 81192);
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

    __FetchIntWithDefault(_T("ExponentialPrecision"), 8);
#if !_CLG_DOUBLEFLOAT
    if (1 != iVaules && iVaules < 8)
    {
        appWarning(_T("Single point float generally does not support quick exponential.\n You may need to set ExponentialPrecision : n with n >= 8\n"));
    }
#endif
    m_InitialCache.constIntegers[ECI_ExponentPrecision] = static_cast<UINT>(iVaules);

    //__FetchIntWithDefault(_T("CacheStaple"), 0);
    //CCommonData::m_bStoreStaple = (0 != iVaules);

    __FetchIntWithDefault(_T("StochasticGaussian"), 0);
    CCommonData::m_bStochasticGaussian = (0 != iVaules);
    
    __FetchIntWithDefault(_T("CacheSolution"), 1);
    CCommonData::m_bStoreLastSolution = (0 != iVaules);

    __FetchIntWithDefault(_T("ActionListLength"), 0);
    m_InitialCache.constIntegers[ECI_ActionListLength] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("FermionFieldCount"), 0);
    m_InitialCache.constIntegers[ECI_FermionFieldLength] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("GaugeFieldCount"), 1);
    m_InitialCache.constIntegers[ECI_GaugeFieldCount] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("BosonFieldCount"), 0);
    m_InitialCache.constIntegers[ECI_BosonFieldCount] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("Tensor2FieldCount"), 0);
    m_InitialCache.constIntegers[ECI_Tensor2FieldCount] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("OtherGaugeFieldCount"), 0);
    m_InitialCache.constIntegers[ECI_OtherGaugeField] = static_cast<UINT>(iVaules);
    if (iVaules > 0)
    {
        appCrucial(_T("OtherGaugeFieldCount is discarded!\n"));
    }

    __FetchStringWithDefault(_T("RandomType"), _T("ER_Schrage"));
    m_InitialCache.eR = __STRING_TO_ENUM(ERandom, sValues);

    __FetchIntWithDefault(_T("MeasureListLength"), 0);
    m_InitialCache.constIntegers[ECI_MeasureListLength] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("UseLogADefinition"), 0);
#if !_CLG_DOUBLEFLOAT
    if (0 != iVaules)
    {
        appCrucial(_T("Single float does not support Log definition!\n"));
    }
#endif
    m_InitialCache.constIntegers[ECI_UseLogADefinition] = static_cast<UINT>(iVaules);

    const UINT iThreadConstraint = m_InitialCache.constIntegers[ECI_ThreadConstaint];
    __FetchIntWithDefault(_T("SummationDecompose"), iThreadConstraint);
    m_InitialCache.constIntegers[ECI_SummationDecompose] = static_cast<UINT>(iVaules);
    appDetailed(_T("Summation decompose: %d\n"), m_InitialCache.constIntegers[ECI_SummationDecompose]);

    __FetchIntWithDefault(_T("Profiler"), 0);
    m_InitialCache.constIntegers[ECI_Profiler] = static_cast<UINT>(iVaules);

    __FetchIntWithDefault(_T("MILC_StaggeredPhase"), 0);
    m_InitialCache.constIntegers[ECI_MILC_StaggeredPhase] = static_cast<UINT>(iVaules);

    memcpy(m_pCudaHelper->m_ConstIntegers, m_InitialCache.constIntegers, sizeof(UINT) * kContentLength);
    memcpy(m_pCudaHelper->m_ConstSignedIntegers, m_InitialCache.constSignedIntegers, sizeof(INT)* kContentLength);
    memcpy(m_pCudaHelper->m_ConstFloats, m_InitialCache.constFloats, sizeof(Real) * kContentLength);
    m_pCudaHelper->CopyConstants();
    m_pCudaHelper->AllocateTemeraryBuffers(_HC_Volume);

#pragma endregion

    m_pCudaHelper->CreateGammaMatrix();
    _CHECKCUDA;
}

void CCLGLibManager::InitialRandom(CParameters &)
{
    //INT iVaules = 0;
    //CCString sValues;

    m_pLatticeData->m_pRandom = new CRandom(m_InitialCache.constIntegers[ECI_RandomSeed], m_InitialCache.eR);
    checkCudaErrors(__cudaMalloc((void**)&(m_pLatticeData->m_pDeviceRandom), sizeof(CRandom)));
    checkCudaErrors(cudaMemcpy(m_pLatticeData->m_pDeviceRandom, m_pLatticeData->m_pRandom, sizeof(CRandom), cudaMemcpyHostToDevice));
    appGeneral(_T("Create the %s random with seed:%d\n"), __ENUM_TO_STRING(ERandom, m_InitialCache.eR).c_str(), m_InitialCache.constIntegers[ECI_RandomSeed]);

    m_pCudaHelper->CopyRandomPointer(m_pLatticeData->m_pDeviceRandom);
    m_pLatticeData->m_uiRandomType = static_cast<UINT>(m_InitialCache.eR);
    m_pLatticeData->m_uiRandomSeed = m_InitialCache.constIntegers[ECI_RandomSeed];
}

CField* CCLGLibManager::CreateGaugeFields(class CParameters& params) const
{
    INT iVaules = 0;
    CCString sValues;

    CCString sGaugeClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldGaugeSU3"));
    sGaugeClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    const EFieldInitialType eGaugeInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    __FetchIntWithDefault(_T("FieldId"), 1);
    BYTE byFieldId = static_cast<BYTE>(iVaules);
    if (m_pLatticeData->m_pFieldMap.Exist(byFieldId))
    {
        appCrucial(_T("Unable to create the gauge field! with wrong field ID %s %d!, using default: %d\n"), sGaugeClassName.c_str(), byFieldId, m_byLoadingFieldId);
        byFieldId = m_byLoadingFieldId;
    }

    CBase* pGaugeField = appCreate(sGaugeClassName);
    CFieldGauge* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldGauge*>(pGaugeField)) : NULL;
    if (NULL == pGauge)
    {
        appCrucial(_T("Unable to create the gauge field! with name %s!\n"), sGaugeClassName.c_str());
        return NULL;
    }

    pGauge->m_byFieldId = byFieldId;
    pGauge->m_pOwner = m_pLatticeData;
    //Improve-1: keep the halo handle's tag id in sync (identity unchanged).
    if (NULL != pGauge->GetHaloBufferHandle()) { pGauge->GetHaloBufferHandle()->SetFieldId(byFieldId); }
    if (EFIT_ReadFromFile == eGaugeInitial)
    {
        CCString sFileType, sFileName;
        if (!params.FetchStringValue(_T("GaugeFileType"), sFileType)
            || !params.FetchStringValue(_T("GaugeFileName"), sFileName))
        {
            appCrucial(_T("Gauge initial type is EFIT_ReadFromFile, but cannot find GaugeFileType or GaugeFileName!\n"));
            _FAIL_EXIT;
        }
        const EFieldFileType eFileType = __STRING_TO_ENUM(EFieldFileType, sFileType);
        pGauge->InitialFieldWithFile(sFileName, eFileType);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    pGauge->InitialOtherParameters(params);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SCHAR>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SCHAR>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SCHAR>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SCHAR>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }
    else
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = 1;
        bc.m_sPeriodic.y = 1;
        bc.m_sPeriodic.z = 1;
        bc.m_sPeriodic.w = 1;
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    m_pLatticeData->m_pGaugeField.AddItem(pGauge);
    m_pLatticeData->m_pFieldMap.SetAt(byFieldId, pGauge);
    m_pLatticeData->m_pOtherFields.AddItem(pGauge);
    m_pLatticeData->m_eFieldInitialTypes.AddItem(eGaugeInitial);
    //m_pLatticeData->CreateFieldPool(byFieldId, 0);
    checkCudaErrors(cudaDeviceSynchronize());
    appGeneral(_T("Create the gauge %s (field id %d) with initial: %s\n"), sGaugeClassName.c_str(), byFieldId, sValues.c_str());
    return pGauge;
}

CField* CCLGLibManager::CreateBosonFields(class CParameters& params) const
{
    INT iVaules = 0;
    CCString sValues;

    CCString sBosonClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldBosonU1"));
    sBosonClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    const EFieldInitialType eFieldInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);
    checkCudaErrors(cudaDeviceSynchronize());
    CBase* pBosonField = appCreate(sBosonClassName);
    CFieldBoson* pBoson = (NULL != pBosonField) ? (dynamic_cast<CFieldBoson*>(pBosonField)) : NULL;
    if (NULL == pBoson)
    {
        appCrucial(_T("Unable to create the boson field! with name %s!"), sBosonClassName.c_str());
        return NULL;
    }

    __FetchIntWithDefault(_T("FieldId"), -1);
    BYTE byFieldId = static_cast<BYTE>(iVaules);
    if (byFieldId >= kMaxFieldCount || byFieldId <= 1 || m_pLatticeData->m_pFieldMap.Exist(byFieldId))
    {
        appCrucial(_T("Unable to create the boson field! with wrong field ID %s %d! Using default: %d\n"), sBosonClassName.c_str(), byFieldId, m_byLoadingFieldId);
        byFieldId = m_byLoadingFieldId;
    }

    pBoson->m_byFieldId = byFieldId;
    pBoson->m_pOwner = m_pLatticeData;
    //Improve-1: keep the halo handle's tag id in sync (identity unchanged).
    if (NULL != pBoson->GetHaloBufferHandle()) { pBoson->GetHaloBufferHandle()->SetFieldId(byFieldId); }
    //checkCudaErrors(cudaDeviceSynchronize());
    //pBoson->InitialField(eFieldInitial);
    //checkCudaErrors(cudaDeviceSynchronize());
    pBoson->InitialOtherParameters(params);
    m_pLatticeData->m_pFieldMap.SetAt(byFieldId, pBoson);
    m_pLatticeData->m_pBosonField.AddItem(pBoson);
    m_pLatticeData->m_pOtherFields.AddItem(pBoson);
    m_pLatticeData->m_eFieldInitialTypes.AddItem(eFieldInitial);
    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SCHAR>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SCHAR>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SCHAR>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SCHAR>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }
    else
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = 1;
        bc.m_sPeriodic.y = 1;
        bc.m_sPeriodic.z = 1;
        bc.m_sPeriodic.w = 1;
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //__FetchIntWithDefault(_T("PoolNumber"), 0);
    //if (iVaules > 0)
    //{
    //m_pLatticeData->CreateFieldPool(byFieldId, iVaules);
    checkCudaErrors(cudaDeviceSynchronize());
    //}
    appGeneral(_T("Create the boson field %s with id %d and initial: %s\n"), sBosonClassName.c_str(), byFieldId, sValues.c_str());
    return pBoson;
}

CField* CCLGLibManager::CreateTensor2Fields(class CParameters& params) const
{
    INT iVaules = 0;
    CCString sValues;

    CCString sTensor2ClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldTensor2SU3"));
    sTensor2ClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    const EFieldInitialType eFieldInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);
    checkCudaErrors(cudaDeviceSynchronize());
    CBase* pTensor2Field = appCreate(sTensor2ClassName);
    CFieldTensor2* pTensor2 = (NULL != pTensor2Field) ? (dynamic_cast<CFieldTensor2*>(pTensor2Field)) : NULL;
    if (NULL == pTensor2)
    {
        appCrucial(_T("Unable to create the tensor2 field! with name %s!"), sTensor2ClassName.c_str());
        return NULL;
    }

    __FetchIntWithDefault(_T("FieldId"), -1);
    BYTE byFieldId = static_cast<BYTE>(iVaules);
    if (byFieldId >= kMaxFieldCount || byFieldId <= 1 || m_pLatticeData->m_pFieldMap.Exist(byFieldId))
    {
        appCrucial(_T("Unable to create the tensor2 field! with wrong field ID %s %d! Using default: %d\n"), sTensor2ClassName.c_str(), byFieldId, m_byLoadingFieldId);
        byFieldId = m_byLoadingFieldId;
    }

    pTensor2->m_byFieldId = byFieldId;
    pTensor2->m_pOwner = m_pLatticeData;
    //Improve-1: keep the halo handle's tag id in sync (identity unchanged).
    if (NULL != pTensor2->GetHaloBufferHandle()) { pTensor2->GetHaloBufferHandle()->SetFieldId(byFieldId); }
    pTensor2->InitialOtherParameters(params);
    m_pLatticeData->m_pFieldMap.SetAt(byFieldId, pTensor2);
    m_pLatticeData->m_pOtherFields.AddItem(pTensor2);
    m_pLatticeData->m_pTensor2Field.AddItem(pTensor2);
    m_pLatticeData->m_eFieldInitialTypes.AddItem(eFieldInitial);
    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SCHAR>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SCHAR>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SCHAR>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SCHAR>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }
    else
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = 1;
        bc.m_sPeriodic.y = 1;
        bc.m_sPeriodic.z = 1;
        bc.m_sPeriodic.w = 1;
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    appGeneral(_T("Create the tensor2 field %s with id %d and initial: %s\n"), sTensor2ClassName.c_str(), byFieldId, sValues.c_str());
    return pTensor2;
}

CField* CCLGLibManager::CreateFermionFields(class CParameters& params) const
{
    INT iVaules = 0;
    CCString sValues;

    CCString sFermionClassName;
    __FetchStringWithDefault(_T("FieldName"), _T("CFieldFermionWilsonSquareSU3"));
    sFermionClassName = sValues;
    __FetchStringWithDefault(_T("FieldInitialType"), _T("EFIT_Random"));
    const EFieldInitialType eFieldInitial = __STRING_TO_ENUM(EFieldInitialType, sValues);

    CBase* pFermionField = appCreate(sFermionClassName);
    CFieldFermion* pFermion = (NULL != pFermionField) ? (dynamic_cast<CFieldFermion*>(pFermionField)) : NULL;
    if (NULL == pFermion)
    {
        appCrucial(_T("Unable to create the fermion field! with name %s!"), sFermionClassName.c_str());
        return NULL;
    }

    __FetchIntWithDefault(_T("FieldId"), -1);
    BYTE byFieldId = static_cast<BYTE>(iVaules);
    if (byFieldId >= kMaxFieldCount || byFieldId <= 1 || m_pLatticeData->m_pFieldMap.Exist(byFieldId))
    {
        appCrucial(_T("Unable to create the fermion field! with wrong field ID %s %d! Using default: %d\n"), sFermionClassName.c_str(), byFieldId, m_byLoadingFieldId);
        byFieldId = m_byLoadingFieldId;
    }

    pFermion->m_byFieldId = byFieldId;
    pFermion->m_pOwner = m_pLatticeData;
    //Improve-1: keep the halo handle's tag id in sync (identity unchanged).
    if (NULL != pFermion->GetHaloBufferHandle()) { pFermion->GetHaloBufferHandle()->SetFieldId(byFieldId); }
    //pFermion->InitialField(eFieldInitial);
    pFermion->InitialOtherParameters(params);
    m_pLatticeData->m_pFieldMap.SetAt(byFieldId, pFermion);
    m_pLatticeData->m_pFermionField.AddItem(pFermion);
    m_pLatticeData->m_pOtherFields.AddItem(pFermion);
    m_pLatticeData->m_eFieldInitialTypes.AddItem(eFieldInitial);
    TArray<INT> periodic;
    if (params.FetchValueArrayINT(_T("Period"), periodic))
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = static_cast<SCHAR>(periodic[0]);
        bc.m_sPeriodic.y = static_cast<SCHAR>(periodic[1]);
        bc.m_sPeriodic.z = static_cast<SCHAR>(periodic[2]);
        bc.m_sPeriodic.w = static_cast<SCHAR>(periodic[3]);
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }
    else
    {
        SBoundCondition bc;
        bc.m_sPeriodic.x = 1;
        bc.m_sPeriodic.y = 1;
        bc.m_sPeriodic.z = 1;
        bc.m_sPeriodic.w = -1;
        m_pLatticeData->SetFieldBoundaryCondition(byFieldId, bc);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    //__FetchIntWithDefault(_T("PoolNumber"), 0);
    //if (iVaules > 0)
    //{
        //m_pLatticeData->CreateFieldPool(byFieldId, iVaules);
    //}
    checkCudaErrors(cudaDeviceSynchronize());
    appGeneral(_T("Create the fermion field %s with id %d and initial: %s\n"), sFermionClassName.c_str(), byFieldId, sValues.c_str());
    return pFermion;
}

void CCLGLibManager::CreateBoundaryFields(class CParameters& params, const CCString& sDefaultName) const
{
    CCString sValues;
    INT iVaules;
    //__FetchStringWithDefault(_T("FieldName"), _T("CFieldBoundaryWilsonSquareSU3"));
    //__FetchStringWithDefault(_T("FieldName"), _T("CFieldBoundaryGaugeSU3"));
    __FetchStringWithDefault(_T("FieldName"), sDefaultName);
    const CCString sFieldClassName = sValues;

    CBase* pBCField = appCreate(sFieldClassName);
    CFieldBoundaryParent* pBC = (NULL != pBCField) ? (dynamic_cast<CFieldBoundaryParent*>(pBCField)) : NULL;

    if (NULL == pBC)
    {
        appCrucial(_T("Unable to create the boundary field! with name %s!"), sFieldClassName.c_str());
    }
    pBC->InitialField(params);
    INT defaultId = -1;
    if (params.GetName() == _T("GaugeBoundary"))
    {
        defaultId = 1;
    }
    __FetchIntWithDefault(_T("FieldId"), defaultId);
    const BYTE byFieldId = static_cast<BYTE>(iVaules);
    m_pLatticeData->m_pAllBoundaryFields.AddItem(pBC);
    m_pLatticeData->m_pBoundaryFieldMap.SetAt(byFieldId, pBC);
    checkCudaErrors(cudaDeviceSynchronize());
    appGeneral(_T("Create the boundary field %s with initial: %s\n"), sFieldClassName.c_str(), sValues.c_str());
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
    else if (NULL != updator && EUT_Heatbath == updator->GetUpdatorType())
    {
        updator->Initial(m_pLatticeData, params);
        m_pLatticeData->m_pUpdator = updator;
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

void CCLGLibManager::CreateMultiShiftSolver(class CParameters& params) const
{
    CCString sSolverName = _T("CMultiShiftGMRES");
    params.FetchStringValue(_T("SolverName"), sSolverName);
    INT byFieldId = 2;
    params.FetchValueINT(_T("SolverForFieldId"), byFieldId);
    CField* pField = m_pLatticeData->GetFieldById(static_cast<BYTE>(byFieldId));
    if (NULL == pField)
    {
        appCrucial(_T("Solver must be created for a specified field!\n"));
    }
    m_pLatticeData->CreateMultiShiftSolver(sSolverName, params, pField, static_cast<BYTE>(byFieldId));
}

void CCLGLibManager::CreateGaugeSmearing(class CParameters& params) const
{
    CCString sSmearingName = _T("CGaugeSmearingAPEStout");
    params.FetchStringValue(_T("SmearingName"), sSmearingName);
    CGaugeSmearing* pSmearing = dynamic_cast<CGaugeSmearing*>(appCreate(sSmearingName));
    if (NULL != pSmearing)
    {
        pSmearing->Initial(m_pLatticeData, params);
    }
}

void CCLGLibManager::CreateGaugeStapleCache(class CParameters& params) const
{
    CCString sStapleCacheName = _T("CStapleCacheSU3");
    params.FetchStringValue(_T("StapleCacheName"), sStapleCacheName);
    CStapleCache* pCache = dynamic_cast<CStapleCache*>(appCreate(sStapleCacheName));
    if (NULL != pCache)
    {
        pCache->Initial(m_pLatticeData, params);
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

void CCLGLibManager::InitialFieldBuffer() const
{
    appAssert(m_pLatticeData->m_eFieldInitialTypes.Num() == m_pLatticeData->m_pOtherFields.Num());
    for (INT i = 0; i < m_pLatticeData->m_pOtherFields.Num(); ++i)
    {
        if (EFIT_ReadFromFile != m_pLatticeData->m_eFieldInitialTypes[i])
        {
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
            m_pLatticeData->m_pOtherFields[i]->InitialField(m_pLatticeData->m_eFieldInitialTypes[i]);
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
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
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    m_pLatticeData->m_pIndex->BakeAllIndexBuffer(m_pLatticeData->m_pIndexCache);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    if (m_pLatticeData->m_pOtherFields.Num() > 0)
    {
        //UBOOL bHasStaggeredFermion = FALSE;
        //UBOOL bPlaqCached = FALSE;
        for (BYTE i = 1; i < kMaxFieldCount; ++i)
        {
            const CField* pf = m_pLatticeData->GetFieldById(i);
            //if (NULL != pf && pf->IsGaugeField() && !bPlaqCached)
            //{
            //    m_pLatticeData->m_pIndex->BakePlaquttes(m_pLatticeData->m_pIndexCache, i);
            //    checkCudaErrors(cudaDeviceSynchronize());
            //    checkCudaErrors(cudaGetLastError());
            //    bPlaqCached = TRUE;
            //}
            //else 
            if (NULL != pf && pf->IsGaugeField())
            {
                m_pLatticeData->m_pIndex->BakePlaquttes(m_pLatticeData->m_pIndexCache, i);
                checkCudaErrors(cudaDeviceSynchronize());
                checkCudaErrors(cudaGetLastError());
            }

            //move index is only used by the hopping of boson and fermion fields
            if (NULL != pf && (pf->IsBosonField() || pf->IsFermionField()))
            {
                m_pLatticeData->m_pIndex->BakeMoveIndex(m_pLatticeData->m_pIndexCache, i);
                checkCudaErrors(cudaDeviceSynchronize());
                checkCudaErrors(cudaGetLastError());
                //if (NULL != dynamic_cast<const CFieldFermionKS*>(pf))
                //{
                //    bHasStaggeredFermion = TRUE;
                //}

                if (NULL != dynamic_cast<const CFieldFermionKS*>(pf))
                {
                    m_pLatticeData->m_pIndex->BakeNaikTable(m_pLatticeData->m_pIndexCache, i);
                    checkCudaErrors(cudaDeviceSynchronize());
                    checkCudaErrors(cudaGetLastError());
                }
            }
        }
        //if (bHasStaggeredFermion)
        {
            //always bake it, so that we have even-odd as eta_5
            m_pLatticeData->m_pIndex->BakeEtaMuTable(m_pLatticeData->m_pIndexCache);
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
    }

    //m_pLatticeData->m_pIndex->BakeEvenOddTable(m_pLatticeData->m_pIndexCache);
    //checkCudaErrors(cudaDeviceSynchronize());
    //checkCudaErrors(cudaGetLastError());

    m_pLatticeData->m_pIndex->CalculateSiteCount(m_pLatticeData->m_pIndexCache);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    m_pCudaHelper->SetDeviceIndex(m_pLatticeData->m_pIndexCache);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

#pragma endregion

UBOOL CCLGLibManager::InitialWithParameter(CParameters &params)
{
    //==================================
    // multi-GPU: MPI must come up before device selection, because the rank
    // decides which device this process takes. On a single-GPU build this is a
    // lone rank with a [1,1,1,1] grid, so nothing below changes behaviour.
    //==================================
    m_pComm = new CLGComm();
    m_pHaloManager = new CHaloManager();
    m_pComm->Initial(NULL, NULL);

    //==================================
    // set device
    //==================================
    m_iDeviceId = 0;
    m_byLoadingFieldId = 1;
    UBOOL bDeviceSet = params.FetchValueINT(_T("DeviceIndex"), m_iDeviceId);

#if _CLG_MULTI_GPU
    if (m_pComm->Size() > 1)
    {
        //One process per GPU. DevicePerNode lets a single-GPU test machine
        //oversubscribe (all ranks share device 0) for 1-vs-N comparison runs --
        //see Docs/MultiGPU-Plan.md section 4.5.
        INT iDevicePerNode = 0;
        if (!params.FetchValueINT(_T("DevicePerNode"), iDevicePerNode) || iDevicePerNode < 1)
        {
            INT iDeviceCount = 1;
            if (cudaSuccess != cudaGetDeviceCount(&iDeviceCount) || iDeviceCount < 1)
            {
                iDeviceCount = 1;
            }
            iDevicePerNode = iDeviceCount;
        }
        m_iDeviceId = static_cast<INT>(m_pComm->Rank() % static_cast<UINT>(iDevicePerNode));
        bDeviceSet = TRUE;
    }
#endif

    if (bDeviceSet)
    {
        checkCudaErrors(cudaSetDevice(m_iDeviceId));
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
        //checkCudaErrors(cudaInitDevice(m_iDeviceId));

        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, m_iDeviceId);
        appGeneral("\nRank %d uses Device %d: \"%s\"\n",
            static_cast<INT>(m_pComm->Rank()), m_iDeviceId, deviceProp.name);
    }

    //==================================

    m_pCudaHelper = new CCudaHelper();
    m_pLatticeData = new CLatticeData();
    m_pFileSystem = new CFileSystem();
    m_pBuffer = new CCudaBuffer();
    m_pFieldPool = new CFieldPool();
    //UBOOL bGaugeBoundaryFieldCreated = FALSE;
    //Allocate Buffer
    Real fBufferSize = F(0.0);
    if (params.FetchValueReal(_T("AllocateBuffer"), fBufferSize))
    {
        if (fBufferSize > F(0.1) && fBufferSize < F(32.0))
        {
            m_pBuffer->Initial(static_cast<FLOAT>(fBufferSize));
            _CHECKCUDA;
        }
    }

    InitialLatticeAndConstant(params);
    _CHECKCUDA;
    InitialRandom(params);
    _CHECKCUDA;
    if (params.Exist(_T("LatticeIndex")))
    {
        CreateIndexAndBoundary(params);
    }
    _CHECKCUDA;
    if (params.Exist(_T("Gauge")))
    {
        CParameters gauge = params.GetParameter(_T("Gauge"));
        const CField* pLoaded = CreateGaugeFields(gauge);
        if (NULL != pLoaded)
        {
            m_byLoadingFieldId = pLoaded->m_byFieldId + 1;
        }
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("GaugeBoundary")))
    {
        CParameters gaugeboundary = params.GetParameter(_T("GaugeBoundary"));
        CreateBoundaryFields(gaugeboundary, _T("CFieldBoundaryGaugeSU3"));
        //bGaugeBoundaryFieldCreated = TRUE;
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    if (m_InitialCache.constIntegers[ECI_GaugeFieldCount] > 1)
    {
        for (UINT i = 2; i <= m_InitialCache.constIntegers[ECI_GaugeFieldCount]; ++i)
        {
            CCString sFermionSubParamName;
            sFermionSubParamName.Format(_T("Gauge%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters fermionField = params.GetParameter(sFermionSubParamName);
                const CField* pLoaded = CreateGaugeFields(fermionField);
                if (NULL != pLoaded)
                {
                    m_byLoadingFieldId = pLoaded->m_byFieldId + 1;
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
            sFermionSubParamName.Format(_T("GaugeBoundary%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters bcfermionField = params.GetParameter(sFermionSubParamName);
                CreateBoundaryFields(bcfermionField, _T("CFieldBoundaryGaugeSU3"));
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
    }

    if (m_InitialCache.constIntegers[ECI_BosonFieldCount] > 0)
    {
        for (UINT i = 1; i <= m_InitialCache.constIntegers[ECI_BosonFieldCount]; ++i)
        {
            CCString sFermionSubParamName;
            sFermionSubParamName.Format(_T("BosonField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters fermionField = params.GetParameter(sFermionSubParamName);
                const CField* pLoaded = CreateBosonFields(fermionField);
                if (NULL != pLoaded)
                {
                    m_byLoadingFieldId = pLoaded->m_byFieldId + 1;
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
            sFermionSubParamName.Format(_T("BoundaryBosonField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters bcfermionField = params.GetParameter(sFermionSubParamName);
                CreateBoundaryFields(bcfermionField, _T("CBosonBoundary"));
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
    }

    if (m_InitialCache.constIntegers[ECI_Tensor2FieldCount] > 0)
    {
        for (UINT i = 1; i <= m_InitialCache.constIntegers[ECI_Tensor2FieldCount]; ++i)
        {
            CCString sTensor2SubParamName;
            sTensor2SubParamName.Format(_T("Tensor2Field%d"), i);
            if (params.Exist(sTensor2SubParamName))
            {
                CParameters tensor2Field = params.GetParameter(sTensor2SubParamName);
                const CField* pLoaded = CreateTensor2Fields(tensor2Field);
                if (NULL != pLoaded)
                {
                    m_byLoadingFieldId = pLoaded->m_byFieldId + 1;
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
    }

    if (m_InitialCache.constIntegers[ECI_FermionFieldLength] > 0)
    {
        for (UINT i = 1; i <= m_InitialCache.constIntegers[ECI_FermionFieldLength]; ++i)
        {
            CCString sFermionSubParamName;
            sFermionSubParamName.Format(_T("FermionField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters fermionField = params.GetParameter(sFermionSubParamName);
                const CField* pLoaded = CreateFermionFields(fermionField);
                if (NULL != pLoaded)
                {
                    m_byLoadingFieldId = pLoaded->m_byFieldId + 1;
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
            sFermionSubParamName.Format(_T("BoundaryFermionField%d"), i);
            if (params.Exist(sFermionSubParamName))
            {
                CParameters bcfermionField = params.GetParameter(sFermionSubParamName);
                CreateBoundaryFields(bcfermionField, _T("CFieldBoundaryWilsonSquareSU3"));
            }
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());
        }
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    //=============================================
    // at last, fill the field pointers
    // and copy the index data to device
    InitialIndexBuffer();
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    m_pCudaHelper->SetFieldPointers();
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    m_pLatticeData->FixAllFieldBoundary();
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    //if (NULL != m_pLatticeData->m_pIndex && m_pLatticeData->m_pIndex->NeedToFixBoundary() && !bGaugeBoundaryFieldCreated)
    //{
    //    appCrucial(_T("Using Dirichlet boundary without specify a gauge boundary!\n"));
    //}

    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    //The buffer can be initialled only after boundary fixed
    InitialFieldBuffer();

    //Other things such as measurement, integrator, solver, gauge fixing, gauge smearing...
    if (m_InitialCache.constIntegers[ECI_ActionListLength] > 0)
    {
        CreateActionList(params);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    if (params.Exist(_T("Solver")))
    {
        CParameters solver = params.GetParameter(_T("Solver"));
        CreateSolver(solver);
    }
    for (INT i = 0; i < kMaxFieldCount; ++i)
    {
        CCString sSolverName = _T("Solver") + appToString(i);
        if (params.Exist(sSolverName))
        {
            CParameters solver = params.GetParameter(sSolverName);
            CreateSolver(solver);
        }
    }

    if (params.Exist(_T("MSSolver")))
    {
        CParameters solver = params.GetParameter(_T("MSSolver"));
        CreateMultiShiftSolver(solver);
    }
    for (INT i = 0; i < kMaxFieldCount; ++i)
    {
        CCString sSolverName = _T("MSSolver") + appToString(i);
        if (params.Exist(sSolverName))
        {
            CParameters solver = params.GetParameter(sSolverName);
            CreateMultiShiftSolver(solver);
        }
    }
    _CHECKCUDA;
    if (params.Exist(_T("GaugeSmearing")))
    {
        CParameters gaugesmearing = params.GetParameter(_T("GaugeSmearing"));
        CreateGaugeSmearing(gaugesmearing);
    }
    for (INT i = 0; i < kMaxFieldCount; ++i)
    {
        CCString smearingname = _T("GaugeSmearing") + appToString(i);
        if (params.Exist(smearingname))
        {
            CParameters smearingparam = params.GetParameter(smearingname);
            CreateGaugeSmearing(smearingparam);
        }
    }
    _CHECKCUDA;
    if (params.Exist(_T("StapleCache")))
    {
        CParameters staplecache = params.GetParameter(_T("StapleCache"));
        CreateGaugeStapleCache(staplecache);
    }
    for (INT i = 0; i < kMaxFieldCount; ++i)
    {
        CCString staplecachename = _T("StapleCache") + appToString(i);
        if (params.Exist(staplecachename))
        {
            CParameters staplecache = params.GetParameter(staplecachename);
            CreateGaugeStapleCache(staplecache);
        }
    }
    _CHECKCUDA;
    if (params.Exist(_T("GaugeFixing")))
    {
        CParameters gaugesmearing = params.GetParameter(_T("GaugeFixing"));
        CreateGaugeFixing(gaugesmearing);
    }
    _CHECKCUDA;
    if (params.Exist(_T("Updator")))
    {
        CParameters updator = params.GetParameter(_T("Updator"));
        CreateUpdator(updator);
    }
    _CHECKCUDA;
    if (m_InitialCache.constIntegers[ECI_MeasureListLength] > 0)
    {
        CreateMeasurement(params);
    }
    _CHECKCUDA;

    appGeneral(_T("\n =========== Initialized ! ==============\n"));
    return TRUE;
}

void CCLGLibManager::Quit()
{
    //for gauge field to be return
    //for (INT i = 0; i < kMaxFieldCount; ++i)
    //{
    //    appSafeDelete(appGetLattice()->m_pGaugeSmearing[i]);
    //}

    appPrintAllErrors();

    GRASet.Quit();
    appSafeDelete(m_pLatticeData);
    appSafeDelete(m_pCudaHelper);
    appSafeDelete(m_pFileSystem);
    for (INT i = 0; i < m_lstBufferCaches.Num(); ++i)
    {
        appSafeDelete(m_lstBufferCaches[i]);
    }
    m_lstBufferCaches.RemoveAll();
    appSafeDelete(m_pFieldPool);
    appSafeDelete(m_pBuffer);

    //checkCudaErrors(cudaSetDevice(m_iDeviceId));
    checkCudaErrors(cudaDeviceReset());

    //Multi-GPU: tear down after the device is released. The CLGComm destructor
    //calls MPI_Finalize when needed; on single-GPU builds both are no-ops.
    appSafeDelete(m_pHaloManager);
    appSafeDelete(m_pComm);
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