//=============================================================================
// FILENAME : CLGComm.cpp
//
// DESCRIPTION:
// See CLGComm.h. <mpi.h> is included only under _CLG_MULTI_GPU so that the
// single-GPU configurations build on machines without MPI installed.
//
// REVISION:
//  [07/29/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"

#if _CLG_MULTI_GPU
#include <mpi.h>
#include <cstdlib>
#endif

__BEGIN_NAMESPACE

#if _CLG_MULTI_GPU
// MPI must be initialized exactly once per process and finalized exactly once
// at process exit. The test harness creates/destroys a CLGComm per test case,
// so tying MPI_Finalize to the CLGComm destructor would finalize MPI after the
// first test and leave every subsequent test calling MPI routines on a
// finalized library ("Attempting to use an MPI routine after finalizing MPI").
// Instead we register a single atexit finalizer on first init.
static void _clgFinalizeMPIAtExit()
{
    INT bFinalized = 0;
    MPI_Finalized(&bFinalized);
    if (!bFinalized)
    {
        MPI_Finalize();
    }
}
#endif

CLGComm::CLGComm()
    : m_uiRank(0)
    , m_uiSize(1)
    , m_bInitialed(FALSE)
{
    for (UINT i = 0; i < 4; ++i)
    {
        m_uiGpuGrid[i] = 1;
        m_uiGridCoord[i] = 0;
        m_uiGlobalOffset[i] = 0;
        m_uiGlobalLattice[i] = 1;
        m_uiLocalLattice[i] = 1;
    }
}

CLGComm::~CLGComm()
{
    Finalize();
}

void CLGComm::Initial(INT* argc, TCHAR*** argv)
{
#if _CLG_MULTI_GPU
    INT bAlready = 0;
    MPI_Initialized(&bAlready);
    if (!bAlready)
    {
        MPI_Init(argc, argv);
        //Finalize once, at process exit, not per CLGComm lifetime.
        atexit(_clgFinalizeMPIAtExit);
    }

    INT iRank = 0;
    INT iSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    MPI_Comm_size(MPI_COMM_WORLD, &iSize);
    m_uiRank = static_cast<UINT>(iRank);
    m_uiSize = static_cast<UINT>(iSize);
#else
    //Single GPU: always a lone rank. Parameters intentionally unused.
    (void)argc;
    (void)argv;
    m_uiRank = 0;
    m_uiSize = 1;
#endif
    m_bInitialed = TRUE;
}

void CLGComm::Finalize()
{
    if (!m_bInitialed)
    {
        return;
    }
    m_bInitialed = FALSE;

    //Note: MPI_Finalize is intentionally NOT called here. It is deferred to a
    //single atexit handler registered in Initial(), so that per-test CLGComm
    //objects can be created and destroyed without finalizing MPI. See the
    //comment on _clgFinalizeMPIAtExit above.
}

UBOOL CLGComm::SetGpuGrid(const UINT* pGrid, const UINT* pGlobalLattice, UINT uiHaloWidth)
{
    if (NULL == pGrid || NULL == pGlobalLattice)
    {
        appCrucial(_T("CLGComm::SetGpuGrid: null argument.\n"));
        return FALSE;
    }

    UINT uiProduct = 1;
    for (UINT i = 0; i < 4; ++i)
    {
        if (pGrid[i] < 1)
        {
            appCrucial(_T("CLGComm::SetGpuGrid: GpuGrid[%d] = %d must be >= 1.\n"), i, pGrid[i]);
            return FALSE;
        }
        uiProduct *= pGrid[i];
    }

    //Hard constraint: the grid must use exactly the ranks we were given.
    if (uiProduct != m_uiSize)
    {
        appCrucial(_T("CLGComm::SetGpuGrid: GpuGrid product %d does not match MPI world size %d.\n"),
            uiProduct, m_uiSize);
        return FALSE;
    }

    for (UINT i = 0; i < 4; ++i)
    {
        //Hard constraint: only even splits are supported for now.
        if (0 != (pGlobalLattice[i] % pGrid[i]))
        {
            appCrucial(_T("CLGComm::SetGpuGrid: lattice length %d in direction %d is not divisible by GpuGrid %d.\n"),
                pGlobalLattice[i], i, pGrid[i]);
            return FALSE;
        }

        const UINT uiLocal = pGlobalLattice[i] / pGrid[i];

        //Hard constraint: a split direction must be at least as long as the halo,
        //otherwise a halo would reach past the immediate neighbour rank.
        if (pGrid[i] > 1 && uiLocal < uiHaloWidth)
        {
            appCrucial(_T("CLGComm::SetGpuGrid: local length %d in split direction %d is smaller than halo width %d.\n"),
                uiLocal, i, uiHaloWidth);
            return FALSE;
        }

        m_uiGpuGrid[i] = pGrid[i];
    }

    //Rank <-> grid coordinate mapping. x is the slowest varying direction, matching
    //the lattice site index convention (see Docs/MultiGPU-Plan.md section 9.2).
    UINT uiRest = m_uiRank;
    for (INT i = 3; i >= 0; --i)
    {
        m_uiGridCoord[i] = uiRest % m_uiGpuGrid[i];
        uiRest = uiRest / m_uiGpuGrid[i];
    }

    for (UINT i = 0; i < 4; ++i)
    {
        m_uiLocalLattice[i] = pGlobalLattice[i] / m_uiGpuGrid[i];
        m_uiGlobalLattice[i] = pGlobalLattice[i];
        m_uiGlobalOffset[i] = m_uiGridCoord[i] * m_uiLocalLattice[i];
    }

    return TRUE;
}

UINT CLGComm::NeighbourRank(UINT dir, INT sign) const
{
    if (dir > 3 || m_uiGpuGrid[dir] <= 1)
    {
        //Direction not split: the neighbour is ourselves (periodic within this rank).
        return m_uiRank;
    }

    UINT uiCoord[4];
    for (UINT i = 0; i < 4; ++i)
    {
        uiCoord[i] = m_uiGridCoord[i];
    }

    const UINT uiExtent = m_uiGpuGrid[dir];
    //Wrap around, so the process grid is a torus like the lattice itself.
    uiCoord[dir] = (uiCoord[dir] + uiExtent + static_cast<UINT>((sign >= 0) ? 1 : -1)) % uiExtent;

    UINT uiRank = 0;
    for (UINT i = 0; i < 4; ++i)
    {
        uiRank = uiRank * m_uiGpuGrid[i] + uiCoord[i];
    }
    return uiRank;
}

#pragma region Gather / scatter

#if _CLG_MULTI_GPU
namespace
{
    //Grid coordinate of an arbitrary rank, using the same x-slowest / t-fastest
    //convention as SetGpuGrid so rank 0 can place any rank's sub-lattice.
    void RankToGridCoord(UINT uiRank, const UINT* pGrid, UINT* pCoordOut)
    {
        UINT uiRest = uiRank;
        for (INT i = 3; i >= 0; --i)
        {
            pCoordOut[i] = uiRest % pGrid[i];
            uiRest = uiRest / pGrid[i];
        }
    }

    //Linear site index inside a lattice of dims [Lx,Ly,Lz,Lt]: x slowest, t fastest.
    inline UINT LocalSiteIndex(const UINT* pDims, UINT x, UINT y, UINT z, UINT w)
    {
        return ((x * pDims[1] + y) * pDims[2] + z) * pDims[3] + w;
    }
}
#endif

BYTE* CLGComm::GatherFieldToRoot(const BYTE* pLocalData, UINT uiBytesPerSite, UINT& uiOutSize) const
{
    const UINT uiLocalVolume = m_uiLocalLattice[0] * m_uiLocalLattice[1]
        * m_uiLocalLattice[2] * m_uiLocalLattice[3];
    const UINT uiLocalBytes = uiLocalVolume * uiBytesPerSite;

#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        //Even splits guarantee every rank owns the same local volume, so a flat
        //MPI_Gather of equal chunks is safe; rank 0 then de-interleaves by geometry.
        BYTE* pRecv = NULL;
        if (IsRoot())
        {
            pRecv = (BYTE*)malloc(static_cast<size_t>(uiLocalBytes) * m_uiSize);
        }

        MPI_Gather(
            const_cast<BYTE*>(pLocalData), static_cast<INT>(uiLocalBytes), MPI_BYTE,
            pRecv, static_cast<INT>(uiLocalBytes), MPI_BYTE,
            0, MPI_COMM_WORLD);

        if (!IsRoot())
        {
            uiOutSize = 0;
            return NULL;
        }

        const UINT uiGlobalVolume = m_uiGlobalLattice[0] * m_uiGlobalLattice[1]
            * m_uiGlobalLattice[2] * m_uiGlobalLattice[3];
        uiOutSize = uiGlobalVolume * uiBytesPerSite;
        BYTE* pGlobal = (BYTE*)malloc(uiOutSize);

        for (UINT r = 0; r < m_uiSize; ++r)
        {
            UINT uiCoord[4];
            RankToGridCoord(r, m_uiGpuGrid, uiCoord);
            UINT uiOffset[4];
            for (UINT d = 0; d < 4; ++d)
            {
                uiOffset[d] = uiCoord[d] * m_uiLocalLattice[d];
            }

            const BYTE* pChunk = pRecv + static_cast<size_t>(r) * uiLocalBytes;
            for (UINT lx = 0; lx < m_uiLocalLattice[0]; ++lx)
            for (UINT ly = 0; ly < m_uiLocalLattice[1]; ++ly)
            for (UINT lz = 0; lz < m_uiLocalLattice[2]; ++lz)
            for (UINT lw = 0; lw < m_uiLocalLattice[3]; ++lw)
            {
                const UINT uiLocalIdx = LocalSiteIndex(m_uiLocalLattice, lx, ly, lz, lw);
                const UINT uiGlobalIdx = LocalSiteIndex(m_uiGlobalLattice,
                    lx + uiOffset[0], ly + uiOffset[1], lz + uiOffset[2], lw + uiOffset[3]);
                memcpy(pGlobal + static_cast<size_t>(uiGlobalIdx) * uiBytesPerSite,
                    pChunk + static_cast<size_t>(uiLocalIdx) * uiBytesPerSite,
                    uiBytesPerSite);
            }
        }

        free(pRecv);
        return pGlobal;
    }
#endif

    //Single rank: local order already is global order.
    uiOutSize = uiLocalBytes;
    BYTE* pGlobal = (BYTE*)malloc(uiOutSize);
    memcpy(pGlobal, pLocalData, uiOutSize);
    return pGlobal;
}

void CLGComm::ScatterFieldFromRoot(const BYTE* pGlobalData, UINT uiBytesPerSite, BYTE* pLocalOut) const
{
    const UINT uiLocalVolume = m_uiLocalLattice[0] * m_uiLocalLattice[1]
        * m_uiLocalLattice[2] * m_uiLocalLattice[3];
    const UINT uiLocalBytes = uiLocalVolume * uiBytesPerSite;

#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        BYTE* pSend = NULL;
        if (IsRoot())
        {
            //Re-interleave the global buffer into per-rank contiguous chunks in
            //the exact order MPI_Scatter will hand them back out.
            pSend = (BYTE*)malloc(static_cast<size_t>(uiLocalBytes) * m_uiSize);
            for (UINT r = 0; r < m_uiSize; ++r)
            {
                UINT uiCoord[4];
                RankToGridCoord(r, m_uiGpuGrid, uiCoord);
                UINT uiOffset[4];
                for (UINT d = 0; d < 4; ++d)
                {
                    uiOffset[d] = uiCoord[d] * m_uiLocalLattice[d];
                }

                BYTE* pChunk = pSend + static_cast<size_t>(r) * uiLocalBytes;
                for (UINT lx = 0; lx < m_uiLocalLattice[0]; ++lx)
                for (UINT ly = 0; ly < m_uiLocalLattice[1]; ++ly)
                for (UINT lz = 0; lz < m_uiLocalLattice[2]; ++lz)
                for (UINT lw = 0; lw < m_uiLocalLattice[3]; ++lw)
                {
                    const UINT uiLocalIdx = LocalSiteIndex(m_uiLocalLattice, lx, ly, lz, lw);
                    const UINT uiGlobalIdx = LocalSiteIndex(m_uiGlobalLattice,
                        lx + uiOffset[0], ly + uiOffset[1], lz + uiOffset[2], lw + uiOffset[3]);
                    memcpy(pChunk + static_cast<size_t>(uiLocalIdx) * uiBytesPerSite,
                        pGlobalData + static_cast<size_t>(uiGlobalIdx) * uiBytesPerSite,
                        uiBytesPerSite);
                }
            }
        }

        MPI_Scatter(
            pSend, static_cast<INT>(uiLocalBytes), MPI_BYTE,
            pLocalOut, static_cast<INT>(uiLocalBytes), MPI_BYTE,
            0, MPI_COMM_WORLD);

        if (IsRoot())
        {
            free(pSend);
        }
        return;
    }
#endif

    //Single rank: global order already is local order.
    memcpy(pLocalOut, pGlobalData, uiLocalBytes);
}

#pragma endregion

#pragma region Collectives

void CLGComm::AllreduceSum(DOUBLE& value) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        DOUBLE result = 0.0;
        MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        value = result;
    }
#else
    (void)value;
#endif
}

void CLGComm::AllreduceSum(cuDoubleComplex& value) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        //Reduce as two contiguous doubles, avoiding a dependency on MPI complex types.
        DOUBLE buffer[2] = { value.x, value.y };
        DOUBLE result[2] = { 0.0, 0.0 };
        MPI_Allreduce(buffer, result, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        value.x = result[0];
        value.y = result[1];
    }
#else
    (void)value;
#endif
}

void CLGComm::AllreduceSum(cuDoubleComplex* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1 && NULL != pValues && uiCount > 0)
    {
        //Reduce as 2*uiCount contiguous doubles (no MPI complex dependency).
        MPI_Allreduce(MPI_IN_PLACE, pValues, static_cast<INT>(uiCount * 2), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#else
    (void)pValues;
    (void)uiCount;
#endif
}

#if !_CLG_DOUBLEFLOAT
void CLGComm::AllreduceSum(Real* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1 && NULL != pValues && uiCount > 0)
    {
        TArray<DOUBLE> staging;
        staging.AddItem(static_cast<DOUBLE>(pValues[0]));
        for (UINT i = 1; i < uiCount; ++i)
        {
            staging.AddItem(static_cast<DOUBLE>(pValues[i]));
        }
        MPI_Allreduce(MPI_IN_PLACE, staging.GetData(), static_cast<INT>(uiCount), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (UINT i = 0; i < uiCount; ++i)
        {
            pValues[i] = static_cast<Real>(staging[i]);
        }
    }
#else
    (void)pValues;
    (void)uiCount;
#endif
}
#endif //#if !_CLG_DOUBLEFLOAT

void CLGComm::AllreduceSum(DOUBLE* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1 && NULL != pValues && uiCount > 0)
    {
        MPI_Allreduce(MPI_IN_PLACE, pValues, static_cast<INT>(uiCount), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#else
    (void)pValues;
    (void)uiCount;
#endif
}

void CLGComm::BroadcastFromRoot(DOUBLE& value) const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
#else
    (void)value;
#endif
}

void CLGComm::Barrier() const
{
#if _CLG_MULTI_GPU
    if (m_uiSize > 1)
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
