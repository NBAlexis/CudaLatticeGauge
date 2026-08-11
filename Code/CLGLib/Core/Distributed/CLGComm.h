//=============================================================================
// FILENAME : CLGComm.h
//
// DESCRIPTION:
// Backend-neutral communication layer for the multi-GPU (one MPI process per
// GPU) build. See Docs/MultiGPU-Plan.md section 3.
//
// When _CLG_MULTI_GPU is 0, every function here degenerates to a single-rank
// no-op so the single-GPU code path is bit-for-bit unchanged. <mpi.h> is only
// included in the .cpp under the _CLG_MULTI_GPU guard, so a machine without
// MPI installed can still build the single-GPU configurations.
//
// REVISION:
//  [07/29/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CLGCOMM_H_
#define _CLGCOMM_H_

__BEGIN_NAMESPACE

/**
* Process topology and collectives for the multi-GPU build.
*
* The process grid is always treated as 4 dimensional (Px, Py, Pz, Pt), even
* when only one direction is actually split. Splitting a single direction is a
* runtime configuration (GpuGrid), never a hard-coded assumption -- see
* Docs/MultiGPU-Plan.md section 8.1-A6.
*/
class CLGAPI CLGComm
{
public:

    CLGComm();
    ~CLGComm();

    /**
    * Initialise MPI and work out this rank's position in the process grid.
    * Safe (and a no-op beyond setting rank 0 / size 1) in single-GPU builds.
    */
    void Initial(INT* argc, TCHAR*** argv);
    void Finalize();

    UINT Rank() const { return m_uiRank; }
    UINT Size() const { return m_uiSize; }
    UBOOL IsRoot() const { return 0 == m_uiRank; }

    /**
    * Process grid shape and this rank's coordinate within it.
    * In single-GPU builds the grid is [1,1,1,1] and the coordinate is [0,0,0,0].
    */
    const UINT* GpuGrid() const { return m_uiGpuGrid; }
    const UINT* GridCoord() const { return m_uiGridCoord; }

    /**
    * Set the process grid, validating it against the global lattice.
    * Returns FALSE (and reports via appCrucial) when a hard constraint is
    * violated -- see Docs/MultiGPU-Plan.md section 9.2:
    *   - product of the grid must equal the MPI world size
    *   - every lattice length must divide evenly by its grid factor
    *   - every split direction's local length must be >= the halo width
    */
    UBOOL SetGpuGrid(const UINT* pGrid, const UINT* pGlobalLattice, UINT uiHaloWidth);

    /**
    * Neighbour rank in direction dir (0..3) with sign +1 / -1.
    * Returns this rank when that direction is not split (grid factor 1).
    */
    UINT NeighbourRank(UINT dir, INT sign) const;

    /** Global offset of this rank's sub-lattice, for local -> global coordinates. */
    const UINT* GlobalOffset() const { return m_uiGlobalOffset; }

    /** Global lattice dimensions (set by SetGpuGrid). */
    const UINT* GlobalLattice() const { return m_uiGlobalLattice; }

    /** This rank's local sub-lattice dimensions (global / grid). */
    const UINT* LocalLattice() const { return m_uiLocalLattice; }

    #pragma region Gather / scatter (whole-field, for IO and 1-vs-N validation)

    /**
    * Gather a decomposed field to rank 0 in GLOBAL site order.
    *
    * Each rank passes its local field bytes exactly as produced by
    * CField::CopyDataOut (local-site order, uiBytesPerSite bytes per site, sites
    * laid out with x slowest / t fastest). Returns a freshly malloc'd buffer on
    * rank 0 holding the whole global lattice in global-site order; returns NULL
    * on every non-root rank. Caller frees. Used to diff -n1 vs -nN (see
    * Docs/MultiGPU-Plan.md section 8.5); no halo involved.
    *
    * Single-rank builds return a plain copy of the input.
    */
    BYTE* GatherFieldToRoot(const BYTE* pLocalData, UINT uiBytesPerSite, UINT& uiOutSize) const;

    /**
    * Inverse of GatherFieldToRoot: rank 0 holds a global-site-ordered buffer,
    * every rank receives its own sub-lattice in local-site order. pGlobalData is
    * read on rank 0 only; pLocalOut must be sized localVolume * uiBytesPerSite on
    * every rank.
    */
    void ScatterFieldFromRoot(const BYTE* pGlobalData, UINT uiBytesPerSite, BYTE* pLocalOut) const;

    #pragma endregion

    #pragma region Collectives

    /** In-place sum across all ranks. No-op when running on a single rank. */
    void AllreduceSum(DOUBLE& value) const;
    void AllreduceSum(cuDoubleComplex& value) const;
    void AllreduceSum(DOUBLE* pValues, UINT uiCount) const;
    void AllreduceSum(cuDoubleComplex* pValues, UINT uiCount) const;
    //P4-3.8: Real (float on float builds) array variant; staged through DOUBLE
    //for a single portable MPI reduce. Declared only on float builds: on double
    //builds Real == DOUBLE, so this signature would collide with the DOUBLE*
    //overload above (M4 double-build fix).
#if !_CLG_DOUBLEFLOAT
    void AllreduceSum(Real* pValues, UINT uiCount) const;
#endif

    /**
    * Barrier across all ranks. No-op on a single rank. Needed around root-only
    * file writes that other ranks subsequently read (e.g. compressed save ->
    * reload round-trips): without it a non-root rank can open the file while
    * the root is mid-write (fopen "wb" truncates), fail the size check, and
    * _FAIL_EXIT -> exit() -> MPI_Finalize then blocks forever against the
    * root's pending collective (I10 gate hang).
    */
    void Barrier() const;

    /**
    * Broadcast a value from rank 0 to every rank, in place. No-op on a single
    * rank. Used to keep host-side stochastic decisions (e.g. the Metropolis
    * accept/reject draw) identical across ranks -- each rank has its own RNG
    * stream, so without this they would diverge even when the global energy is
    * the same on all ranks.
    */
    void BroadcastFromRoot(DOUBLE& value) const;

    #pragma endregion

protected:

    UINT m_uiRank;
    UINT m_uiSize;
    UINT m_uiGpuGrid[4];
    UINT m_uiGridCoord[4];
    UINT m_uiGlobalOffset[4];
    UINT m_uiGlobalLattice[4];
    UINT m_uiLocalLattice[4];
    UBOOL m_bInitialed;
};

//Defined as an inline accessor over GCLGManager in Core/CLGLibManager.h,
//matching the existing GetBuffer() / appGetCudaHelper() convention.
inline class CLGComm* appGetComm();

/**
* I9: unified global-sum entry point for host-side scalars outside CMeasure
* (e.g. action energies). Thin wrapper over appGetComm()->AllreduceSum so call
* sites no longer spell the collective (or the _CLG_MULTI_GPU guard) directly;
* no-op when _CLG_MULTI_GPU is 0 or the comm is not initialised (single rank),
* so single-GPU behaviour is bit-for-bit unchanged.
*/
inline void appGlobalSum(DOUBLE& value)
{
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(value);
    }
#endif
}

inline void appGlobalSum(cuDoubleComplex& value)
{
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(value);
    }
#endif
}

/**
* I9: cuDoubleComplex array variant (e.g. PolyakovXY local-slice density
* arrays, whose element semantics -- local slice vs global slice -- stay with
* the call site and its GpuGrid guards; this helper only changes HOW the sum
* is done). Always double precision, so it is valid on both float and double
* builds. No-op when _CLG_MULTI_GPU is 0 or the comm is not initialised.
*/
inline void appGlobalSum(cuDoubleComplex* pValues, UINT uiCount)
{
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && NULL != pValues && uiCount > 0)
    {
        appGetComm()->AllreduceSum(pValues, uiCount);
    }
#endif
}

/** Only rank 0 should emit ordinary progress logs; every rank reports errors. */
#define _CLG_IS_LOG_RANK (NULL == appGetComm() || appGetComm()->IsRoot())

__END_NAMESPACE

#endif //#ifndef _CLGCOMM_H_

//=============================================================================
// END OF FILE
//=============================================================================
