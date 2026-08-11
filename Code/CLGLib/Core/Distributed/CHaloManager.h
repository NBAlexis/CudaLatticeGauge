//=============================================================================
// FILENAME : CHaloManager.h
//
// DESCRIPTION:
// Global halo manager for the multi-GPU build. See Docs/MultiGPU-Plan.md
// sections 3.1 / 3.2.
//
// Responsibility: make a field's halo valid to a requested width before a
// stencil kernel reads it, and mark it stale after a kernel writes it. The
// public _LAUNCH_KERNEL family drives this automatically (Improve-1 launch
// guard); call sites never talk to MPI directly.
//
// Phase 0: skeleton only -- validity tracking without exchange (there is no
// sub-lattice decomposition until Phase 1).
//
// REVISION:
//  [07/29/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CHALOMANAGER_H_
#define _CHALOMANAGER_H_

//Improve-1 (multi-GPU-improve1.md 3.1): buffer-extent halo identity; the
//manager keeps a debug registry of bound handles below. Included BEFORE the
//namespace block -- including it inside would nest a CLGLib::CLGLib namespace.
#include "Core/Distributed/CHaloBufferHandle.h"

__BEGIN_NAMESPACE

/**
* Improve-1 (multi-GPU-improve1.md 3.3): the per-launch guard produced by
* CHaloManager::BeginFromArguments. It records the DISTINCT managed buffer
* handles a launch's arguments hit (after the pre-launch Ensure of every
* HaloCapable hit) so that, once the backend launch has been submitted
* successfully, Commit() can conservatively invalidate every recorded buffer
* (3.3 rule 3: over-invalidating only costs an extra exchange).
*
* Host-management only; never part of the device kernel ABI.
*/
class CLGAPI CHaloLaunchGuard
{
public:

    CHaloLaunchGuard() : m_bCommitted(FALSE) {}

    /** Dedup: a handle hit by several arguments (or aliases) is processed once. */
    UBOOL AlreadyRecorded(const CHaloBufferHandle* pHandle) const;

    /** Record a managed hit; the handle must be bound (Info() fail-fasts otherwise). */
    void Record(const CHaloBufferHandle* pHandle);

    /**
    * Called by the public launch macro AFTER the backend launch's existing
    * immediate error check: NotifyWritten every recorded handle. Idempotent
    * (a second Commit is a no-op); the macro calls it exactly once. Adds no
    * device synchronize of its own.
    */
    void Commit();

    /** Number of distinct recorded handles (tests / diagnostics). */
    UINT RecordedCount() const { return static_cast<UINT>(m_lstRecorded.Num()); }

protected:

    TArray<const CHaloBufferHandle*> m_lstRecorded;
    UBOOL m_bCommitted;
};

/**
* Improve-1 (multi-GPU-improve1.md 3.4.1): host-side handle for a device
* POINTER ARRAY whose elements resolve to managed buffers (e.g. rational-
* approximation field arrays). The OWNER of the array binds the outer byte
* extent (the array allocation itself) together with the FIXED set of member
* handles its elements point to. When a launch argument hits the outer
* extent, the guard expands the set: every HaloCapable member is Ensured and
* every member is recorded for post-launch invalidation. After members
* reallocate/reorder the owner must Unbind/re-Bind (3.4.1: generation and
* rebind on membership change).
*
* Host-management only; never part of the device kernel ABI.
*/
class CLGAPI CHaloBufferSetHandle
{
public:

    CHaloBufferSetHandle()
        : m_pOuter(NULL)
        , m_uiOuterBytes(0)
        , m_ullGeneration(0)
        , m_bBound(FALSE)
    {
    }

    /**
    * Bind the outer extent [pOuter, pOuter + uiOuterBytes) to a fixed member
    * set. Every member must be an already-bound CHaloBufferHandle. Fails
    * (logs CRUCIAL) when the outer extent overlaps a registered handle or
    * another set's outer extent, or when a member is unbound.
    */
    UBOOL Bind(const BYTE* pOuter, size_t uiOuterBytes, const TArray<const CHaloBufferHandle*>& members);
    void Unbind();
    UBOOL IsBound() const { return m_bBound; }

    const BYTE* OuterBegin() const { return m_pOuter; }
    size_t OuterBytes() const { return m_uiOuterBytes; }
    UINT MemberCount() const { return static_cast<UINT>(m_Members.Num()); }
    const CHaloBufferHandle* Member(UINT uiIndex) const { return m_Members[static_cast<INT>(uiIndex)]; }
    ULONGLONG Generation() const { return m_ullGeneration; }

protected:

    TArray<const CHaloBufferHandle*> m_Members;
    const BYTE* m_pOuter;
    size_t m_uiOuterBytes;
    ULONGLONG m_ullGeneration;
    UBOOL m_bBound;
};

/**
* Tracks, PER BUFFER HANDLE, how wide a halo is currently valid
* (multi-GPU-improve1.md 3.1/3.6): the version state lives in each
* CHaloBufferHandle (interior version / halo version / valid width), never in
* a per-field-id table, so two buffers sharing a field id hold independent
* halos.
*
* Correctness rule (Docs/MultiGPU-Plan.md section 3.2): a kernel that writes a
* field invalidates that field's halo. Over-marking dirty only costs an extra
* exchange; under-marking is a correctness bug. Prefer conservative marking.
*/
class CLGAPI CHaloManager
{
public:

    CHaloManager();
    ~CHaloManager();

    /**
    * Ensure the buffer's halo is valid to at least uiWidth sites. Cheap no-op
    * when the handle already records a valid halo of that width for the
    * CURRENT interior version, so putting this in front of every stencil
    * launch is fine. Validates uiWidth <= configured halo width; a failed or
    * unsupported refill never updates the version state.
    */
    void Ensure(CHaloBufferHandle& buffer, UINT uiWidth);

    /** Invalidate every tracked halo, e.g. after loading a configuration. */
    void MarkAllDirty();

    /**
    * Handle-targeted halo refill (multi-GPU-improve1.md 3.6): refill the halo
    * tail of the buffer described by this exact handle -- this is the primary
    * entry point; pointer/stride/counts come ONLY from the bound
    * SHaloBufferInfo. Returns FALSE (and never touches version state) when the
    * handle is unbound, not halo-capable, carries a stale layout generation,
    * or has an unusable descriptor; returns TRUE when the halo is in sync
    * (including the vacuous nothing-split case). Does NOT itself update the
    * version state -- Ensure does that on success.
    */
    UBOOL RefillHalo(CHaloBufferHandle& buffer);

    #pragma region Improve-1 handle registry

    /**
    * Debug registry of bound CHaloBufferHandle extents (multi-GPU-improve1.md
    * 3.1): two simultaneously bound handles may NOT declare overlapping byte
    * extents -- identical/overlapping views must share one handle, so a second
    * Bind over the same bytes is always a bug. Returns FALSE (and reports) on
    * overlap; the caller's Bind then fails without changing state. Called from
    * CHaloBufferHandle::Bind/Unbind only.
    */
    UBOOL OnHandleBound(const class CHaloBufferHandle* pHandle, const SHaloBufferInfo& info);
    void OnHandleUnbound(const class CHaloBufferHandle* pHandle);

    /** Number of currently registered handles (tests / diagnostics). */
    UINT RegisteredHandleCount() const { return static_cast<UINT>(m_RegisteredHandles.Num()); }

    #pragma endregion

    #pragma region Improve-1 launch guard

    /**
    * Registry range lookup (multi-GPU-improve1.md 3.3 rule 4): the bound
    * handle whose byte extent CONTAINS pArg -- a pointer into the interior of
    * an extent belongs to that extent. NULL when pArg is not a managed
    * buffer, so unregistered pointers stay ordinary kernel arguments.
    */
    CHaloBufferHandle* LookupByPointer(const void* pArg);

    /**
    * Set range lookup (3.4.1): the bound CHaloBufferSetHandle whose OUTER byte
    * extent contains pArg (a pointer array allocation), NULL when no set
    * matches. Consulted only after LookupByPointer misses.
    */
    const CHaloBufferSetHandle* LookupSetByPointer(const void* pArg) const;

    /** Configured full halo width (reads _HC_HaloWidth in the .cpp, which is
    * included after CudaHelper.h -- this header is not). */
    UINT ConfiguredHaloWidth() const;

    /** Set registration; CHaloBufferSetHandle::Bind/Unbind only. */
    UBOOL OnSetBound(const CHaloBufferSetHandle* pSet);
    void OnSetUnbound(const CHaloBufferSetHandle* pSet);

    /** Number of currently registered set handles (tests / diagnostics). */
    UINT RegisteredSetCount() const { return static_cast<UINT>(m_RegisteredSets.Num()); }

    /**
    * Improve-1 (3.3): inspect every launch argument exactly once, BEFORE the
    * backend launch. A pointer argument hitting a registered HaloCapable
    * handle is Ensured to the FULL configured halo width (Improve 1 never
    * tries to shrink the exchange); every managed hit (HaloCapable or
    * LocalOnly) is recorded in the returned guard exactly once, so the
    * successful launch can conservatively invalidate it. Non-pointer
    * arguments are ordinary kernel parameters, ignored at compile time.
    */
    template<typename... TArgs>
    CHaloLaunchGuard BeginFromArguments(const TArgs&... args)
    {
        if (m_uiRefillInFlight > 0)
        {
            //A guarded launch from inside a refill means an infrastructure
            //kernel used the public macro instead of _CLG_LAUNCH_KERNEL_RAW.
            appCrucial(_T("CHaloManager::BeginFromArguments: entered while a halo refill is in flight (guarded launch recursion). Use the private raw macro for infrastructure kernels.\n"));
        }
        CHaloLaunchGuard guard;
        using TExpand = INT[];
        (void)TExpand{ 0, (CollectLaunchArgument(args, guard), 0)... };
        return guard;
    }

    /**
    * Improve-1 (3.4.4): rank-local variant for _LAUNCH_KERNEL_RANK_LOCAL.
    * Same single-evaluation traversal, but a HaloCapable hit is a caller bug
    * (fail-fast: a rank-asymmetric launch may never touch a distributed
    * halo). LocalOnly hits are recorded for post-launch invalidation only --
    * NO Ensure and NO MPI ever happens here.
    */
    template<typename... TArgs>
    CHaloLaunchGuard BeginFromArgumentsRankLocal(const TArgs&... args)
    {
        CHaloLaunchGuard guard;
        using TExpand = INT[];
        (void)TExpand{ 0, (CollectLaunchArgumentRankLocal(args, guard), 0)... };
        return guard;
    }

    #pragma endregion

protected:

    /**
    * Exchange engine behind RefillHalo(CHaloBufferHandle&) -- its only caller.
    * Gathers this rank's boundary sites, exchanges the FACE/EDGE/CORNER/HYPER
    * halo blocks with the process-torus neighbours (per-block-kind disjoint
    * MPI tag spaces derived from byFieldId), and writes the results into the
    * halo tail of pFieldBytes. byFieldId only tags the MPI messages.
    */
    void RefillHaloBuffer(BYTE* pFieldBytes, UINT uiBytesPerSite, BYTE byFieldId);

    /** Per-handle action shared by direct hits and set expansion: Ensure a
    * HaloCapable handle to the FULL configured width, then record it once. */
    void CollectHandle(CHaloBufferHandle* pHandle, CHaloLaunchGuard& guard)
    {
        if (NULL == pHandle || guard.AlreadyRecorded(pHandle))
        {
            return;
        }
        if (pHandle->Info().m_bHaloCapable)
        {
            Ensure(*pHandle, ConfiguredHaloWidth());
        }
        guard.Record(pHandle);
    }

    /** Set expansion (3.4.1): every member of a hit set is processed. */
    void CollectSet(const CHaloBufferSetHandle* pSet, CHaloLaunchGuard& guard)
    {
        for (UINT i = 0; i < pSet->MemberCount(); ++i)
        {
            //const_cast: set membership is stored const for read-back; the
            //guard drives Ensure/Record on the member just like a direct hit.
            CollectHandle(const_cast<CHaloBufferHandle*>(pSet->Member(i)), guard);
        }
    }

    /** Per-argument dispatch: pointers go through the registry (direct
    * handle hit first, then set outer extent), everything else is an
    * ordinary kernel argument (compile-time no-op). */
    template<typename T>
    void CollectLaunchArgument(T* pArg, CHaloLaunchGuard& guard)
    {
        if (NULL == pArg)
        {
            return;
        }
        CHaloBufferHandle* pHandle = LookupByPointer(pArg);
        if (NULL != pHandle)
        {
            CollectHandle(pHandle, guard);
            return;
        }
        const CHaloBufferSetHandle* pSet = LookupSetByPointer(pArg);
        if (NULL != pSet)
        {
            CollectSet(pSet, guard);
        }
    }

    template<typename T>
    void CollectLaunchArgument(const T&, CHaloLaunchGuard&)
    {
        //Non-pointer kernel argument: nothing to do.
    }

    /** Function-pointer kernel argument (device callbacks): never a managed
    * data extent; skipped so the traversal never casts it to const void*. */
    template<typename TRet, typename... TParams>
    void CollectLaunchArgument(TRet(*)(TParams...), CHaloLaunchGuard&) {}

    /** Rank-local per-handle action: HaloCapable hit fails fast (caller bug,
    * NOT recorded), LocalOnly is recorded for post-launch invalidation. */
    void CollectHandleRankLocal(CHaloBufferHandle* pHandle, CHaloLaunchGuard& guard)
    {
        if (NULL == pHandle || guard.AlreadyRecorded(pHandle))
        {
            return;
        }
        if (pHandle->Info().m_bHaloCapable)
        {
            appCrucial(_T("CHaloManager: _LAUNCH_KERNEL_RANK_LOCAL argument hits a HaloCapable buffer; rank-asymmetric launches must use collective/scatter semantics instead.\n"));
            return;
        }
        guard.Record(pHandle);
    }

    /** Rank-local variant: HaloCapable hit fails fast, never Ensures. */
    template<typename T>
    void CollectLaunchArgumentRankLocal(T* pArg, CHaloLaunchGuard& guard)
    {
        if (NULL == pArg)
        {
            return;
        }
        CHaloBufferHandle* pHandle = LookupByPointer(pArg);
        if (NULL != pHandle)
        {
            CollectHandleRankLocal(pHandle, guard);
            return;
        }
        const CHaloBufferSetHandle* pSet = LookupSetByPointer(pArg);
        if (NULL != pSet)
        {
            for (UINT i = 0; i < pSet->MemberCount(); ++i)
            {
                CollectHandleRankLocal(const_cast<CHaloBufferHandle*>(pSet->Member(i)), guard);
            }
        }
    }

    template<typename T>
    void CollectLaunchArgumentRankLocal(const T&, CHaloLaunchGuard&)
    {
        //Non-pointer kernel argument: nothing to do.
    }

    /** Function-pointer kernel argument: never a managed extent, skipped. */
    template<typename TRet, typename... TParams>
    void CollectLaunchArgumentRankLocal(TRet(*)(TParams...), CHaloLaunchGuard&) {}

    //Improve-1: bound handle registry (pointers only; extents are read back
    //through each handle's Info()).
    TArray<const class CHaloBufferHandle*> m_RegisteredHandles;

    //Improve-1 (3.4.1): bound set-handle registry (device pointer arrays).
    TArray<const class CHaloBufferSetHandle*> m_RegisteredSets;

    //Improve-1: anti-recursion diagnostic; non-zero while a handle refill is
    //in flight (BeginFromArguments must never be re-entered from there).
    UINT m_uiRefillInFlight;

    //Temporary host staging buffers for MPI exchange (MS-MPI is not CUDA-aware).
    //Per (side) send/recv buffers so a non-blocking Isend/Irecv/Waitall exchange of
    //both faces of a split direction can be in flight simultaneously without
    //deadlock (blocking Sendrecv deadlocks on a size-2 torus where both neighbours
    //are the same rank). Sized to the largest single-face block; reused across
    //fields. NULL / 0-sized when nothing is split.
    BYTE* m_pCommBufferSend[2];
    BYTE* m_pCommBufferRecv[2];
    UINT m_uiCommBufferBytes;
};

//Defined as an inline accessor over GCLGManager in Core/CLGLibManager.h,
//matching the existing GetBuffer() / appGetCudaHelper() convention.
inline class CHaloManager* appGetHaloManager();

__END_NAMESPACE

#endif //#ifndef _CHALOMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================
