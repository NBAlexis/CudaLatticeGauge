//=============================================================================
// FILENAME : CHaloManager.cpp
//
// DESCRIPTION:
// See CHaloManager.h.
//
// Phase 0 scope: validity bookkeeping only. The actual pack / MPI exchange /
// unpack lands in Phase 2, once Phase 1 has given each rank a sub-lattice and
// pointed out-of-range neighbour SIndex entries at halo storage. Deliberately
// synchronous (no streams) -- see Docs/MultiGPU-Plan.md section 3.
//
// REVISION:
//  [07/29/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"

#if _CLG_MULTI_GPU
#include <mpi.h>
#endif

__BEGIN_NAMESPACE

CHaloManager::CHaloManager()
    : m_uiCommBufferBytes(0)
    , m_uiRefillInFlight(0)
{
    m_pCommBufferSend[0] = NULL;
    m_pCommBufferSend[1] = NULL;
    m_pCommBufferRecv[0] = NULL;
    m_pCommBufferRecv[1] = NULL;
    MarkAllDirty();
}

CHaloManager::~CHaloManager()
{
    //Host staging buffers (MS-MPI is not CUDA-aware): plain free, not cudaFree.
    for (UINT s = 0; s < 2; ++s)
    {
        if (NULL != m_pCommBufferSend[s]) { free(m_pCommBufferSend[s]); m_pCommBufferSend[s] = NULL; }
        if (NULL != m_pCommBufferRecv[s]) { free(m_pCommBufferRecv[s]); m_pCommBufferRecv[s] = NULL; }
    }
}

void CHaloManager::Ensure(CHaloBufferHandle& buffer, UINT uiWidth)
{
#if _CLG_MULTI_GPU
    if (!buffer.IsBound())
    {
        appCrucial(_T("CHaloManager::Ensure: unbound handle.\n"));
        return;
    }
    const SHaloBufferInfo& info = buffer.Info();
    if (!info.m_bHaloCapable)
    {
        //A LocalOnly buffer must never need a halo; Ensuring one is a caller bug.
        if (uiWidth > 0)
        {
            appCrucial(_T("CHaloManager::Ensure: buffer is LocalOnly (not halo-capable) but width %d requested.\n"), uiWidth);
        }
        return;
    }
    if (info.m_ullLayoutGeneration != appGetLayoutGeneration())
    {
        appCrucial(_T("CHaloManager::Ensure: stale layout generation %llu (current %llu); the owner must re-Bind after a layout rebuild.\n"),
            info.m_ullLayoutGeneration, appGetLayoutGeneration());
        return;
    }
    if (uiWidth > static_cast<UINT>(_HC_HaloWidth))
    {
        appCrucial(_T("CHaloManager::Ensure: required width %d exceeds configured halo width %d.\n"),
            uiWidth, static_cast<UINT>(_HC_HaloWidth));
        return;
    }

    //Already valid to at least this width for the CURRENT interior: nothing to
    //do. Two Ensures without an intermediate write never change any version.
    if (buffer.ValidWidth() >= uiWidth && buffer.HaloVersion() == buffer.InteriorVersion())
    {
        return;
    }

    //Improve-1 always refills the FULL configured halo (3.6); the recorded
    //valid width is the actually filled width. A failed refill leaves the
    //version state untouched.
    if (RefillHalo(buffer))
    {
        buffer.NotifyHaloSynchronized(static_cast<UINT>(_HC_HaloWidth));
    }
#else
    (void)buffer;
    (void)uiWidth;
#endif
}

void CHaloManager::MarkAllDirty()
{
    for (INT i = 0; i < m_RegisteredHandles.Num(); ++i)
    {
        // const_cast: the registry stores const pointers for read-back
        // diagnostics; NotifyWritten is the one state transition MarkAllDirty
        // is allowed to drive on every tracked extent.
        const_cast<CHaloBufferHandle*>(m_RegisteredHandles[i])->NotifyWritten();
    }
}

#pragma region Improve-1 handle registry

UBOOL CHaloManager::OnHandleBound(const CHaloBufferHandle* pHandle, const SHaloBufferInfo& info)
{
    //Reject a second bound handle over any overlapping byte extent: identical
    //or overlapping views must SHARE one handle (multi-GPU-improve1.md 3.1).
    //Note a zero-capacity extent [p, p) never overlaps.
    const BYTE* pNewBegin = info.m_pDeviceData;
    const BYTE* pNewEnd = info.m_pDeviceData + info.m_uiCapacityBytes;
    for (INT i = 0; i < m_RegisteredHandles.Num(); ++i)
    {
        const CHaloBufferHandle* pOther = m_RegisteredHandles[i];
        if (pOther == pHandle)
        {
            continue;
        }
        const SHaloBufferInfo& other = pOther->Info();
        const BYTE* pOtherBegin = other.m_pDeviceData;
        const BYTE* pOtherEnd = other.m_pDeviceData + other.m_uiCapacityBytes;
        if (pNewBegin < pOtherEnd && pOtherBegin < pNewEnd)
        {
            appCrucial(_T("CHaloManager::OnHandleBound: extent [%p, +%llu) overlaps a bound handle [%p, +%llu); overlapping views must share one handle.\n"),
                pNewBegin, static_cast<ULONGLONG>(info.m_uiCapacityBytes),
                pOtherBegin, static_cast<ULONGLONG>(other.m_uiCapacityBytes));
            return FALSE;
        }
    }
#if _CLG_MULTI_GPU
    //Improve-1 (multi-GPU-improve1.md 3.8): a halo-capable extent's site
    //counts must be exactly consistent with the CURRENT layout. The owner's
    //per-site element multiplier (sites vs links) is shared by both counts,
    //so cross-multiplication checks consistency without knowing it. The
    //generation must be current too, or _HC_Volume/_HC_HaloSiteCount below
    //would describe a layout the owner no longer has. (LocalOnly extents may
    //be any scratch shape and are not checked.)
    if (info.m_bHaloCapable)
    {
        if (info.m_ullLayoutGeneration != appGetLayoutGeneration())
        {
            appCrucial(_T("CHaloManager::OnHandleBound: halo-capable extent bound with stale layout generation %llu (current %llu); re-Bind after the latest layout bake.\n"),
                info.m_ullLayoutGeneration, appGetLayoutGeneration());
            return FALSE;
        }
        const ULONGLONG ullLayoutVolume = static_cast<ULONGLONG>(_HC_Volume);
        const ULONGLONG ullLayoutHalo = static_cast<ULONGLONG>(_HC_HaloSiteCount());
        if (static_cast<ULONGLONG>(info.m_uiHaloSiteCount) * ullLayoutVolume
            != static_cast<ULONGLONG>(info.m_uiLocalSiteCount) * ullLayoutHalo)
        {
            appCrucial(_T("CHaloManager::OnHandleBound: counts (local %d, halo %d) are inconsistent with the current layout (volume %d, halo %d).\n"),
                static_cast<INT>(info.m_uiLocalSiteCount), static_cast<INT>(info.m_uiHaloSiteCount),
                static_cast<INT>(_HC_Volume), static_cast<INT>(_HC_HaloSiteCount()));
            return FALSE;
        }
    }
#endif

    m_RegisteredHandles.AddItem(pHandle);
    return TRUE;
}

void CHaloManager::OnHandleUnbound(const CHaloBufferHandle* pHandle)
{
    for (INT i = 0; i < m_RegisteredHandles.Num(); ++i)
    {
        if (m_RegisteredHandles[i] == pHandle)
        {
            m_RegisteredHandles.RemoveAt(i);
            return;
        }
    }
}

#pragma endregion

#pragma region Improve-1 launch guard

UBOOL CHaloLaunchGuard::AlreadyRecorded(const CHaloBufferHandle* pHandle) const
{
    for (INT i = 0; i < m_lstRecorded.Num(); ++i)
    {
        if (m_lstRecorded[i] == pHandle)
        {
            return TRUE;
        }
    }
    return FALSE;
}

void CHaloLaunchGuard::Record(const CHaloBufferHandle* pHandle)
{
    if (NULL != pHandle && !AlreadyRecorded(pHandle))
    {
        m_lstRecorded.AddItem(pHandle);
    }
}

void CHaloLaunchGuard::Commit()
{
    if (m_bCommitted)
    {
        return;
    }
    m_bCommitted = TRUE;
    //Improve-1 (3.3 rule 3): the launch has been submitted; conservatively
    //invalidate every managed buffer its arguments touched, even the ones it
    //only read. Over-invalidating costs one extra exchange; under-marking is
    //a correctness bug.
    for (INT i = 0; i < m_lstRecorded.Num(); ++i)
    {
        //const_cast: the guard records const pointers for read-back; the
        //post-launch NotifyWritten is the one state transition it must drive.
        const_cast<CHaloBufferHandle*>(m_lstRecorded[i])->NotifyWritten();
    }
}

CHaloBufferHandle* CHaloManager::LookupByPointer(const void* pArg)
{
    if (NULL == pArg)
    {
        return NULL;
    }
    const BYTE* p = static_cast<const BYTE*>(pArg);
    for (INT i = 0; i < m_RegisteredHandles.Num(); ++i)
    {
        const SHaloBufferInfo& info = m_RegisteredHandles[i]->Info();
        if (p >= info.m_pDeviceData && p < info.m_pDeviceData + info.m_uiCapacityBytes)
        {
            //const_cast: registry stores const pointers for read-back; the
            //guard needs the mutable handle for Ensure/Record.
            return const_cast<CHaloBufferHandle*>(m_RegisteredHandles[i]);
        }
    }
    return NULL;
}

UINT CHaloManager::ConfiguredHaloWidth() const
{
    return static_cast<UINT>(_HC_HaloWidth);
}

const CHaloBufferSetHandle* CHaloManager::LookupSetByPointer(const void* pArg) const
{
    if (NULL == pArg)
    {
        return NULL;
    }
    const BYTE* p = static_cast<const BYTE*>(pArg);
    for (INT i = 0; i < m_RegisteredSets.Num(); ++i)
    {
        const CHaloBufferSetHandle* pSet = m_RegisteredSets[i];
        if (p >= pSet->OuterBegin() && p < pSet->OuterBegin() + pSet->OuterBytes())
        {
            return pSet;
        }
    }
    return NULL;
}

UBOOL CHaloManager::OnSetBound(const CHaloBufferSetHandle* pSet)
{
    //The outer extent (a device pointer array) must not overlap any managed
    //data extent or another set's outer extent.
    const BYTE* pNewBegin = pSet->OuterBegin();
    const BYTE* pNewEnd = pSet->OuterBegin() + pSet->OuterBytes();
    for (INT i = 0; i < m_RegisteredHandles.Num(); ++i)
    {
        const SHaloBufferInfo& other = m_RegisteredHandles[i]->Info();
        const BYTE* pOtherBegin = other.m_pDeviceData;
        const BYTE* pOtherEnd = other.m_pDeviceData + other.m_uiCapacityBytes;
        if (pNewBegin < pOtherEnd && pOtherBegin < pNewEnd)
        {
            appCrucial(_T("CHaloManager::OnSetBound: outer extent [%p, +%llu) overlaps a bound handle [%p, +%llu).\n"),
                pNewBegin, static_cast<ULONGLONG>(pSet->OuterBytes()),
                pOtherBegin, static_cast<ULONGLONG>(other.m_uiCapacityBytes));
            return FALSE;
        }
    }
    for (INT i = 0; i < m_RegisteredSets.Num(); ++i)
    {
        const CHaloBufferSetHandle* pOther = m_RegisteredSets[i];
        if (pOther == pSet)
        {
            continue;
        }
        const BYTE* pOtherBegin = pOther->OuterBegin();
        const BYTE* pOtherEnd = pOther->OuterBegin() + pOther->OuterBytes();
        if (pNewBegin < pOtherEnd && pOtherBegin < pNewEnd)
        {
            appCrucial(_T("CHaloManager::OnSetBound: outer extent [%p, +%llu) overlaps another set [%p, +%llu).\n"),
                pNewBegin, static_cast<ULONGLONG>(pSet->OuterBytes()),
                pOtherBegin, static_cast<ULONGLONG>(pOther->OuterBytes()));
            return FALSE;
        }
    }
    m_RegisteredSets.AddItem(pSet);
    return TRUE;
}

void CHaloManager::OnSetUnbound(const CHaloBufferSetHandle* pSet)
{
    for (INT i = 0; i < m_RegisteredSets.Num(); ++i)
    {
        if (m_RegisteredSets[i] == pSet)
        {
            m_RegisteredSets.RemoveAt(i);
            return;
        }
    }
}

UBOOL CHaloBufferSetHandle::Bind(const BYTE* pOuter, size_t uiOuterBytes, const TArray<const CHaloBufferHandle*>& members)
{
    if (m_bBound)
    {
        appCrucial(_T("CHaloBufferSetHandle::Bind: already bound (Unbind first).\n"));
        return FALSE;
    }
    if (NULL == pOuter || 0 == uiOuterBytes || members.Num() <= 0)
    {
        appCrucial(_T("CHaloBufferSetHandle::Bind: empty outer extent or empty member set.\n"));
        return FALSE;
    }
    for (INT i = 0; i < members.Num(); ++i)
    {
        if (NULL == members[i] || !members[i]->IsBound())
        {
            appCrucial(_T("CHaloBufferSetHandle::Bind: member %d is NULL or unbound.\n"), i);
            return FALSE;
        }
    }
    m_Members = members;
    m_pOuter = pOuter;
    m_uiOuterBytes = uiOuterBytes;
    m_ullGeneration = appGetLayoutGeneration();
    if (NULL == appGetHaloManager() || !appGetHaloManager()->OnSetBound(this))
    {
        m_Members.RemoveAll();
        m_pOuter = NULL;
        m_uiOuterBytes = 0;
        m_ullGeneration = 0;
        return FALSE;
    }
    m_bBound = TRUE;
    return TRUE;
}

void CHaloBufferSetHandle::Unbind()
{
    if (!m_bBound)
    {
        return;
    }
    if (NULL != appGetHaloManager())
    {
        appGetHaloManager()->OnSetUnbound(this);
    }
    m_Members.RemoveAll();
    m_pOuter = NULL;
    m_uiOuterBytes = 0;
    m_ullGeneration = 0;
    m_bBound = FALSE;
}

#pragma endregion

UBOOL CHaloManager::RefillHalo(CHaloBufferHandle& buffer)
{
#if _CLG_MULTI_GPU
    //Improve-1 (multi-GPU-improve1.md 3.6, I4): the ONLY source of
    //pointer/stride/counts is this handle's bound descriptor. IsBound() is
    //checked before Info() because Info() fail-fasts on an unbound handle.
    if (!buffer.IsBound())
    {
        return FALSE;
    }
    const SHaloBufferInfo& info = buffer.Info();
    if (!info.m_bHaloCapable)
    {
        return FALSE;
    }
    if (info.m_ullLayoutGeneration != appGetLayoutGeneration())
    {
        appCrucial(_T("CHaloManager::RefillHalo: handle carries stale layout generation %llu (current %llu); the owner must re-Bind after a layout rebuild. Refill refused.\n"),
            info.m_ullLayoutGeneration, appGetLayoutGeneration());
        return FALSE;
    }
    if (NULL == info.m_pDeviceData || 0 == info.m_uiBytesPerSite)
    {
        return FALSE;
    }

    const CLGComm* pComm = appGetComm();
    if (NULL == pComm || pComm->Size() <= 1)
    {
        return TRUE; //No MPI or lone rank: nothing splits, halo vacuously in sync.
    }

    //Arm the anti-recursion diagnostic around the actual gather/exchange:
    //a guarded launch re-entering from here is a bug (3.3).
    ++m_uiRefillInFlight;
    RefillHaloBuffer(info.m_pDeviceData, info.m_uiBytesPerSite, info.m_byFieldId);
    --m_uiRefillInFlight;
    return TRUE;
#else
    (void)buffer;
    return TRUE;
#endif
}

void CHaloManager::RefillHaloBuffer(BYTE* pFieldBytes, UINT uiBytesPerSite, BYTE byFieldId)
{
#if _CLG_MULTI_GPU
    const CLGComm* pComm = appGetComm();
    if (NULL == pComm || pComm->Size() <= 1)
    {
        return; //No MPI or lone rank: nothing splits, halo is empty.
    }

    if (0 == uiBytesPerSite || NULL == pFieldBytes)
    {
        return;
    }
    const UINT uiVolume = _HC_Volume;

    //Gather this rank's boundary sites (halo-slot order) directly into the field
    //buffer's halo tail. For the self-neighbor case this IS the final halo content;
    //for off-rank neighbors it's the send buffer and we'll overwrite with the recv.
    BYTE* pHaloTail = pFieldBytes + uiVolume * uiBytesPerSite;
    const UINT uiHaloSites = appHaloGatherLocal(pFieldBytes, pHaloTail, uiBytesPerSite);
    if (0 == uiHaloSites)
    {
        return; //Nothing split.
    }

    //Phase 1: face-only exchange (split directions, both sides). Each split
    //direction exchanges 2 face blocks (neg/pos) independently. When the
    //process-grid neighbor is this rank itself (periodic grid with size==1 in
    //that direction), the gather above already put the correct data in place
    //(self-copy), so skip MPI for that face.
    const UINT uiGrid[4] = { _HC_GpuGridX, _HC_GpuGridY, _HC_GpuGridZ, _HC_GpuGridT };
    const UINT uiLocalL[4] = { _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt };
    const UINT uiHaloWidth = _HC_HaloWidth;

    for (UINT dir = 0; dir < 4; ++dir)
    {
        if (uiGrid[dir] <= 1)
        {
            continue; //Not split in this direction.
        }

        const UINT uiFaceBlockSites = uiHaloWidth * (uiVolume / uiLocalL[dir]);
        const UINT uiFaceBlockBytes = uiFaceBlockSites * uiBytesPerSite;

        //(Re)size the four host staging buffers to this direction's face block.
        //MS-MPI is NOT CUDA-aware, so device face data must be staged host-side.
        if (m_uiCommBufferBytes < uiFaceBlockBytes)
        {
            for (UINT s = 0; s < 2; ++s)
            {
                if (NULL != m_pCommBufferSend[s]) { free(m_pCommBufferSend[s]); }
                if (NULL != m_pCommBufferRecv[s]) { free(m_pCommBufferRecv[s]); }
                m_pCommBufferSend[s] = static_cast<BYTE*>(malloc(uiFaceBlockBytes));
                m_pCommBufferRecv[s] = static_cast<BYTE*>(malloc(uiFaceBlockBytes));
            }
            m_uiCommBufferBytes = uiFaceBlockBytes;
        }

        //Non-blocking exchange of BOTH faces of this split direction at once. A
        //blocking Sendrecv per side deadlocks on a size-2 torus (both neighbours are
        //the same rank, so each side's recv waits on data the peer only sends on its
        //other side). Post all Irecv + Isend, then Waitall. Self-neighbour faces
        //(size-1 grid dim) already hold correct data from the gather -> skip.
        BYTE* pFace[2] = { NULL, NULL };
        MPI_Request reqs[4];
        UINT uiReqCount = 0;
        for (UINT side = 0; side < 2; ++side)
        {
            const UINT uiBlockOffset = _haloFaceBlockOffsetSites(uiLocalL, uiGrid, uiVolume, uiHaloWidth, dir, side);
            pFace[side] = pHaloTail + uiBlockOffset * uiBytesPerSite;

            const INT iSign = (0 == side) ? -1 : +1;
            const UINT uiNbrRank = pComm->NeighbourRank(dir, iSign);
            if (uiNbrRank == pComm->Rank())
            {
                pFace[side] = NULL; //Self-neighbour: gathered data already correct.
                continue;
            }

            appSimpleCopyDH(m_pCommBufferSend[side], pFace[side], uiFaceBlockBytes);

            //Tag pairs my (dir,side) send with the peer's (dir,1-side) send: my -dir
            //face lands in the peer's +dir halo and vice versa. Both directions of a
            //pair therefore carry distinct tags, so the two concurrent exchanges with
            //the SAME rank (size-2 torus) never collide.
            const INT iSendTag = static_cast<INT>(byFieldId * 10 + dir * 2 + side);
            const INT iRecvTag = static_cast<INT>(byFieldId * 10 + dir * 2 + (1 - side));

            MPI_Irecv(m_pCommBufferRecv[side], static_cast<INT>(uiFaceBlockBytes), MPI_BYTE,
                uiNbrRank, iRecvTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
            MPI_Isend(m_pCommBufferSend[side], static_cast<INT>(uiFaceBlockBytes), MPI_BYTE,
                uiNbrRank, iSendTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
        }

        if (uiReqCount > 0)
        {
            MPI_Waitall(static_cast<INT>(uiReqCount), reqs, MPI_STATUSES_IGNORE);
            for (UINT side = 0; side < 2; ++side)
            {
                if (NULL != pFace[side])
                {
                    appSimpleCopyHD(pFace[side], m_pCommBufferRecv[side], uiFaceBlockBytes);
                }
            }
        }
    }

    //P4-5.3: EDGE blocks -- a cell out in exactly TWO split directions. Each
    //(dir1,dir2) pair has 4 side-pairs; the peer is the rank reached by moving
    //side1 in dir1 and side2 in dir2 on the process torus (an off-diagonal
    //neighbour, never a face neighbour). Layout per CLGHaloLayout.h: dir-pair
    //(d1<d2) major, side-pair 2*s1+s2, layer1-major then layer2, edgeIdx.
    //The gather already staged this rank's near-boundary EDGE in every edge
    //slot (P4-5.2), so the exchange is the same shape as the face loop: stage
    //host-side, non-blocking send/recv, one Waitall for the whole kind.
    //Tag space: edge tags live at byFieldId*100 + 100 + pairIdx*4 + sp, above
    //the face tags (byFieldId*10 + dir*2 + side) so the same field's face and
    //edge exchanges with the SAME rank (size-2 torus) never collide.
    {
        struct SEdgeExch
        {
            BYTE* pSend;
            BYTE* pRecv;
            BYTE* pBlock;
            UINT uiBytes;
        };
        SEdgeExch exch[24]; //max 6 dir-pairs * 4 side-pairs
        UINT uiExchCount = 0;
        MPI_Request reqs[48]; //2 per exchange
        UINT uiReqCount = 0;

        UINT uiPairIdx = 0;
        for (UINT d1 = 0; d1 < 4; ++d1)
        {
            for (UINT d2 = d1 + 1; d2 < 4; ++d2)
            {
                if (uiGrid[d1] > 1 && uiGrid[d2] > 1)
                {
                    const UINT uiEdgeVol = uiVolume / (uiLocalL[d1] * uiLocalL[d2]);
                    const UINT uiBlockSites = uiHaloWidth * uiHaloWidth * uiEdgeVol;
                    const UINT uiBlockBytes = uiBlockSites * uiBytesPerSite;
                    for (UINT sp = 0; sp < 4; ++sp)
                    {
                        const UINT s1 = sp / 2;
                        const UINT s2 = sp % 2;
                        const UINT uiBlockOffset = _haloEdgeBlockOffsetSites(
                            uiLocalL, uiGrid, uiVolume, uiHaloWidth, d1, d2, s1, s2);
                        BYTE* pBlock = pHaloTail + uiBlockOffset * uiBytesPerSite;

                        //Peer = move side1 along d1 and side2 along d2 on the
                        //process-grid torus (x slowest, matching CLGComm).
                        UINT uiCoord[4];
                        for (UINT i = 0; i < 4; ++i)
                        {
                            uiCoord[i] = pComm->GridCoord()[i];
                        }
                        uiCoord[d1] = (uiCoord[d1] + uiGrid[d1] + static_cast<UINT>((0 == s1) ? -1 : 1)) % uiGrid[d1];
                        uiCoord[d2] = (uiCoord[d2] + uiGrid[d2] + static_cast<UINT>((0 == s2) ? -1 : 1)) % uiGrid[d2];
                        UINT uiNbrRank = 0;
                        for (UINT i = 0; i < 4; ++i)
                        {
                            uiNbrRank = uiNbrRank * uiGrid[i] + uiCoord[i];
                        }
                        if (uiNbrRank == pComm->Rank())
                        {
                            continue; //Self-neighbour: gathered data already correct.
                        }

                        //Temporary staging buffers per block: up to 24 concurrent
                        //edge exchanges, too many for the 2-slot member buffers.
                        SEdgeExch& e = exch[uiExchCount++];
                        e.pSend = static_cast<BYTE*>(malloc(uiBlockBytes));
                        e.pRecv = static_cast<BYTE*>(malloc(uiBlockBytes));
                        e.pBlock = pBlock;
                        e.uiBytes = uiBlockBytes;
                        appSimpleCopyDH(e.pSend, pBlock, uiBlockBytes);

                        //Tag mirrors the face rule: my (s1,s2) send pairs with the
                        //peer's (1-s1,1-s2) send. Both directions carry distinct
                        //tags, so concurrent exchanges with the same rank never
                        //collide.
                        const INT iSendTag = static_cast<INT>(byFieldId * 100 + 100 + uiPairIdx * 4 + sp);
                        const INT iRecvTag = static_cast<INT>(byFieldId * 100 + 100 + uiPairIdx * 4 + (3 - sp));
                        MPI_Irecv(e.pRecv, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                            uiNbrRank, iRecvTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
                        MPI_Isend(e.pSend, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                            uiNbrRank, iSendTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
                    }
                }
                ++uiPairIdx;
            }
        }

        if (uiReqCount > 0)
        {
            MPI_Waitall(static_cast<INT>(uiReqCount), reqs, MPI_STATUSES_IGNORE);
            for (UINT i = 0; i < uiExchCount; ++i)
            {
                appSimpleCopyHD(exch[i].pBlock, exch[i].pRecv, exch[i].uiBytes);
                free(exch[i].pSend);
                free(exch[i].pRecv);
            }
        }
    }

    //P4-5.3: CORNER blocks -- a cell out in exactly THREE split directions.
    //Each (dir1,dir2,dir3) triple has 8 side-triples; the peer is the rank
    //reached by moving side1/2/3 along d1/d2/d3. Layout: dir-triple major,
    //side-triple 4*s1+2*s2+s3, layer2-major, layer1, layer0, cornerIdx.
    //Tag space: corner tags at byFieldId*100 + 200 + tripleIdx*8 + st, above
    //both face and edge tags.
    {
        struct SCornerExch
        {
            BYTE* pSend;
            BYTE* pRecv;
            BYTE* pBlock;
            UINT uiBytes;
        };
        SCornerExch exch[32]; //max 4 dir-triples * 8 side-triples
        UINT uiExchCount = 0;
        MPI_Request reqs[64]; //2 per exchange
        UINT uiReqCount = 0;

        UINT uiTripleIdx = 0;
        for (UINT d1 = 0; d1 < 4; ++d1)
        {
            for (UINT d2 = d1 + 1; d2 < 4; ++d2)
            {
                for (UINT d3 = d2 + 1; d3 < 4; ++d3)
                {
                    if (uiGrid[d1] > 1 && uiGrid[d2] > 1 && uiGrid[d3] > 1)
                    {
                        const UINT uiCornerVol = uiVolume / (uiLocalL[d1] * uiLocalL[d2] * uiLocalL[d3]);
                        const UINT uiBlockSites = uiHaloWidth * uiHaloWidth * uiHaloWidth * uiCornerVol;
                        const UINT uiBlockBytes = uiBlockSites * uiBytesPerSite;
                        for (UINT st = 0; st < 8; ++st)
                        {
                            const UINT s1 = st / 4;
                            const UINT s2 = (st / 2) % 2;
                            const UINT s3 = st % 2;
                            const UINT uiBlockOffset = _haloCornerBlockOffsetSites(
                                uiLocalL, uiGrid, uiVolume, uiHaloWidth, d1, d2, d3, s1, s2, s3);
                            BYTE* pBlock = pHaloTail + uiBlockOffset * uiBytesPerSite;

                            UINT uiCoord[4];
                            for (UINT i = 0; i < 4; ++i)
                            {
                                uiCoord[i] = pComm->GridCoord()[i];
                            }
                            uiCoord[d1] = (uiCoord[d1] + uiGrid[d1] + static_cast<UINT>((0 == s1) ? -1 : 1)) % uiGrid[d1];
                            uiCoord[d2] = (uiCoord[d2] + uiGrid[d2] + static_cast<UINT>((0 == s2) ? -1 : 1)) % uiGrid[d2];
                            uiCoord[d3] = (uiCoord[d3] + uiGrid[d3] + static_cast<UINT>((0 == s3) ? -1 : 1)) % uiGrid[d3];
                            UINT uiNbrRank = 0;
                            for (UINT i = 0; i < 4; ++i)
                            {
                                uiNbrRank = uiNbrRank * uiGrid[i] + uiCoord[i];
                            }
                            if (uiNbrRank == pComm->Rank())
                            {
                                continue;
                            }

                            SCornerExch& e = exch[uiExchCount++];
                            e.pSend = static_cast<BYTE*>(malloc(uiBlockBytes));
                            e.pRecv = static_cast<BYTE*>(malloc(uiBlockBytes));
                            e.pBlock = pBlock;
                            e.uiBytes = uiBlockBytes;
                            appSimpleCopyDH(e.pSend, pBlock, uiBlockBytes);

                            const INT iSendTag = static_cast<INT>(byFieldId * 100 + 200 + uiTripleIdx * 8 + st);
                            const INT iRecvTag = static_cast<INT>(byFieldId * 100 + 200 + uiTripleIdx * 8 + (7 - st));
                            MPI_Irecv(e.pRecv, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                                uiNbrRank, iRecvTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
                            MPI_Isend(e.pSend, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                                uiNbrRank, iSendTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
                        }
                    }
                    ++uiTripleIdx;
                }
            }
        }

        if (uiReqCount > 0)
        {
            MPI_Waitall(static_cast<INT>(uiReqCount), reqs, MPI_STATUSES_IGNORE);
            for (UINT i = 0; i < uiExchCount; ++i)
            {
                appSimpleCopyHD(exch[i].pBlock, exch[i].pRecv, exch[i].uiBytes);
                free(exch[i].pSend);
                free(exch[i].pRecv);
            }
        }
    }

    //Improve-1 (3.8): HYPER-CORNER blocks -- a cell out in ALL FOUR split
    //directions. 16 side-quads (8*s0+4*s1+2*s2+s3, side of the smallest dir in
    //the MSB, continuing the edge/corner rule); the peer is the rank reached by
    //moving along all four dirs. Block = haloWidth^4 sites (no free axis).
    //Tag space: hyper tags at byFieldId*100 + 300 + sq, above face/edge/corner
    //tags; recv mirrors the side-quad (15 - sq).
    if (uiGrid[0] > 1 && uiGrid[1] > 1 && uiGrid[2] > 1 && uiGrid[3] > 1)
    {
        struct SHyperExch
        {
            BYTE* pSend;
            BYTE* pRecv;
            BYTE* pBlock;
            UINT uiBytes;
        };
        SHyperExch exch[16]; //max 16 side-quads
        UINT uiExchCount = 0;
        MPI_Request reqs[32]; //2 per exchange
        UINT uiReqCount = 0;

        const UINT uiBlockSites = _haloHyperBlockSites(uiLocalL, uiGrid, uiVolume, uiHaloWidth);
        const UINT uiBlockBytes = uiBlockSites * uiBytesPerSite;
        for (UINT sq = 0; sq < 16; ++sq)
        {
            const UINT s0 = sq / 8;
            const UINT s1 = (sq / 4) % 2;
            const UINT s2 = (sq / 2) % 2;
            const UINT s3 = sq % 2;
            const UINT uiBlockOffset = _haloHyperBlockOffsetSites(
                uiLocalL, uiGrid, uiVolume, uiHaloWidth, s0, s1, s2, s3);
            BYTE* pBlock = pHaloTail + uiBlockOffset * uiBytesPerSite;

            UINT uiCoord[4];
            for (UINT i = 0; i < 4; ++i)
            {
                uiCoord[i] = pComm->GridCoord()[i];
            }
            uiCoord[0] = (uiCoord[0] + uiGrid[0] + static_cast<UINT>((0 == s0) ? -1 : 1)) % uiGrid[0];
            uiCoord[1] = (uiCoord[1] + uiGrid[1] + static_cast<UINT>((0 == s1) ? -1 : 1)) % uiGrid[1];
            uiCoord[2] = (uiCoord[2] + uiGrid[2] + static_cast<UINT>((0 == s2) ? -1 : 1)) % uiGrid[2];
            uiCoord[3] = (uiCoord[3] + uiGrid[3] + static_cast<UINT>((0 == s3) ? -1 : 1)) % uiGrid[3];
            UINT uiNbrRank = 0;
            for (UINT i = 0; i < 4; ++i)
            {
                uiNbrRank = uiNbrRank * uiGrid[i] + uiCoord[i];
            }
            if (uiNbrRank == pComm->Rank())
            {
                continue; //Self-neighbour: gathered data already correct.
            }

            SHyperExch& e = exch[uiExchCount++];
            e.pSend = static_cast<BYTE*>(malloc(uiBlockBytes));
            e.pRecv = static_cast<BYTE*>(malloc(uiBlockBytes));
            e.pBlock = pBlock;
            e.uiBytes = uiBlockBytes;
            appSimpleCopyDH(e.pSend, pBlock, uiBlockBytes);

            const INT iSendTag = static_cast<INT>(byFieldId * 100 + 300 + static_cast<INT>(sq));
            const INT iRecvTag = static_cast<INT>(byFieldId * 100 + 300 + static_cast<INT>(15 - sq));
            MPI_Irecv(e.pRecv, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                uiNbrRank, iRecvTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
            MPI_Isend(e.pSend, static_cast<INT>(uiBlockBytes), MPI_BYTE,
                uiNbrRank, iSendTag, MPI_COMM_WORLD, &reqs[uiReqCount++]);
        }

        if (uiReqCount > 0)
        {
            MPI_Waitall(static_cast<INT>(uiReqCount), reqs, MPI_STATUSES_IGNORE);
            for (UINT i = 0; i < uiExchCount; ++i)
            {
                appSimpleCopyHD(exch[i].pBlock, exch[i].pRecv, exch[i].uiBytes);
                free(exch[i].pSend);
                free(exch[i].pRecv);
            }
        }
    }
#else
    (void)byFieldId;
#endif
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
