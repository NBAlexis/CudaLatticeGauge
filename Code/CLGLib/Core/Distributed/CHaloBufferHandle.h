//=============================================================================
// FILENAME : CHaloBufferHandle.h
//
// DESCRIPTION:
// Buffer-extent halo identity for the multi-GPU FIELD protocol
// (multi-GPU-improve1.md section 3.1).
//
// Every device buffer extent managed by the FIELD halo protocol owns exactly
// ONE CHaloBufferHandle; views over the same/overlapping extent must share it
// (or fail to Bind), while disjoint, independently writable pool slices each
// hold their own. The handle is owned by the field (or temporary buffer
// owner) that owns the extent and must not outlive it.
//
// The handle carries:
//   - SHaloBufferInfo: pointer/capacity/stride/counts/layout generation --
//     the ONLY source the halo manager may use for pack/unpack bounds
//     (capacity is validated at Bind, never trusted as "expected").
//   - version state: interior version (bumped on every write), halo version
//     (set equal to interior ONLY after a fully successful refill) and the
//     currently valid halo width. Consumers compare versions for EQUALITY
//     only; they must not rely on initial values or increment sizes.
//
// CHaloManager keeps a debug registry of all bound handles and rejects two
// simultaneously bound handles declaring overlapping extents.
//
// REVISION:
//  [08/09/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CHALOBUFFERHANDLE_H_
#define _CHALOBUFFERHANDLE_H_

__BEGIN_NAMESPACE

/**
* Static description of one managed device buffer extent. Filled by the owner
* at Bind time; m_byFieldId is retained for index/BC/MPI tagging only and is
* never used as buffer identity.
*/
struct CLGAPI SHaloBufferInfo
{
    BYTE* m_pDeviceData;
    size_t m_uiCapacityBytes;
    UINT  m_uiBytesPerSite;
    UINT  m_uiLocalSiteCount;
    UINT  m_uiHaloSiteCount;
    ULONGLONG m_ullLayoutGeneration;
    BYTE  m_byFieldId;
    UBOOL m_bHaloCapable;
};

class CLGAPI CHaloBufferHandle
{
public:

    CHaloBufferHandle();
    ~CHaloBufferHandle();

    UBOOL IsBound() const { return m_bBound; }

    /**
    * The validated extent description. Fail-fast when unbound: after Unbind
    * any read/write access through the handle is a bug.
    */
    const SHaloBufferInfo& Info() const;

    /**
    * Attach to an extent. Validates
    *   capacityBytes >= (localSiteCount + haloSiteCount) * bytesPerSite
    * (with an overflow-safe product) and registers with the halo manager,
    * which rejects extents overlapping any other bound handle. Returns FALSE
    * on rejection (no state changed). On success interior version is non-zero
    * and halo version / valid width are invalid.
    */
    UBOOL Bind(const SHaloBufferInfo& info);

    /** Detach; deregisters. Access after Unbind fails fast. */
    void Unbind();

    ULONGLONG InteriorVersion() const;
    ULONGLONG HaloVersion() const;
    UINT ValidWidth() const;

    /**
    * Monotonically bump the interior version and invalidate the halo
    * (valid width -> 0). A counter that wraps to 0 is a hard error
    * (fail-fast) rather than silently colliding with an old "equal" state.
    */
    void NotifyWritten();

    /**
    * May only be called after a refill has COMPLETED successfully; a failed
    * refill must never update the halo version.
    */
    void NotifyHaloSynchronized(UINT uiWidth);

    /**
    * The field id is assigned after construction (CLGLibManager); keep the
    * info block in sync. Identity (pointer/extent/version) is unchanged.
    */
    void SetFieldId(BYTE byFieldId);

protected:

    SHaloBufferInfo m_sInfo;
    ULONGLONG m_ullInteriorVersion;
    ULONGLONG m_ullHaloVersion;
    UINT m_uiValidWidth;
    UBOOL m_bBound;
};

//Declared here (defined as an inline accessor over GCLGManager in
//Core/CLGLibManager.h) so field template bases can snapshot the layout
//generation at Bind time -- the same convention as appGetHaloManager() in
//CHaloManager.h. Template definitions need a visible declaration at the
//point of use.
inline ULONGLONG appGetLayoutGeneration();

__END_NAMESPACE

#endif //#ifndef _CHALOBUFFERHANDLE_H_

//=============================================================================
// END OF FILE
//=============================================================================
