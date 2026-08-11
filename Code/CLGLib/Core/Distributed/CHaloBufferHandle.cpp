//=============================================================================
// FILENAME : CHaloBufferHandle.cpp
//
// DESCRIPTION:
// See CHaloBufferHandle.h.
//
// REVISION:
//  [08/09/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CHaloBufferHandle::CHaloBufferHandle()
    : m_ullInteriorVersion(0)
    , m_ullHaloVersion(0)
    , m_uiValidWidth(0)
    , m_bBound(FALSE)
{
    memset(&m_sInfo, 0, sizeof(m_sInfo));
}

CHaloBufferHandle::~CHaloBufferHandle()
{
    if (m_bBound)
    {
        Unbind();
    }
}

const SHaloBufferInfo& CHaloBufferHandle::Info() const
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::Info: access after Unbind (or before Bind).\n"));
    }
    return m_sInfo;
}

UBOOL CHaloBufferHandle::Bind(const SHaloBufferInfo& info)
{
    if (m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::Bind: already bound, Unbind first.\n"));
        return FALSE;
    }

    //capacity >= (localSiteCount + haloSiteCount) * bytesPerSite, with an
    //overflow-safe product (size_t here; guard before multiplying).
    const size_t uiTotalSites = static_cast<size_t>(info.m_uiLocalSiteCount) + info.m_uiHaloSiteCount;
    if (0 != info.m_uiBytesPerSite
        && uiTotalSites > static_cast<size_t>(0xFFFFFFFFFFFFFFFFULL) / info.m_uiBytesPerSite)
    {
        appCrucial(_T("CHaloBufferHandle::Bind: site count * bytesPerSite overflows.\n"));
        return FALSE;
    }
    if (info.m_uiCapacityBytes < uiTotalSites * info.m_uiBytesPerSite)
    {
        appCrucial(_T("CHaloBufferHandle::Bind: capacity %llu < required %llu bytes.\n"),
            static_cast<ULONGLONG>(info.m_uiCapacityBytes),
            static_cast<ULONGLONG>(uiTotalSites * info.m_uiBytesPerSite));
        return FALSE;
    }
    if (NULL == info.m_pDeviceData && uiTotalSites > 0)
    {
        appCrucial(_T("CHaloBufferHandle::Bind: NULL data pointer with non-zero extent.\n"));
        return FALSE;
    }

    //Debug registry: reject overlapping extents (multi-GPU-improve1.md 3.1).
    if (NULL != appGetHaloManager()
        && !appGetHaloManager()->OnHandleBound(this, info))
    {
        return FALSE;
    }

    m_sInfo = info;
    m_ullInteriorVersion = 1;
    m_ullHaloVersion = 0;
    m_uiValidWidth = 0;
    m_bBound = TRUE;
    return TRUE;
}

void CHaloBufferHandle::Unbind()
{
    if (!m_bBound)
    {
        return;
    }
    if (NULL != appGetHaloManager())
    {
        appGetHaloManager()->OnHandleUnbound(this);
    }
    m_bBound = FALSE;
    m_ullInteriorVersion = 0;
    m_ullHaloVersion = 0;
    m_uiValidWidth = 0;
}

ULONGLONG CHaloBufferHandle::InteriorVersion() const
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::InteriorVersion: access while unbound.\n"));
    }
    return m_ullInteriorVersion;
}

ULONGLONG CHaloBufferHandle::HaloVersion() const
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::HaloVersion: access while unbound.\n"));
    }
    return m_ullHaloVersion;
}

UINT CHaloBufferHandle::ValidWidth() const
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::ValidWidth: access while unbound.\n"));
    }
    return m_uiValidWidth;
}

void CHaloBufferHandle::NotifyWritten()
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::NotifyWritten: access while unbound.\n"));
        return;
    }
    ++m_ullInteriorVersion;
    if (0 == m_ullInteriorVersion)
    {
        //Wrapped onto 0 / an old value: fail fast rather than silently
        //colliding with a stale "equal" comparison (2^64 writes -- a guard
        //for the unreachable).
        appCrucial(_T("CHaloBufferHandle::NotifyWritten: interior version overflow.\n"));
        _FAIL_EXIT;
    }
    m_uiValidWidth = 0;
}

void CHaloBufferHandle::NotifyHaloSynchronized(UINT uiWidth)
{
    if (!m_bBound)
    {
        appCrucial(_T("CHaloBufferHandle::NotifyHaloSynchronized: access while unbound.\n"));
        return;
    }
    m_ullHaloVersion = m_ullInteriorVersion;
    m_uiValidWidth = uiWidth;
}

void CHaloBufferHandle::SetFieldId(BYTE byFieldId)
{
    if (m_bBound)
    {
        m_sInfo.m_byFieldId = byFieldId;
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
