//=============================================================================
// FILENAME : CLGHaloLayoutRuntime.h
//
// DESCRIPTION:
// Host convenience wrappers over CLGHaloLayout.h that read this rank's lattice
// and process grid straight from the const-integer macros (_HC_Volume, _HC_Lx,
// _HC_GpuGridX, _HC_HaloWidth). Must be included AFTER Data/CCommonData.h.
//
// On single-GPU builds every grid factor is 1, so every count here is 0 and
// callers allocate exactly the pre-halo size -- bit-identical to before.
//
// REVISION:
//  [07/30/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CLGHALOLAYOUTRUNTIME_H_
#define _CLGHALOLAYOUTRUNTIME_H_

__BEGIN_NAMESPACE

//Extra SITE slots to append after the local volume for a site-indexed field.
inline UINT _HC_HaloSiteCount()
{
    const UINT uiLocalL[4] = { _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt };
    const UINT uiGrid[4] = { _HC_GpuGridX, _HC_GpuGridY, _HC_GpuGridZ, _HC_GpuGridT };
    return _haloTotalSites(uiLocalL, uiGrid, _HC_Volume, _HC_HaloWidth);
}

//Extra LINK slots for a gauge field = halo sites * Dir (Dir links per site).
inline UINT _HC_HaloLinkCount()
{
    return _HC_HaloSiteCount() * _HC_Dir;
}

__END_NAMESPACE

#endif //#ifndef _CLGHALOLAYOUTRUNTIME_H_

//=============================================================================
// END OF FILE
//=============================================================================
