//=============================================================================
// FILENAME : CLGHaloLayout.h
//
// DESCRIPTION:
// The ONE place that defines how halo storage is laid out. Both the index
// baking (which writes an out-of-lattice neighbour's SIndex.m_uiSiteIndex to a
// halo slot) and the halo pack/unpack (which fills those slots from a
// neighbour rank) include this header, so the two can never disagree on the
// numbering. See Docs/MultiGPU-Plan.md section 10, "Phase 1 halo storage".
//
// Layout (Design B -- extended flat field buffer):
//   [ 0 .. Volume-1 ]                : the local sub-lattice, unchanged.
//   [ Volume .. Volume+haloSites-1 ] : halo sites, one slot per block as
//                                      numbered below.
//
// Block families, in order (Improve-1 3.8: codimension 1..4 under one rule):
//   FACE  (1 split dir out):  dir-major, side, layer, faceIdx.
//   EDGE  (2 split dirs out): dir-pair-major, side-pair 2*s1+s2, layer1-major,
//                             layer2, edgeIdx.
//   CORNER(3 split dirs out): dir-triple-major, side-triple 4*s1+2*s2+s3,
//                             layer2-major, layer1, layer0, cornerIdx.
//   HYPER (4 split dirs out): side-quad 8*s0+4*s1+2*s2+s3, layer3-major,
//                             layer2, layer1, layer0 (no free axis).
//
// Only SPLIT directions (GpuGrid[dir] > 1) get a face block. A non-split
// direction contributes zero slots (its out-of-lattice neighbour periodic-wraps
// to a local site, exactly as single-GPU). This keeps single-GPU builds at
// haloSites == 0, so the field buffer size and every index are bit-identical.
//
// A site-field halo slot is a site index; a gauge (link) field multiplies by
// _?C_Dir and adds the link direction, matching _deviceGetLinkIndex.
//
// REVISION:
//  [07/30/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CLGHALOLAYOUT_H_
#define _CLGHALOLAYOUT_H_

__BEGIN_NAMESPACE

//Face volume in direction dir = local Volume / local L[dir]: the number of sites
//on one (dir=const) hyperplane of this rank's sub-lattice.
//pLocalL is the 4 local lattice lengths [Lx,Ly,Lz,Lt]; uiLocalVolume their product.
inline UINT _haloFaceVolume(const UINT* pLocalL, UINT uiLocalVolume, UINT dir)
{
    return (0 == pLocalL[dir]) ? 0 : (uiLocalVolume / pLocalL[dir]);
}

//Edge volume in directions (dir1,dir2) = local Volume / (L[dir1]*L[dir2]): the
//number of sites on the (dir1=const,dir2=const) co-dimension-2 hyperplane.
inline UINT _haloEdgeVolume(const UINT* pLocalL, UINT uiLocalVolume, UINT dir1, UINT dir2)
{
    return (0 == pLocalL[dir1] || 0 == pLocalL[dir2])
        ? 0 : (uiLocalVolume / (pLocalL[dir1] * pLocalL[dir2]));
}

//Corner volume in directions (dir1,dir2,dir3) = local Volume / (L1*L2*L3): the
//number of sites on the co-dimension-3 hyperplane.
inline UINT _haloCornerVolume(const UINT* pLocalL, UINT uiLocalVolume, UINT dir1, UINT dir2, UINT dir3)
{
    return (0 == pLocalL[dir1] || 0 == pLocalL[dir2] || 0 == pLocalL[dir3])
        ? 0 : (uiLocalVolume / (pLocalL[dir1] * pLocalL[dir2] * pLocalL[dir3]));
}

//Number of halo SITE slots for one (dir,side) face block:
//  haloWidth layers, each holding one face volume of sites.
//Zero when the direction is not split.
inline UINT _haloFaceBlockSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir)
{
    if (pGrid[dir] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * _haloFaceVolume(pLocalL, uiLocalVolume, dir);
}

//Number of halo SITE slots for one (dir1,dir2,side1,side2) edge block:
//  haloWidth^2 layers, each holding one edge volume of sites.
//Zero unless BOTH directions are split.
inline UINT _haloEdgeBlockSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2)
{
    if (pGrid[dir1] <= 1 || pGrid[dir2] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * _haloEdgeVolume(pLocalL, uiLocalVolume, dir1, dir2);
}

//Number of halo SITE slots for one (dir1,dir2,dir3,side1,side2,side3) corner
//block: haloWidth^3 layers, each holding one corner volume of sites.
//Zero unless ALL THREE directions are split.
inline UINT _haloCornerBlockSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2, UINT dir3)
{
    if (pGrid[dir1] <= 1 || pGrid[dir2] <= 1 || pGrid[dir3] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * uiHaloWidth
        * _haloCornerVolume(pLocalL, uiLocalVolume, dir1, dir2, dir3);
}

//Starting SITE slot (relative to the first halo slot, i.e. relative to Volume)
//of face block (dir, side). Blocks are ordered dir-major then side (0=neg,1=pos).
//This is the exact prefix sum the bake and the pack/unpack must share.
inline UINT _haloFaceBlockOffsetSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir, UINT side)
{
    UINT uiOffset = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        const UINT uiBlock = _haloFaceBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d);
        if (d < dir)
        {
            uiOffset += 2 * uiBlock; //both sides of every earlier direction
        }
        else if (d == dir)
        {
            uiOffset += side * uiBlock; //neg side precedes pos side within this direction
            break;
        }
    }
    return uiOffset;
}

//Total halo SITE count across all split directions (both sides), FACE blocks only.
//Kept separate so the P4-5 edge/corner blocks can be appended after all faces
//without touching the face numbering.
inline UINT _haloTotalFaceSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        uiTotal += 2 * _haloFaceBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d);
    }
    return uiTotal;
}

//Starting SITE slot (relative to Volume) of edge block (dir1,dir2,side1,side2).
//Edge blocks are ordered after ALL face blocks, then dir-pair-major
//(d1 < d2 ascending), then (side1,side2) with side-pair = 2*side1+side2.
//This is the exact prefix sum the bake and the pack/unpack must share.
inline UINT _haloEdgeBlockOffsetSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2, UINT side1, UINT side2)
{
    UINT uiOffset = _haloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            const UINT uiBlock = _haloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
            if (d1 < dir1 || (d1 == dir1 && d2 < dir2))
            {
                uiOffset += 4 * uiBlock; //both sides of both directions
            }
            else if (d1 == dir1 && d2 == dir2)
            {
                uiOffset += (2 * side1 + side2) * uiBlock;
                break;
            }
        }
        if (d1 >= dir1)
        {
            break;
        }
    }
    return uiOffset;
}

//Total halo SITE count of all edge blocks (corner blocks are ordered after
//ALL face + edge blocks, so this is part of the corner prefix sum).
inline UINT _haloTotalEdgeSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = 0;
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            uiTotal += 4 * _haloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
        }
    }
    return uiTotal;
}

//Starting SITE slot (relative to Volume) of corner block
//(dir1,dir2,dir3,side1,side2,side3). Corner blocks are ordered after ALL face +
//edge blocks, then dir-triple-major (d1 < d2 < d3 ascending), then
//(side1,side2,side3) with side-triple = 4*side1 + 2*side2 + side3.
//This is the exact prefix sum the bake and the pack/unpack must share.
inline UINT _haloCornerBlockOffsetSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2, UINT dir3,
    UINT side1, UINT side2, UINT side3)
{
    UINT uiOffset = _haloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + _haloTotalEdgeSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            for (UINT d3 = d2 + 1; d3 < 4; ++d3)
            {
                const UINT uiBlock = _haloCornerBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2, d3);
                if (d1 < dir1 || (d1 == dir1 && d2 < dir2) || (d1 == dir1 && d2 == dir2 && d3 < dir3))
                {
                    uiOffset += 8 * uiBlock; //both sides of all three directions
                }
                else if (d1 == dir1 && d2 == dir2 && d3 == dir3)
                {
                    uiOffset += (4 * side1 + 2 * side2 + side3) * uiBlock;
                    break;
                }
            }
            if (d1 == dir1 && d2 >= dir2)
            {
                break;
            }
        }
        if (d1 >= dir1)
        {
            break;
        }
    }
    return uiOffset;
}

//Improve-1 (3.8): hyper-corner volume = sites on the co-dimension-4
//hyperplane -- exactly one (no free axis), but keep the formula explicit for
//symmetry with the lower codimensions.
inline UINT _haloHyperVolume(const UINT* pLocalL, UINT uiLocalVolume)
{
    return (0 == pLocalL[0] || 0 == pLocalL[1] || 0 == pLocalL[2] || 0 == pLocalL[3])
        ? 0 : (uiLocalVolume / (pLocalL[0] * pLocalL[1] * pLocalL[2] * pLocalL[3]));
}

//Number of halo SITE slots for one (side0,side1,side2,side3) hyper-corner
//block: haloWidth^4 layers of one site each (hyper-corner volume is 1).
//Zero unless ALL FOUR directions are split.
inline UINT _haloHyperBlockSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth)
{
    if (pGrid[0] <= 1 || pGrid[1] <= 1 || pGrid[2] <= 1 || pGrid[3] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * uiHaloWidth * uiHaloWidth
        * _haloHyperVolume(pLocalL, uiLocalVolume);
}

//Total halo SITE count of face + edge + corner blocks -- everything before the
//hyper-corner blocks. _haloHyperBlockOffsetSites uses this as its base.
inline UINT _haloTotalFaceEdgeCornerSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = _haloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            uiTotal += 4 * _haloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
        }
    }
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            for (UINT d3 = d2 + 1; d3 < 4; ++d3)
            {
                uiTotal += 8 * _haloCornerBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2, d3);
            }
        }
    }
    return uiTotal;
}

//Starting SITE slot (relative to Volume) of hyper-corner block
//(side0,side1,side2,side3). Hyper blocks are ordered after ALL face + edge +
//corner blocks, then side-quad = 8*side0 + 4*side1 + 2*side2 + side3 (side0 is
//the side of the SMALLEST direction, continuing the edge/corner rule).
//This is the exact prefix sum the bake and the pack/unpack must share.
inline UINT _haloHyperBlockOffsetSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth,
    UINT side0, UINT side1, UINT side2, UINT side3)
{
    return _haloTotalFaceEdgeCornerSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + (8 * side0 + 4 * side1 + 2 * side2 + side3)
        * _haloHyperBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
}

//Total halo SITE count across all split directions (both sides), including the
//P4-5 edge + corner blocks and the Improve-1 (3.8) hyper-corner blocks
//appended after the corners.
inline UINT _haloTotalSites(const UINT* pLocalL, const UINT* pGrid,
    UINT uiLocalVolume, UINT uiHaloWidth)
{
    return _haloTotalFaceEdgeCornerSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + 16 * _haloHyperBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
}

//Improve-1 (3.8): three-state result of the halo redirect decision.
//  EHR_LocalOrNonSplitWrap -- interior cell, or out of range only in NON-split
//      directions: the caller applies this field's normal local boundary
//      handling (periodic wrap / boundary condition), exactly as single-GPU.
//  EHR_Halo -- out of range in one or more SPLIT directions, every crossed
//      layer within haloWidth: uiSiteSlotOut receives the halo SITE slot
//      (relative to Volume); codimension 1..4 share one numbering rule.
//  EHR_Invalid -- out of range in a SPLIT direction BEYOND haloWidth: no legal
//      target. This must NEVER be degraded to the silent local wrap; bakes
//      write the invalid SIndex sentinel so any later use fails loudly in
//      debug instead of reading wrong data.
enum EHaloRedirectResult
{
    EHR_LocalOrNonSplitWrap,
    EHR_Halo,
    EHR_Invalid,
};

//---------------------------------------------------------------------------
// Device mirrors. Identical formulas to the host inlines above, but callable
// from bake / pack / unpack kernels. Kept pure (all inputs explicit) so this
// header has no include-order dependency on the constant-memory macros.
//
// PHASE 1 SCOPE: face-only halo (a cell out of the local lattice in EXACTLY one
// split direction). Corner/edge cells (out in two or more split directions) are
// only read by staple/plaquette/force kernels, which Phase 1 does not exercise
// (acceptance = decomposition + gather/scatter + self-exchange, no force/solver
// -- Docs/MultiGPU-Plan.md Phase 1). Those cells keep the periodic wrap for now;
// filling them (sequential face shifts) is a Phase 2 item.
//
// P4-5 SCOPE: _deviceHaloRedirectSite and the pack/unpack are extended to cover
// EDGE (exactly two split dirs) and CORNER (exactly three split dirs) cells,
// whose slots live after the face blocks per the host layout above.
//
// Improve-1 (3.8) SCOPE: the redirect decision is three-state
// (EHaloRedirectResult) and covers the HYPER-CORNER (all four split dirs out)
// as a fourth block family after the corners, so codimension 1/2/3/4 share one
// numbering rule. A split-direction crossing BEYOND haloWidth is EHR_Invalid
// and must never be silently wrapped.
//---------------------------------------------------------------------------

#ifdef __CUDACC__

__device__ __inline__ static UINT _deviceHaloFaceVolume(
    const UINT* pLocalL, UINT uiLocalVolume, UINT dir)
{
    return (0 == pLocalL[dir]) ? 0 : (uiLocalVolume / pLocalL[dir]);
}

__device__ __inline__ static UINT _deviceHaloEdgeVolume(
    const UINT* pLocalL, UINT uiLocalVolume, UINT dir1, UINT dir2)
{
    return (0 == pLocalL[dir1] || 0 == pLocalL[dir2])
        ? 0 : (uiLocalVolume / (pLocalL[dir1] * pLocalL[dir2]));
}

__device__ __inline__ static UINT _deviceHaloCornerVolume(
    const UINT* pLocalL, UINT uiLocalVolume, UINT dir1, UINT dir2, UINT dir3)
{
    return (0 == pLocalL[dir1] || 0 == pLocalL[dir2] || 0 == pLocalL[dir3])
        ? 0 : (uiLocalVolume / (pLocalL[dir1] * pLocalL[dir2] * pLocalL[dir3]));
}

__device__ __inline__ static UINT _deviceHaloFaceBlockSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth, UINT dir)
{
    if (pGrid[dir] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * _deviceHaloFaceVolume(pLocalL, uiLocalVolume, dir);
}

__device__ __inline__ static UINT _deviceHaloEdgeBlockSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2)
{
    if (pGrid[dir1] <= 1 || pGrid[dir2] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * _deviceHaloEdgeVolume(pLocalL, uiLocalVolume, dir1, dir2);
}

__device__ __inline__ static UINT _deviceHaloCornerBlockSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth, UINT dir1, UINT dir2, UINT dir3)
{
    if (pGrid[dir1] <= 1 || pGrid[dir2] <= 1 || pGrid[dir3] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * uiHaloWidth
        * _deviceHaloCornerVolume(pLocalL, uiLocalVolume, dir1, dir2, dir3);
}

__device__ __inline__ static UINT _deviceHaloFaceBlockOffsetSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    UINT dir, UINT side)
{
    UINT uiOffset = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        const UINT uiBlock = _deviceHaloFaceBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d);
        if (d < dir)
        {
            uiOffset += 2 * uiBlock;
        }
        else if (d == dir)
        {
            uiOffset += side * uiBlock;
            break;
        }
    }
    return uiOffset;
}

__device__ __inline__ static UINT _deviceHaloTotalFaceSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        uiTotal += 2 * _deviceHaloFaceBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d);
    }
    return uiTotal;
}

__device__ __inline__ static UINT _deviceHaloEdgeBlockOffsetSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    UINT dir1, UINT dir2, UINT side1, UINT side2)
{
    UINT uiOffset = _deviceHaloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            const UINT uiBlock = _deviceHaloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
            if (d1 < dir1 || (d1 == dir1 && d2 < dir2))
            {
                uiOffset += 4 * uiBlock;
            }
            else if (d1 == dir1 && d2 == dir2)
            {
                uiOffset += (2 * side1 + side2) * uiBlock;
                break;
            }
        }
        if (d1 >= dir1)
        {
            break;
        }
    }
    return uiOffset;
}

__device__ __inline__ static UINT _deviceHaloTotalEdgeSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = 0;
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            uiTotal += 4 * _deviceHaloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
        }
    }
    return uiTotal;
}

__device__ __inline__ static UINT _deviceHaloCornerBlockOffsetSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    UINT dir1, UINT dir2, UINT dir3, UINT side1, UINT side2, UINT side3)
{
    UINT uiOffset = _deviceHaloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + _deviceHaloTotalEdgeSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            for (UINT d3 = d2 + 1; d3 < 4; ++d3)
            {
                const UINT uiBlock = _deviceHaloCornerBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2, d3);
                if (d1 < dir1 || (d1 == dir1 && d2 < dir2) || (d1 == dir1 && d2 == dir2 && d3 < dir3))
                {
                    uiOffset += 8 * uiBlock;
                }
                else if (d1 == dir1 && d2 == dir2 && d3 == dir3)
                {
                    uiOffset += (4 * side1 + 2 * side2 + side3) * uiBlock;
                    break;
                }
            }
            if (d1 == dir1 && d2 >= dir2)
            {
                break;
            }
        }
        if (d1 >= dir1)
        {
            break;
        }
    }
    return uiOffset;
}

//Improve-1 (3.8) device mirrors of the hyper-corner helpers.
__device__ __inline__ static UINT _deviceHaloHyperVolume(
    const UINT* pLocalL, UINT uiLocalVolume)
{
    return (0 == pLocalL[0] || 0 == pLocalL[1] || 0 == pLocalL[2] || 0 == pLocalL[3])
        ? 0 : (uiLocalVolume / (pLocalL[0] * pLocalL[1] * pLocalL[2] * pLocalL[3]));
}

__device__ __inline__ static UINT _deviceHaloHyperBlockSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth)
{
    if (pGrid[0] <= 1 || pGrid[1] <= 1 || pGrid[2] <= 1 || pGrid[3] <= 1)
    {
        return 0;
    }
    return uiHaloWidth * uiHaloWidth * uiHaloWidth * uiHaloWidth
        * _deviceHaloHyperVolume(pLocalL, uiLocalVolume);
}

__device__ __inline__ static UINT _deviceHaloTotalFaceEdgeCornerSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth)
{
    UINT uiTotal = _deviceHaloTotalFaceSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            uiTotal += 4 * _deviceHaloEdgeBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2);
        }
    }
    for (UINT d1 = 0; d1 < 4; ++d1)
    {
        for (UINT d2 = d1 + 1; d2 < 4; ++d2)
        {
            for (UINT d3 = d2 + 1; d3 < 4; ++d3)
            {
                uiTotal += 8 * _deviceHaloCornerBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth, d1, d2, d3);
            }
        }
    }
    return uiTotal;
}

__device__ __inline__ static UINT _deviceHaloHyperBlockOffsetSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    UINT side0, UINT side1, UINT side2, UINT side3)
{
    return _deviceHaloTotalFaceEdgeCornerSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + (8 * side0 + 4 * side1 + 2 * side2 + side3)
        * _deviceHaloHyperBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
}

//Total halo SITE count across all split directions (device mirror), including
//the P4-5 edge + corner blocks and the Improve-1 (3.8) hyper-corner blocks
//appended after the corners.
__device__ __inline__ static UINT _deviceHaloTotalSites(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth)
{
    return _deviceHaloTotalFaceEdgeCornerSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth)
        + 16 * _deviceHaloHyperBlockSites(pLocalL, pGrid, uiLocalVolume, uiHaloWidth);
}

//Face-internal linear index: the position of a site within a (tdir = const)
//hyperplane of the local sub-lattice, enumerating the OTHER three axes in
//ascending axis order (tdir skipped), earliest axis slowest. Both the bake and
//the pack/unpack compute this identically so a halo slot means the same site on
//both ends. pLocalL is the 4 local lengths; pCoord the 4 in-range local coords.
__device__ __inline__ static UINT _deviceHaloFaceIndex(
    const UINT* pLocalL, const INT* pCoord, UINT tdir)
{
    UINT uiIdx = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        if (d == tdir)
        {
            continue;
        }
        uiIdx = uiIdx * pLocalL[d] + static_cast<UINT>(pCoord[d]);
    }
    return uiIdx;
}

//Edge-internal linear index: the position of a site within the co-dimension-2
//hyperplane (tdir1=const, tdir2=const), enumerating the OTHER two axes in
//ascending order (both skipped), earliest axis slowest.
__device__ __inline__ static UINT _deviceHaloEdgeIndex(
    const UINT* pLocalL, const INT* pCoord, UINT tdir1, UINT tdir2)
{
    UINT uiIdx = 0;
    for (UINT d = 0; d < 4; ++d)
    {
        if (d == tdir1 || d == tdir2)
        {
            continue;
        }
        uiIdx = uiIdx * pLocalL[d] + static_cast<UINT>(pCoord[d]);
    }
    return uiIdx;
}

//Corner-internal linear index: the position of a site within the co-dimension-3
//hyperplane (tdir1=tdir2=tdir3=const), i.e. along the single remaining axis.
__device__ __inline__ static UINT _deviceHaloCornerIndex(
    const UINT* pLocalL, const INT* pCoord, UINT tdir1, UINT tdir2, UINT tdir3)
{
    for (UINT d = 0; d < 4; ++d)
    {
        if (d == tdir1 || d == tdir2 || d == tdir3)
        {
            continue;
        }
        return static_cast<UINT>(pCoord[d]);
    }
    return 0;
}

//Decide what an out-of-lattice neighbour cell resolves to (Improve-1 3.8,
//three-state):
//  - EHR_LocalOrNonSplitWrap: interior, or out of range only in NON-split
//    directions -- the caller keeps this field's normal local wrap/boundary.
//  - EHR_Halo: out of range in one or more SPLIT directions, every crossed
//    layer within haloWidth -- uiSiteSlotOut receives the halo SITE slot
//    (relative to Volume). Codimension 1 (face), 2 (edge), 3 (corner) and 4
//    (hyper-corner) all share the one numbering rule of this header.
//  - EHR_Invalid: a SPLIT direction crossed BEYOND haloWidth -- no legal
//    target, must NEVER be degraded to the silent local wrap (the historical
//    FALSE return did exactly that for >haloWidth cells).
//pRawCoord is the signed neighbour coordinate BEFORE wrapping; pWrappedCoord
//is every axis periodic-wrapped into local range.
__device__ __inline__ static EHaloRedirectResult _deviceHaloRedirectSite(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    const INT* pRawCoord, const INT* pWrappedCoord, UINT& uiSiteSlotOut)
{
    //Collect every split direction that is out of range, with side + layer.
    //(The scan walks d ascending, so the collected dirs are already sorted.
    //The fourth dir needs no slot of its own: the hyper-corner branch below is
    //only reached when ALL FOUR dirs are out, i.e. dirs are 0,1,2,3.)
    UINT tdir0 = 0, tdir1 = 0, tdir2 = 0;
    UINT tside0 = 0, tside1 = 0, tside2 = 0, tside3 = 0;
    UINT tlayer0 = 0, tlayer1 = 0, tlayer2 = 0, tlayer3 = 0;
    UINT uiSplitOutCount = 0;

    for (UINT d = 0; d < 4; ++d)
    {
        const INT iL = static_cast<INT>(pLocalL[d]);
        const UBOOL bSplit = (pGrid[d] > 1);
        if (pRawCoord[d] < 0)
        {
            if (bSplit)
            {
                if (0 == uiSplitOutCount) { tdir0 = d; tside0 = 0; tlayer0 = static_cast<UINT>(-1 - pRawCoord[d]); }
                else if (1 == uiSplitOutCount) { tdir1 = d; tside1 = 0; tlayer1 = static_cast<UINT>(-1 - pRawCoord[d]); }
                else if (2 == uiSplitOutCount) { tdir2 = d; tside2 = 0; tlayer2 = static_cast<UINT>(-1 - pRawCoord[d]); }
                else { tside3 = 0; tlayer3 = static_cast<UINT>(-1 - pRawCoord[d]); }
                ++uiSplitOutCount;
            }
        }
        else if (pRawCoord[d] >= iL)
        {
            if (bSplit)
            {
                if (0 == uiSplitOutCount) { tdir0 = d; tside0 = 1; tlayer0 = static_cast<UINT>(pRawCoord[d] - iL); }
                else if (1 == uiSplitOutCount) { tdir1 = d; tside1 = 1; tlayer1 = static_cast<UINT>(pRawCoord[d] - iL); }
                else if (2 == uiSplitOutCount) { tdir2 = d; tside2 = 1; tlayer2 = static_cast<UINT>(pRawCoord[d] - iL); }
                else { tside3 = 1; tlayer3 = static_cast<UINT>(pRawCoord[d] - iL); }
                ++uiSplitOutCount;
            }
        }
    }

    //A cell is representable only when every out-axis layer is within haloWidth.
    UBOOL bInWidth = TRUE;
    if (uiSplitOutCount >= 1 && tlayer0 >= uiHaloWidth) { bInWidth = FALSE; }
    if (uiSplitOutCount >= 2 && tlayer1 >= uiHaloWidth) { bInWidth = FALSE; }
    if (uiSplitOutCount >= 3 && tlayer2 >= uiHaloWidth) { bInWidth = FALSE; }
    if (uiSplitOutCount >= 4 && tlayer3 >= uiHaloWidth) { bInWidth = FALSE; }

    if (0 == uiSplitOutCount)
    {
        return EHR_LocalOrNonSplitWrap;
    }
    if (!bInWidth)
    {
        //Improve-1 (3.8): split-direction crossing beyond haloWidth has NO
        //legal target -- the caller writes the invalid sentinel, never a wrap.
        return EHR_Invalid;
    }

    if (1 == uiSplitOutCount)
    {
        const UINT uiFaceVol = _deviceHaloFaceVolume(pLocalL, uiLocalVolume, tdir0);
        const UINT uiFaceIdx = _deviceHaloFaceIndex(pLocalL, pWrappedCoord, tdir0);
        const UINT uiBlockOffset = _deviceHaloFaceBlockOffsetSites(
            pLocalL, pGrid, uiLocalVolume, uiHaloWidth, tdir0, tside0);
        uiSiteSlotOut = uiBlockOffset + tlayer0 * uiFaceVol + uiFaceIdx;
        return EHR_Halo;
    }

    if (2 == uiSplitOutCount)
    {
        const UINT uiEdgeVol = _deviceHaloEdgeVolume(pLocalL, uiLocalVolume, tdir0, tdir1);
        const UINT uiEdgeIdx = _deviceHaloEdgeIndex(pLocalL, pWrappedCoord, tdir0, tdir1);
        const UINT uiBlockOffset = _deviceHaloEdgeBlockOffsetSites(
            pLocalL, pGrid, uiLocalVolume, uiHaloWidth, tdir0, tdir1, tside0, tside1);
        //layer1-major within the block (layer1 changes slowest), matching the
        //pack/unpack convention chosen for edge blocks.
        uiSiteSlotOut = uiBlockOffset + tlayer1 * uiHaloWidth * uiEdgeVol + tlayer0 * uiEdgeVol + uiEdgeIdx;
        return EHR_Halo;
    }

    if (3 == uiSplitOutCount)
    {
        const UINT uiCornerVol = _deviceHaloCornerVolume(pLocalL, uiLocalVolume, tdir0, tdir1, tdir2);
        const UINT uiCornerIdx = _deviceHaloCornerIndex(pLocalL, pWrappedCoord, tdir0, tdir1, tdir2);
        const UINT uiBlockOffset = _deviceHaloCornerBlockOffsetSites(
            pLocalL, pGrid, uiLocalVolume, uiHaloWidth, tdir0, tdir1, tdir2, tside0, tside1, tside2);
        //layer2-major, then layer1, then layer0 within the block.
        uiSiteSlotOut = uiBlockOffset
            + tlayer2 * uiHaloWidth * uiHaloWidth * uiCornerVol
            + tlayer1 * uiHaloWidth * uiCornerVol
            + tlayer0 * uiCornerVol
            + uiCornerIdx;
        return EHR_Halo;
    }

    //4 == uiSplitOutCount (hyper-corner): side-quad 8*s0+4*s1+2*s2+s3,
    //layer3-major then layer2, layer1, layer0; the hyperplane has no free axis.
    {
        const UINT uiBlockOffset = _deviceHaloHyperBlockOffsetSites(
            pLocalL, pGrid, uiLocalVolume, uiHaloWidth, tside0, tside1, tside2, tside3);
        uiSiteSlotOut = uiBlockOffset
            + ((tlayer3 * uiHaloWidth + tlayer2) * uiHaloWidth + tlayer1) * uiHaloWidth
            + tlayer0;
        return EHR_Halo;
    }
}

//Multi-GPU, boundary-aware halo redirect: like _deviceHaloRedirectSite, but a
//crossing at the GLOBAL edge of a split direction (this rank owns the physical
//boundary there) is NOT redirected -- the boundary condition (Dirichlet clamp /
//projective reflection) applies at that edge instead. Pass pGlobalOffset /
//pGlobalLen (ECI_GlobalOffsetX+i / ECI_GlobalLx+i) for the rank's offset and the
//global lattice. Used by the non-torus boundary bakes so a decomposition never
//silently wraps the physical boundary (halo exchange is torus-only).
//Improve-1 (3.8): three-state result; EHR_Invalid still means "split-direction
//crossing beyond haloWidth" after the global-edge clamp and must never wrap.
__device__ __inline__ static EHaloRedirectResult _deviceHaloRedirectSiteBoundaryAware(
    const UINT* pLocalL, const UINT* pGrid, UINT uiLocalVolume, UINT uiHaloWidth,
    const INT* pRawCoord, const INT* pWrappedCoord,
    const INT* pGlobalOffset, const INT* pGlobalLen,
    const BYTE* pBoundaryDirMask, UINT& uiSiteSlotOut)
{
    //pBoundaryDirMask[d] != 0 marks directions whose boundary condition applies
    //at the GLOBAL edge (Dirichlet clamp / projective reflection). Only those
    //directions must NOT be redirected to halo when this rank owns the edge;
    //periodic/anti-periodic directions still redirect (the neighbour rank is
    //the torus wrap), so a decomposition never loses the wrap data.
    INT iRawAdj[4] = { pRawCoord[0], pRawCoord[1], pRawCoord[2], pRawCoord[3] };
    for (UINT d = 0; d < 4; ++d)
    {
        if (0 == pBoundaryDirMask[d])
        {
            continue;
        }
        //Negative-side crossing on the rank that owns the global negative edge.
        if (iRawAdj[d] < 0 && 0 == pGlobalOffset[d])
        {
            iRawAdj[d] = 0;
        }
        //Positive-side crossing on the rank that owns the global positive edge.
        if (iRawAdj[d] >= static_cast<INT>(pLocalL[d])
            && (pGlobalOffset[d] + static_cast<INT>(pLocalL[d])) == pGlobalLen[d])
        {
            iRawAdj[d] = static_cast<INT>(pLocalL[d]) - 1;
        }
    }
    return _deviceHaloRedirectSite(pLocalL, pGrid, uiLocalVolume, uiHaloWidth,
        iRawAdj, pWrappedCoord, uiSiteSlotOut);
}


#endif //#ifdef __CUDACC__

//The host convenience wrappers that read the const-integer macros (_HC_Volume
//etc.) live in CLGHaloLayoutRuntime.h, included after CCommonData.h. They can't
//live here because this header must precede CudaHelper.h (include-order needs of
//the launch guard) while CCommonData.h comes later.

__END_NAMESPACE

#endif //#ifndef _CLGHALOLAYOUT_H_

//=============================================================================
// END OF FILE
//=============================================================================
