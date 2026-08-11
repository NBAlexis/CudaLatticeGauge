//=============================================================================
// FILENAME : CGaugeFixing.cpp
//
// DESCRIPTION:
// Base class helpers for the multi-GPU gauge-fixing path (P4-2.1/P4-2.2):
// rank0 global-lattice context switch for the gather -> fix -> scatter flow.
//
// REVISION:
//  [08/03/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CGaugeFixing.h"
#include "Data/Lattice/CIndexData.h"
#include "Data/Lattice/CIndex.h"
#include "Data/Lattice/CIndexSquare.h"
#include "Data/Boundary/CBoundaryConditionTorusSquare.h"
#include "Core/Distributed/CLGComm.h"

__BEGIN_NAMESPACE

#if _CLG_MULTI_GPU
namespace
{
    //Recompute the block/thread decomposition constants for a GLOBAL lattice
    //(single rank), mirroring CLGLibManager's auto-decompose block
    //(CLGLibManager.cpp, "if (bAutoDecompose)") so kernels launched under the
    //temporary global context use exactly the grid/thread config of a single-GPU
    //run on the whole lattice. Call AFTER the global Lx/Ly/Lz/Lt are written and
    //BEFORE CopyConstants(). P4-2.2: without this the fixing kernels would launch
    //with the LOCAL decomposition and only sweep the local sub-lattice.
    void SetupGlobalDecompose(CCudaHelper* pHelper, const UINT* pG)
    {
        const TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
        const UINT uiLx = pG[0];
        const UINT uiLy = pG[1];
        const UINT uiLz = pG[2];
        const UINT uiLt = pG[3];
        const UINT uiVol = uiLx * uiLy * uiLz * uiLt;
        const UINT uiDir = static_cast<UINT>(_DC_Dir);

        //4D (x,y,z) decomposition.
        TArray<UINT> lat;
        lat.AddItem(uiLx);
        lat.AddItem(uiLy);
        lat.AddItem(uiLz);
        const TArray<UINT> d = _getDecompose(deviceConstraints, lat);
        pHelper->m_ConstIntegers[ECI_DecompX] = d[0];
        pHelper->m_ConstIntegers[ECI_DecompY] = d[1];
        pHelper->m_ConstIntegers[ECI_DecompZ] = d[2];
        pHelper->m_ConstIntegers[ECI_DecompLx] = d[3];
        pHelper->m_ConstIntegers[ECI_DecompLy] = d[4];
        pHelper->m_ConstIntegers[ECI_DecompLz] = d[5];
        pHelper->m_ConstIntegers[ECI_ThreadCountPerBlock] = d[3] * d[4] * d[5];

        //Flat 1D decomposition over the whole volume (+dir / +half variants).
        const UINT constrainx = min(deviceConstraints[0], deviceConstraints[1]);
        UINT uib = uiVol > constrainx ? appCeil(uiVol, constrainx) : 1;
        UINT uit = uiVol > constrainx ? appCeil(uiVol, uib) : uiVol;
        UINT uitd = uiDir * (constrainx / uiDir);
        UINT uibd = (uiVol * uiDir + uitd - 1) / uitd;
        pHelper->m_ConstIntegers[ECI_DecompAllBlock] = uib;
        pHelper->m_ConstIntegers[ECI_DecompAllThread] = uit;
        pHelper->m_ConstIntegers[ECI_DecompAllBlockDir] = uibd;
        pHelper->m_ConstIntegers[ECI_DecompAllThreadDir] = uitd;

        const UINT uiVolHalf = uiVol / 2;
        uib = uiVolHalf > constrainx ? appCeil(uiVolHalf, constrainx) : 1;
        uit = uiVolHalf > constrainx ? appCeil(uiVolHalf, uib) : uiVolHalf;
        uitd = uiDir * (constrainx / uiDir);
        uibd = (uiVolHalf * uiDir + uitd - 1) / uitd;
        pHelper->m_ConstIntegers[ECI_DecompAllBlockHalf] = uib;
        pHelper->m_ConstIntegers[ECI_DecompAllThreadHalf] = uit;
        pHelper->m_ConstIntegers[ECI_DecompAllBlockDirHalf] = uibd;
        pHelper->m_ConstIntegers[ECI_DecompAllThreadDirHalf] = uitd;

        //3D slice decompositions (x,y,z), (x,y,t), (x,z,t), (y,z,t).
        lat.RemoveAll();
        lat.AddItem(uiLx);
        lat.AddItem(uiLy);
        lat.AddItem(uiLz);
        TArray<UINT> d3 = _getDecompose(deviceConstraints, lat);
        pHelper->m_ConstIntegers[ECI_DecompX3D] = d3[0];
        pHelper->m_ConstIntegers[ECI_DecompY3D] = d3[1];
        pHelper->m_ConstIntegers[ECI_DecompZ3D] = d3[2];
        pHelper->m_ConstIntegers[ECI_DecompLx3D] = d3[3];
        pHelper->m_ConstIntegers[ECI_DecompLy3D] = d3[4];
        pHelper->m_ConstIntegers[ECI_DecompLz3D] = d3[5];

        lat.RemoveAll();
        lat.AddItem(uiLx);
        lat.AddItem(uiLy);
        lat.AddItem(uiLt);
        d3 = _getDecompose(deviceConstraints, lat);
        pHelper->m_ConstIntegers[ECI_DecompX3DXYT] = d3[0];
        pHelper->m_ConstIntegers[ECI_DecompY3DXYT] = d3[1];
        pHelper->m_ConstIntegers[ECI_DecompZ3DXYT] = d3[2];
        pHelper->m_ConstIntegers[ECI_DecompLx3DXYT] = d3[3];
        pHelper->m_ConstIntegers[ECI_DecompLy3DXYT] = d3[4];
        pHelper->m_ConstIntegers[ECI_DecompLz3DXYT] = d3[5];

        lat.RemoveAll();
        lat.AddItem(uiLx);
        lat.AddItem(uiLz);
        lat.AddItem(uiLt);
        d3 = _getDecompose(deviceConstraints, lat);
        pHelper->m_ConstIntegers[ECI_DecompX3DXZT] = d3[0];
        pHelper->m_ConstIntegers[ECI_DecompY3DXZT] = d3[1];
        pHelper->m_ConstIntegers[ECI_DecompZ3DXZT] = d3[2];
        pHelper->m_ConstIntegers[ECI_DecompLx3DXZT] = d3[3];
        pHelper->m_ConstIntegers[ECI_DecompLy3DXZT] = d3[4];
        pHelper->m_ConstIntegers[ECI_DecompLz3DXZT] = d3[5];

        lat.RemoveAll();
        lat.AddItem(uiLy);
        lat.AddItem(uiLz);
        lat.AddItem(uiLt);
        d3 = _getDecompose(deviceConstraints, lat);
        pHelper->m_ConstIntegers[ECI_DecompX3DYZT] = d3[0];
        pHelper->m_ConstIntegers[ECI_DecompY3DYZT] = d3[1];
        pHelper->m_ConstIntegers[ECI_DecompZ3DYZT] = d3[2];
        pHelper->m_ConstIntegers[ECI_DecompLx3DYZT] = d3[3];
        pHelper->m_ConstIntegers[ECI_DecompLy3DYZT] = d3[4];
        pHelper->m_ConstIntegers[ECI_DecompLz3DYZT] = d3[5];
    }
}
#endif

#if _CLG_MULTI_GPU
UBOOL CGaugeFixing::MGEnterGlobalFixerContext()
{
    if (NULL == appGetComm())
    {
        //Single GPU: the lattice is already the whole lattice. No-op.
        return TRUE;
    }
    if (!appGetComm()->IsRoot())
    {
        //Non-root ranks have nothing to prepare; they block on the scatter
        //call issued by the caller after the root finished fixing.
        return FALSE;
    }

    CCudaHelper* pHelper = appGetCudaHelper();

    //Save the local (decomposed) lattice constants.
    memcpy(m_uiSavedConstIntegers, pHelper->m_ConstIntegers, sizeof(UINT) * 128);

    //Switch to the GLOBAL lattice so the single-GPU fixing loop and CheckRes
    //run over the whole gathered lattice. Identity on single-GPU (no comm).
    const UINT* pG = appGetComm()->GlobalLattice();
    const UINT uiLx = pG[0];
    const UINT uiLy = pG[1];
    const UINT uiLz = pG[2];
    const UINT uiLt = pG[3];
    const UINT uiVol = uiLx * uiLy * uiLz * uiLt;
    const UINT uiDir = static_cast<UINT>(_DC_Dir);
    pHelper->m_ConstIntegers[ECI_Lx] = uiLx;
    pHelper->m_ConstIntegers[ECI_Ly] = uiLy;
    pHelper->m_ConstIntegers[ECI_Lz] = uiLz;
    pHelper->m_ConstIntegers[ECI_Lt] = uiLt;
    pHelper->m_ConstIntegers[ECI_Volume] = uiVol;
    pHelper->m_ConstIntegers[ECI_VolumeHalf] = uiVol / 2;
    pHelper->m_ConstIntegers[ECI_Volume_xyz] = uiLx * uiLy * uiLz;
    pHelper->m_ConstIntegers[ECI_Volume_xyt] = uiLx * uiLy * uiLt;
    pHelper->m_ConstIntegers[ECI_Volume_xzt] = uiLx * uiLz * uiLt;
    pHelper->m_ConstIntegers[ECI_Volume_yzt] = uiLy * uiLz * uiLt;
    pHelper->m_ConstIntegers[ECI_PlaqutteCount] = uiVol * uiDir * (uiDir - 1) / 2;
    pHelper->m_ConstIntegers[ECI_LinkCount] = uiVol * uiDir;
    pHelper->m_ConstIntegers[ECI_MultX] = uiLy * uiLz * uiLt;
    pHelper->m_ConstIntegers[ECI_MultY] = uiLz * uiLt;
    pHelper->m_ConstIntegers[ECI_MultZ] = uiLt;
    pHelper->m_ConstIntegers[ECI_GridDimZT] = uiLz * uiLt;

    //A global (single-rank) lattice context: no process-grid split, no offset.
    //The whole table is saved above and restored on exit, so the local
    //(decomposed) values come back unchanged.
    pHelper->m_ConstIntegers[ECI_GpuGridX] = 1;
    pHelper->m_ConstIntegers[ECI_GpuGridY] = 1;
    pHelper->m_ConstIntegers[ECI_GpuGridZ] = 1;
    pHelper->m_ConstIntegers[ECI_GpuGridT] = 1;
    pHelper->m_ConstIntegers[ECI_GlobalOffsetX] = 0;
    pHelper->m_ConstIntegers[ECI_GlobalOffsetY] = 0;
    pHelper->m_ConstIntegers[ECI_GlobalOffsetZ] = 0;
    pHelper->m_ConstIntegers[ECI_GlobalOffsetT] = 0;

    //Block/thread decomposition for the GLOBAL lattice (P4-2.2: fixing kernels
    //launch with these constants).
    SetupGlobalDecompose(pHelper, pG);
    pHelper->CopyConstants();

    //Build a GLOBAL index cache (alloc + bake under the global constants).
    m_pSavedGlobalIndexCache = new CIndexData();
    m_pSavedGlobalIndex = new CIndexSquare();
    m_pSavedGlobalIndex->SetBoundaryCondition(new CBoundaryConditionTorusSquare());
    m_pSavedGlobalIndex->BakeAllIndexBuffer(m_pSavedGlobalIndexCache);
    checkCudaErrors(cudaGetLastError());
    pHelper->SetDeviceIndex(m_pSavedGlobalIndexCache);
    return TRUE;
}

void CGaugeFixing::MGExitGlobalFixerContext()
{
    if (NULL == appGetComm() || !appGetComm()->IsRoot())
    {
        return;
    }

    CCudaHelper* pHelper = appGetCudaHelper();

    //Restore the local (decomposed) lattice constants.
    memcpy(pHelper->m_ConstIntegers, m_uiSavedConstIntegers, sizeof(UINT) * 128);
    pHelper->CopyConstants();

    //Restore the local index cache and free the temporary global one.
    if (NULL != appGetLattice() && NULL != appGetLattice()->m_pIndexCache)
    {
        pHelper->SetDeviceIndex(appGetLattice()->m_pIndexCache);
    }
    appSafeDelete(m_pSavedGlobalIndex);
    appSafeDelete(m_pSavedGlobalIndexCache);
    m_pSavedGlobalIndex = NULL;
    m_pSavedGlobalIndexCache = NULL;
}
#endif

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
