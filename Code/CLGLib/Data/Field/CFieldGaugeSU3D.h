//=============================================================================
// FILENAME : CFieldGaugeSU3D.h
// 
// DESCRIPTION:
// There is simplifications for periodic boundary condition
// which is invalid for Dirichlet.
// This is only for Dirichlet.
//
// REVISION:
//  [05/17/2019 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU3D_H_
#define _CFIELDGAUGE_SU3D_H_

#define gaugeSU3KernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeSU3KernelFuncionEnd \
    } 



__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3D)

class CLGAPI CFieldGaugeSU3D : public CFieldGaugeSU3
{
    __CLGDECLARE_FIELD(CFieldGaugeSU3D)

public:

    CFieldGaugeSU3D() : CFieldGaugeSU3() {}

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;

    void MakeRandomGenerator() override;
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;

    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }

    Real CalculateKinematicEnergy() const override;

#pragma endregion

#pragma region BLAS

    void FixBoundary() override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;
    CCString GetInfos(const CCString &tab) const override;
};

#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3(
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[1])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3Dir(
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[1])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ]
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3D_H_

//=============================================================================
// END OF FILE
//=============================================================================