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

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

#pragma endregion

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

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirOne(
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Id()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx].m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirZero(
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Zero()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx].m_uiSiteIndex, byDir)];
}

/**
 * calculate D_mu A _nu = Delta _mu + [A_mu, A _nu]
 * Use U now to calculate A pure
 * me will be changed, so, if me is A phys, copy me first
 */
static __device__ __inline__ deviceSU3 _deviceDPureMu(
    const deviceSU3* __restrict__ piA, 
    const deviceSU3* __restrict__ piApure,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->m_pWalkingTable[uiBigIdx * _DC_Dir * 2 + byMu];

    deviceSU3 res = _deviceGetGaugeBCSU3DirZero(piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCSU3DirZero(piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU3DirZero(piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZero(piA, uiBigIdx, byNu));
    res.Sub(_deviceGetGaugeBCSU3DirZero(piA, uiSiteBig_m_mu, byNu));
    return res;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3D_H_

//=============================================================================
// END OF FILE
//=============================================================================