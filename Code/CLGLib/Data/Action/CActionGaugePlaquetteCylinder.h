//=============================================================================
// FILENAME : CActionGaugePlaquetteCylinder.h
//
// DESCRIPTION:
// This is the gauge action in cylindrical coordinates.
// The lattice directions 0,1,2,3 are interpreted as (r, phi, z, t).
// The coupling is a position dependent function K(r) = beta(r) * w(r),
// where w = 1/r for planes containing phi, and w = r for other planes.
// The coupling of one plaquette is the average of its 4 corners.
// See Docs/Applications/Cylinder.tex for details.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTECYLINDER_H_
#define _CACTIONGAUGEPLAQUETTECYLINDER_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteCylinder)

class CLGAPI CActionGaugePlaquetteCylinder : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteCylinder)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquetteCylinder();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    DOUBLE* m_pDeviceBetaArray;
    TArray<DOUBLE> m_fBetaArray;

    //r = m_fRStart + m_fDeltaR * n0, r = 0 is a singularity, so m_fRStart must be positive
    Real m_fRStart;
    Real m_fREnd;
    Real m_fDeltaR;

    //same convention as CActionGaugePlaquette, whether the energy is calculated
    //by clover (per site) or by plaquette with explicit 4-corner coupling
    UBOOL m_bCloverEnergy;

    UINT m_uiPlaqutteCount;

    DOUBLE CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const;

    DOUBLE CalculatePlaqutteEnergyUsePlaqutte(const CFieldGaugeSU3* pGauge) const;

    void CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTECYLINDER_H_

//=============================================================================
// END OF FILE
//=============================================================================
