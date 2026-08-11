//=============================================================================
// FILENAME : CActionGaugePlaquettePSU3WithBoundary.h
//
// DESCRIPTION:
// The fundamental-lift Wilson/Villain PSU(3) action with a Z3 2-form field B.
//
//   S(U,B) = -(Beta / 3) * sum_p Re[B_p Tr(U_p)]
//
// B is a CFieldTensor2Z3 (a Z3 tensor 2-form defined on every plaquette) that
// carries the SU(3) lift center transition function, the global boundary twist
// omega in H^2(T^4,Z3)=Z3^6, and (when AllowMonopole=1) Z3 monopole defects.
//
// Alternating update (Markov composite kernel):
//   fixed B -> HMC evolves the SU(3) links U (B kept fixed in the trajectory)
//   fixed U -> Z3 heatbath sweep evolves B (after every HMC attempt, accepted
//              or rejected)
//
// Modes:
//   AllowMonopole = 0 : flat ensemble, dB=0 enforced by link-star + global
//                       coclosed-sheet sweeps; sums over all H2 twist sectors.
//   AllowMonopole = 1 : unrestricted Villain ensemble, independent plaquette
//                       heatbath; Z3 monopoles allowed.
//
// REVISION:
//  [08/09/26]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTE_PSU3WITHBOUNDARY_H_
#define _CACTIONGAUGEPLAQUETTE_PSU3WITHBOUNDARY_H_

#include "CAction.h"
#include "Data/Field/Tensor2/CFieldTensor2Z3.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquettePSU3WithBoundary)

class CLGAPI CActionGaugePlaquettePSU3WithBoundary : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquettePSU3WithBoundary)
public:

    CActionGaugePlaquettePSU3WithBoundary();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    /**
    * Full Energy override: the base CAction::Energy only forwards gauge fields,
    * while this action needs the working tensor2 (B) of the integrator copy.
    */
    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields) override;

    /**
    * Called before the fields are copied back (bWillBeAccept=TRUE) or
    * discarded (FALSE). Sweeps B on the state that will actually be kept,
    * then records the post-sweep energy for the next trajectory.
    */
    void OnFinishTrajectory(UBOOL bWillBeAccept, INT gaugeNum, INT bosonNum, INT tensor2Num,
        CFieldGauge* const* gaugeFields, CFieldBoson* const* bosonFields,
        CFieldTensor2* const* tensor2Fields) override;

    /**
    * Second callback after the copy back: publish the post-sweep energy as the
    * new current energy (do NOT fall back to the pre-sweep m_fNewEnergy).
    */
    void OnFinishTrajectory(UBOOL bAccepted) override;

    CCString GetInfos(const CCString& tab) const override;

protected:

    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;

    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    /**
    * Full energy of (U, B). The kernel reads the working B pointer passed in.
    */
    DOUBLE CalculateEnergyNow(const CFieldGaugeSU3* pGauge, const CFieldTensor2Z3* pBoundary);

private:

    void SweepZ3(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary);
    void SweepZ3Unconstrained(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary);
    void SweepZ3LinkStar(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary, BYTE byAlpha, UBOOL bEven);
    void SweepZ3GlobalTwist(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary, BYTE byPlaqIndex);
    UINT CountMonopoles(const CFieldTensor2Z3* pBoundary) const;

    static INT FindGaugeIndexById(INT num, const CFieldGauge* const* gaugeFields, BYTE byFieldId);
    static INT FindTensor2IndexById(INT num, const CFieldTensor2* const* tensor2Fields, BYTE byFieldId);

    BYTE m_byBoundaryTensor2FieldId;
    UBOOL m_bAllowMonopole;
    UINT m_uiZ3SweepsPerTrajectory;
    UINT m_uiGlobalTwistSweepsPerTrajectory;
    UBOOL m_bCheckMonopole;

    // lattice B, fixed during an HMC trajectory (used by force)
    const CFieldTensor2Z3* m_pBoundaryField;

    DOUBLE m_fPostSweepEnergy;
    UBOOL m_bPostSweepEnergyValid;

    // sheet position for the six global twist sheets: component mu/nu of the
    // orientation pair (mu,nu) is taken from this array
    UINT m_aiTwistSheetPosition[4];
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_PSU3WITHBOUNDARY_H_

//=============================================================================
// END OF FILE
//=============================================================================
