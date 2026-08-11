//=============================================================================
// FILENAME : CActionGaugePlaquettePSU3.h
// 
// DESCRIPTION:
// PSU(N) adjoint plaquette action for SU2 and SU3.
// S_dyn = -beta_tilde * sum_P |Tr(U_P)|^2
// with beta_tilde = beta_A / (N^2 - 1) (SU(3): beta_A/8, SU(2): beta_A/3),
// where the YAML parameter Beta is beta_A.
//
// The force kernel returns only the raw matrix
//   Y = -beta_tilde * sum_P Tr(U_P) * staple,
// accumulated into the shared force field. The common integrator performs
// the U * Y^dag + TA projection, yielding
//   F_code = -beta_tilde/2 * sum_P [Tr(U_P^dag) U_P - Tr(U_P) U_P^dag].
//
// REVISION:
//  [07/02/26]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTE_PSU3_H_
#define _CACTIONGAUGEPLAQUETTE_PSU3_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquettePSU3)

class CLGAPI CActionGaugePlaquettePSU3 : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquettePSU3)
public:

    CActionGaugePlaquettePSU3();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    UINT m_uiPlaqutteCount;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_PSU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
