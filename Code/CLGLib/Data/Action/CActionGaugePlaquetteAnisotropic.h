//=============================================================================
// FILENAME : CActionGaugePlaquetteAnisotropic.h
//
// DESCRIPTION:
// Anisotropic tree-level Symanzik gauge action.
// S = sum_{planes} w(mu,nu) * beta * [plaq + RectOverPlaq * rect],
// with w = 1/Xi for spatial planes and w = Xi for planes containing the
// temporal direction. Xi is the bare gauge anisotropy, not the measured one.
// The field level anisotropy kernels use the inverse convention (spatial
// planes are multiplied by xi, temporal planes by 1/xi), so 1/Xi is passed.
//
// REVISION:
//  [mm/dd/yy]
//  [07/25/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTEANISOTROPIC_H_
#define _CACTIONGAUGEPLAQUETTEANISOTROPIC_H_

#include "CActionGaugePlaquette.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteAnisotropic)

class CLGAPI CActionGaugePlaquetteAnisotropic : public CActionGaugePlaquette
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteAnisotropic)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquetteAnisotropic();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;

    //the bare gauge anisotropy
    DOUBLE m_fXi;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const class CFieldGauge* pGauge, UINT uiUpdateIterate) override;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTEANISOTROPIC_H_

//=============================================================================
// END OF FILE
//=============================================================================
