//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating3D.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [27/10/2022 nbale]
//=============================================================================
#include "CActionGaugePlaquetteRotating.h"

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING3D_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATING3D_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotating3D)

class CLGAPI CActionGaugePlaquetteRotating3D : public CActionGaugePlaquetteRotating
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotating3D)
public:

    CActionGaugePlaquetteRotating3D();
    
    CCString GetInfos(const CCString &tab) const override;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_

//=============================================================================
// END OF FILE
//=============================================================================