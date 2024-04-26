//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingU1_3D.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGU1_3D_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATINGU1_3D_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingU1_3D)

class CLGAPI CActionGaugePlaquetteRotatingU1_3D : public CActionGaugePlaquetteRotatingU1
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingU1_3D)
public:

    CActionGaugePlaquetteRotatingU1_3D();

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;

    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    CCString GetInfos(const CCString &tab) const override;

};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGU1_3D_H_

//=============================================================================
// END OF FILE
//=============================================================================