//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [05/07/2019 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATING_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotating)

class CLGAPI CActionGaugePlaquetteRotating : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotating)

public:

    CActionGaugePlaquetteRotating();
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString& tab) const override;

    void SetGaugeOmega(DOUBLE fOmega);
    DOUBLE GetOmega() const { return m_fOmega; }

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    void EnergyDirichlet(const class CFieldGaugeSU3* pGauge);
    void EnergyProjectivePlane(const class CFieldGaugeSU3* pGauge);
    void EnergyTorus(const class CFieldGaugeSU3* pGauge);

    void CalculateForceOnGaugeDirichlet(const class CFieldGaugeSU3* pGauge, class CFieldGaugeSU3* pForce) const;
    void CalculateForceOnGaugeProjectivePlane(const class CFieldGaugeSU3* pGauge, class CFieldGaugeSU3* pForce) const;
    void CalculateForceOnGaugeTorus(const class CFieldGaugeSU3* pGauge, class CFieldGaugeSU3* pForce) const;


#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fOmega;
#else
    Real m_fOmega;
#endif
    UBOOL m_bCloverEnergy;
    UBOOL m_bShiftHalfCoord;
    UBOOL m_bTorus;

public:

    //===== test functions ======
    DOUBLE XYTerm1(const class CFieldGauge* pGauge);
    DOUBLE XYTerm2(const class CFieldGauge* pGauge);

protected:

    UINT m_uiPlaqutteCount;
};


__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_

//=============================================================================
// END OF FILE
//=============================================================================