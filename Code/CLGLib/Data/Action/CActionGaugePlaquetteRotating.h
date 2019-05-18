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

    virtual Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable);
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId);

    virtual UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const;
    virtual void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate);
    virtual void OnFinishTrajectory(UBOOL bAccepted);
    virtual CCString GetInfos(const CCString &tab) const;

    void SetBeta(Real fBeta);
    void SetOmega(Real fOmega) { m_fOmega = fOmega; }
    //Real GetEnergyPerPlaqutte() const;

    Real m_fOmega;
    SSmallInt4 m_sCenter;

protected:

    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_

//=============================================================================
// END OF FILE
//=============================================================================