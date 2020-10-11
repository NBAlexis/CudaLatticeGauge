//=============================================================================
// FILENAME : CActionGaugePlaquette.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_H_
#define _CACTIONGAUGEPLAQUETTE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquette)

class CLGAPI CActionGaugePlaquette : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquette)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquette();

#if !_CLG_DOUBLEFLOAT
    DOUBLE Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#else
    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#endif
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

#if !_CLG_DOUBLEFLOAT
    void SetBeta(DOUBLE fBeta);
#else
    void SetBeta(Real fBeta);
#endif
    //Real GetEnergyPerPlaqutte() const;
    UBOOL m_bCloverEnergy;

protected:

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fLastEnergy;
    DOUBLE m_fNewEnergy;
    DOUBLE m_fBetaOverN;
#else
    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
#endif

    //Not using it
    //UBOOL m_bUsing4PlaqutteEnergy;
    UINT m_uiPlaqutteCount;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_H_

//=============================================================================
// END OF FILE
//=============================================================================