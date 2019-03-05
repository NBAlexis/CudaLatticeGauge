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

    virtual Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable);
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId);

    virtual UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple = NULL) const;
    virtual void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate);
    virtual void OnFinishTrajectory(UBOOL bAccepted);
    virtual CCString GetInfos(const CCString &tab) const;

    void SetBeta(Real fBeta);
    Real GetEnergyPerPlaqutte() const;

protected:

    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_H_

//=============================================================================
// END OF FILE
//=============================================================================