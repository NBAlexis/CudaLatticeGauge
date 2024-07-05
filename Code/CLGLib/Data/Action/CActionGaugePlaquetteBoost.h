//=============================================================================
// FILENAME : CActionGaugePlaquetteBoost.h
// 
// DESCRIPTION:
// 
// Periodic boundary is assumed
//
// REVISION:
//  [08/03/2020 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_BOOST_H_
#define _CACTIONGAUGEPLAQUETTE_BOOST_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteBoost)

class CLGAPI CActionGaugePlaquetteBoost : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteBoost)
public:

    CActionGaugePlaquetteBoost();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;
    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    UINT m_uiPlaqutteCount;
};


__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_

//=============================================================================
// END OF FILE
//=============================================================================