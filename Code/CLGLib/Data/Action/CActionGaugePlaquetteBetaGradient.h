//=============================================================================
// FILENAME : CActionGaugePlaquetteGradient.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
// It always use clover energy so don't need to set clover energy
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTEBETAGRADIENT_H_
#define _CACTIONGAUGEPLAQUETTEBETAGRADIENT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteGradient)

class CLGAPI CActionGaugePlaquetteGradient : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteGradient)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquetteGradient();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;
    void SetBetaList(const TArray<DOUBLE>& fBeta);

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    DOUBLE* m_pDeviceBetaArray;
    TArray<DOUBLE> m_fBetaArray;

    //Not using it
    //UBOOL m_bUsing4PlaqutteEnergy;
    UINT m_uiPlaqutteCount;

    DOUBLE CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const;

    void CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTEBETAGRADIENT_H_

//=============================================================================
// END OF FILE
//=============================================================================