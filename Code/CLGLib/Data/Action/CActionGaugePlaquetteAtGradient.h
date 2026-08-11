//=============================================================================
// FILENAME : CActionGaugePlaquetteAtGradient.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
// It always use clover energy so don't need to set clover energy
//
// REVISION:
//  [mm/dd/yy]
//  [07/27/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONGAUGEPLAQUETTEATGRADIENT_H_
#define _CACTIONGAUGEPLAQUETTEATGRADIENT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteAtGradient)

class CLGAPI CActionGaugePlaquetteAtGradient : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteAtGradient)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquetteAtGradient();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;
    void SetXiList(const TArray<DOUBLE>& fXi);

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    DOUBLE* m_pDeviceXiArray;
    TArray<DOUBLE> m_fXiArray;

    //Not using it
    //UBOOL m_bUsing4PlaqutteEnergy;
    UINT m_uiPlaqutteCount;

    DOUBLE CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const;

    void CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTEATGRADIENT_H_

//=============================================================================
// END OF FILE
//=============================================================================