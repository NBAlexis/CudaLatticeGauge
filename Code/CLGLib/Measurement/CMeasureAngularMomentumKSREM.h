//=============================================================================
// FILENAME : CMeasureAngularMomentumKSREM.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================
#include "CMeasureAngularMomentumKS.h"

#ifndef _CMEASUREANGULARMOMENTUMKSREM_H_
#define _CMEASUREANGULARMOMENTUMKSREM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAngularMomentumKSREM)

class CLGAPI CMeasureAngularMomentumKSREM : public CMeasureAngularMomentumKS
{
    __CLGDECLARE_CLASS(CMeasureAngularMomentumKSREM)
public:

    CMeasureAngularMomentumKSREM()
        : CMeasureAngularMomentumKS()
    {
        
    }

    void OnConfigurationAcceptedZ4SingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;

protected:

    //void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void ApplyOrbitalMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge, BYTE byGaugeFieldId) const override;
    void ApplySpinMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge, BYTE byGaugeFieldId) const override;

    void ApplyPotentialMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge, BYTE byGaugeFieldId) const override;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKSREM_H_

//=============================================================================
// END OF FILE
//=============================================================================