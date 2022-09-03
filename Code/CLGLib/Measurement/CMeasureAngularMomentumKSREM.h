//=============================================================================
// FILENAME : CMeasureAngularMomentumKSREM.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================

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
        , m_bTwistedBoundary(TRUE)
    {
        
    }

protected:

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void ApplyOrbitalMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const override;
    void ApplySpinMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const override;

    //Potential angular momentum uses gamma4, which do not have magnetic vector
    //void ApplyPotentialMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const override;
    BYTE m_byGaugeType;
    UBOOL m_bTwistedBoundary;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKSREM_H_

//=============================================================================
// END OF FILE
//=============================================================================