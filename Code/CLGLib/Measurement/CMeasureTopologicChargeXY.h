//=============================================================================
// FILENAME : CMeasureTopologicChargeXY.h
// 
// DESCRIPTION:
// NOTE: This is not ready yet!! I think the definition of topological charge maybe different in rotating frame
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURETOPOLOGICALCHARGEXY_H_
#define _CMEASURETOPOLOGICALCHARGEXY_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureTopologicChargeXY)

class CLGAPI CMeasureTopologicChargeXY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureTopologicChargeXY)
public:
    CMeasureTopologicChargeXY()
        : CMeasure()
        , m_pXYHostDensity(NULL)
        , m_pXYDeviceDensity(NULL)
    {
    }

    ~CMeasureTopologicChargeXY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    Real* m_pXYHostDensity;
    Real* m_pXYDeviceDensity;
    
    TArray<Real> m_lstXYDensity;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURETOPOLOGICALCHARGEXY_H_

//=============================================================================
// END OF FILE
//=============================================================================