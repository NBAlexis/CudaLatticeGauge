//=============================================================================
// FILENAME : CMeasureAMomemtumJG.h
// 
// DESCRIPTION:
// This is measurement for angular momentum of rotating frame
// The angular momentum of each site is calculated, average is taken over z and t directions
// Then the average is taken over all configurations (result in angular momentum on a x-y plane)
//
// We assume lx * ly < max-thread
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#ifndef _CMEASUREAMOMENTUMJG_H_
#define _CMEASUREAMOMENTUMJG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAMomentumJG)

//=================== make them static functions because they are used also in CMeasureAMomentumJF ===================
/**
* array[x, y] = array[x, y] / (lz * lt)
*/
extern CLGAPI void _AverageXYPlane(Real* pDeviceRes);

/**
* array[x, y] = 0
*/
extern CLGAPI void _ZeroXYPlane(Real* pDeviceRes);

class CLGAPI CMeasureAMomentumJG : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumJG)
public:
    CMeasureAMomentumJG() 
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBufferOneConfig(NULL)
        , m_byFieldId(1)

        , m_pDistributionR(NULL)
        , m_pDistributionJG(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJG(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bMeasureDistribution(FALSE)

        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasureAMomentumJG();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    Real * m_pHostDataBuffer;
    Real * m_pDeviceDataBufferOneConfig;
    BYTE m_byFieldId;

    UINT* m_pDistributionR;
    Real* m_pDistributionJG;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJG;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;

    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    TArray<Real> m_lstRes;

public:

    TArray<Real> m_lstJGAll;
    TArray<Real> m_lstJGInner;

    TArray<UINT> m_lstR;
    TArray<Real> m_lstJG;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMJG_H_

//=============================================================================
// END OF FILE
//=============================================================================