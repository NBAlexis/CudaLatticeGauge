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

extern CLGAPI void _AverageXYPlaneC(CLGComplex* pDeviceRes);

class CLGAPI CMeasureAMomentumJG : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumJG)
public:
    CMeasureAMomentumJG() 
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBufferOneConfig(NULL)
        , m_pHostSpinBuffer(NULL)
        , m_pDeviceSpinBuffer(NULL)
        , m_byFieldId(1)

        , m_pDistributionR(NULL)
        , m_pDistributionJG(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJG(NULL)
        , m_pDistributionJGS(NULL)
        , m_pHostDistributionJGS(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bMeasureDistribution(FALSE)

        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)

        , m_bMeasureSpin(FALSE)
        , m_pE(NULL)
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

    //those buffer are reused in calculation of JG chen
    Real* m_pHostSpinBuffer;
    Real* m_pDeviceSpinBuffer;
    BYTE m_byFieldId;

    UINT* m_pDistributionR;
    Real* m_pDistributionJG;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJG;

    //those buffer are reused in calculation of JG chen
    Real* m_pDistributionJGS;
    Real* m_pHostDistributionJGS;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;

    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    TArray<Real> m_lstRes;
    TArray<Real> m_lstResJGS;
    TArray<Real> m_lstResJGChen;

    UBOOL m_bMeasureSpin;
    CFieldGauge* m_pE;

public:

    TArray<Real> m_lstJGAll;
    TArray<Real> m_lstJGInner;

    TArray<UINT> m_lstR;
    //jg as function of R
    TArray<Real> m_lstJG;

    TArray<Real> m_lstJGSAll;
    TArray<Real> m_lstJGSInner;
    //jgs as function of R
    TArray<Real> m_lstJGS;

    TArray<Real> m_lstJGChenAll;
    TArray<Real> m_lstJGChenInner;
    //jg chen as function of R
    TArray<Real> m_lstJGChen;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMJG_H_

//=============================================================================
// END OF FILE
//=============================================================================