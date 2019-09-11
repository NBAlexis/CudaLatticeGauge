//=============================================================================
// FILENAME : CMeasurePolyakovXY.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop
//
// REVISION:
//  [05/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVXY_H_
#define _CMEASUREPOLYAKOVXY_H_

__BEGIN_NAMESPACE

/**
* array[x, y] = 0
*/
extern CLGAPI void _ZeroXYPlaneC(CLGComplex* pDeviceRes);

/**
* No need to initial pRes = 0
*/
extern CLGAPI void _PolyakovAtSite(const deviceSU3* __restrict__ pDeviceBuffer, deviceSU3* pRes);

__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovXY)

class CLGAPI CMeasurePolyakovXY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakovXY)

public:

    enum { _kGammaInInterests = 8, };

    CMeasurePolyakovXY()
        : CMeasure()
          , m_pXYHostLoopDensity(NULL)
          , m_pTmpDeviceSum(NULL)
          , m_pXYDeviceLoopDensity(NULL)
          , m_pTmpLoop(NULL)

          , m_pDistributionR(NULL)
          , m_pDistributionP(NULL)
          , m_pHostDistributionR(NULL)
          , m_pHostDistributionP(NULL)

          , m_uiConfigurationCount(0)
          , m_uiMaxR(1)
          , m_uiEdgeR(1)
          , m_bMeasureDistribution(FALSE)
          , m_bShowResult(TRUE)
          , m_cAverageLoop()
    {
    }

    ~CMeasurePolyakovXY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    CLGComplex* m_pXYHostLoopDensity;
    CLGComplex* m_pTmpDeviceSum;
    CLGComplex* m_pXYDeviceLoopDensity;
    deviceSU3* m_pTmpLoop;

    //The count of points with x^2+y^2=r^2
    UINT* m_pDistributionR;
    //<P>(R^2)
    CLGComplex* m_pDistributionP;

    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistributionP;

    UINT m_uiConfigurationCount;
    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    TArray<CLGComplex> m_lstLoop;
    TArray<CLGComplex> m_lstLoopInner;
    TArray<CLGComplex> m_lstLoopDensity;
    
    CLGComplex m_cAverageLoop;
    TArray<CLGComplex> m_lstAverageLoopDensity;

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstP;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVXY_H_

//=============================================================================
// END OF FILE
//=============================================================================