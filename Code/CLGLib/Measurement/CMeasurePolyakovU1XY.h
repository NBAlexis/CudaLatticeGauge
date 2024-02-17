//=============================================================================
// FILENAME : CMeasurePolyakovU1XY.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVU1XY_H_
#define _CMEASUREPOLYAKOVU1XY_H_

__BEGIN_NAMESPACE


/**
* No need to initial pRes = 0
*/
extern CLGAPI void _PolyakovAtSiteU1(const CLGComplex* __restrict__ pDeviceBuffer, CLGComplex* pRes);

__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovU1XY)

class CLGAPI CMeasurePolyakovU1XY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakovU1XY)

public:

    enum { _kGammaInInterests = 8, };

    CMeasurePolyakovU1XY()
        : CMeasure()
          , m_pXYHostLoopDensity(NULL)
          , m_pZHostLoopDensity(NULL)
          , m_pTmpDeviceSum(NULL)
          , m_pXYDeviceLoopDensity(NULL)
          , m_pZDeviceLoopDensity(NULL)
          , m_pTmpLoop(NULL)

          , m_pTmpLoopZ(NULL)

          , m_pXYHostLoopDensityAbs(NULL)
          , m_pZHostLoopDensityAbs(NULL)
          , m_pXYDeviceLoopDensityAbs(NULL)
          , m_pZDeviceLoopDensityAbs(NULL)

          , m_pDistributionR(NULL)
          , m_pDistributionP(NULL)
          , m_pDistributionPAbs(NULL)
          , m_pHostDistributionR(NULL)
          , m_pHostDistributionP(NULL)
          , m_pHostDistributionPAbs(NULL)

          , m_uiMaxR(1)
          , m_uiEdgeR(1)
          , m_bMeasureDistribution(FALSE)
          , m_bShowResult(TRUE)
          , m_bMeasureLoopZ(FALSE)
          , m_bMeasureZSlice(FALSE)
          , m_bShiftCenter(FALSE)
          , m_cAverageLoop()
    {
    }

    ~CMeasurePolyakovU1XY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    CLGComplex* m_pXYHostLoopDensity;
    CLGComplex* m_pZHostLoopDensity;
    CLGComplex* m_pTmpDeviceSum;
    CLGComplex* m_pXYDeviceLoopDensity;
    CLGComplex* m_pZDeviceLoopDensity;
    CLGComplex* m_pTmpLoop;
    CLGComplex* m_pTmpLoopZ;

    Real* m_pXYHostLoopDensityAbs;
    Real* m_pZHostLoopDensityAbs;
    Real* m_pXYDeviceLoopDensityAbs;
    Real* m_pZDeviceLoopDensityAbs;

    //The count of points with x^2+y^2=r^2
    UINT* m_pDistributionR;
    //<P>(R^2)
    CLGComplex* m_pDistributionP;
    Real* m_pDistributionPAbs;

    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistributionP;
    Real* m_pHostDistributionPAbs;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    UBOOL m_bMeasureLoopZ;
    UBOOL m_bMeasureZSlice;
    UBOOL m_bShiftCenter;

    //all
    TArray<CLGComplex> m_lstLoop;
    //inner
    TArray<CLGComplex> m_lstLoopInner;

    TArray<Real> m_lstLoopAbs;
    TArray<Real> m_lstLoopAbsInner;

    //not using
    //TArray<CLGComplex> m_lstLoopDensity;

    //all
    TArray<CLGComplex> m_lstLoopZ;
    //inner
    TArray<CLGComplex> m_lstLoopZInner;

    //not using
    TArray<CLGComplex> m_lstLoopZDensity;

    CLGComplex m_cAverageLoop;
    //TArray<CLGComplex> m_lstAverageLoopDensity;

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstP;
    TArray<Real> m_lstPAbs;
    TArray<CLGComplex> m_lstPZ;
    TArray<CLGComplex> m_lstPZSlice;
    TArray<Real> m_lstPZSliceAbs;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVU1XY_H_

//=============================================================================
// END OF FILE
//=============================================================================