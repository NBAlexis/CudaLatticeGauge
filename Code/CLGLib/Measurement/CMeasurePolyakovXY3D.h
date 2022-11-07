//=============================================================================
// FILENAME : CMeasurePolyakovXY3D.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVXY3D_H_
#define _CMEASUREPOLYAKOVXY3D_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovXY3D)

class CLGAPI CMeasurePolyakovXY3D : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakovXY3D)

public:

    enum { _kGammaInInterests = 8, };

    CMeasurePolyakovXY3D()
        : CMeasure()
          , m_pXYHostLoopDensity(NULL)
          , m_pXYDeviceLoopDensity(NULL)
          , m_pXYHostLoopDensityAbs(NULL)
          , m_pXYDeviceLoopDensityAbs(NULL)

          , m_pDistributionR(NULL)
          , m_pDistributionP(NULL)
          , m_pDistributionPAbs(NULL)
          , m_pHostDistributionR(NULL)
          , m_pHostDistributionP(NULL)
          , m_pHostDistributionPAbs(NULL)

          , m_uiConfigurationCount(0)
          , m_uiMaxR(1)
          , m_uiEdgeR(1)
          , m_bMeasureDistribution(FALSE)
          , m_bShowResult(TRUE)
          , m_bU1(FALSE)
          , m_bShiftCenter(FALSE)
          , m_cAverageLoop()
    {
    }

    ~CMeasurePolyakovXY3D();

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
    CLGComplex* m_pXYDeviceLoopDensity;
    Real* m_pXYHostLoopDensityAbs;
    Real* m_pXYDeviceLoopDensityAbs;

    //The count of points with x^2+y^2=r^2
    UINT* m_pDistributionR;
    //<P>(R^2)
    CLGComplex* m_pDistributionP;
    Real* m_pDistributionPAbs;

    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistributionP;
    Real* m_pHostDistributionPAbs;

    UINT m_uiConfigurationCount;
    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

    UBOOL m_bU1;

public:

    UBOOL m_bShiftCenter;

    //all
    TArray<CLGComplex> m_lstLoop;
    //inner
    TArray<CLGComplex> m_lstLoopInner;

    TArray<Real> m_lstLoopAbs;
    TArray<Real> m_lstLoopAbsInner;

    //not using
    TArray<CLGComplex> m_lstLoopDensity;

    CLGComplex m_cAverageLoop;
    TArray<CLGComplex> m_lstAverageLoopDensity;

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstP;
    TArray<Real> m_lstPAbs;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVXY3D_H_

//=============================================================================
// END OF FILE
//=============================================================================