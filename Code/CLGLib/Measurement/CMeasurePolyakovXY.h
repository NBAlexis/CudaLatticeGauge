//=============================================================================
// FILENAME : CMeasurePolyakovXY.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop
// 
// Now support all gauge fields
// and support all U_mu, and all slice
// and support r-distribution for Polyakov and Z loop
//
// REVISION:
//  [05/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVXY_H_
#define _CMEASUREPOLYAKOVXY_H_

__BEGIN_NAMESPACE


/**
* No need to initial pRes = 0
*/
//extern CLGAPI void _PolyakovAtSite(const deviceSU3* __restrict__ pDeviceBuffer, deviceSU3* pRes, BYTE byFieldId);

__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakovXY)

class CLGAPI CMeasurePolyakovXY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakovXY)

public:

    //enum { _kGammaInInterests = 8, };

    CMeasurePolyakovXY()
        : CMeasure()
          , m_pXYHostLoopDensity(NULL)
          , m_pXHostLoopDensity(NULL)
          , m_pYHostLoopDensity(NULL)
          , m_pZHostLoopDensity(NULL)
          , m_pTHostLoopDensity(NULL)
          , m_pXYDeviceLoopDensity(NULL)
          , m_pXDeviceLoopDensity(NULL)
          , m_pYDeviceLoopDensity(NULL)
          , m_pZDeviceLoopDensity(NULL)
          , m_pTDeviceLoopDensity(NULL)

          , m_pXYHostLoopDensityAbs(NULL)
          , m_pXHostLoopDensityAbs(NULL)
          , m_pYHostLoopDensityAbs(NULL)
          , m_pZHostLoopDensityAbs(NULL)
          , m_pTHostLoopDensityAbs(NULL)
          , m_pXYDeviceLoopDensityAbs(NULL)
          , m_pXDeviceLoopDensityAbs(NULL)
          , m_pYDeviceLoopDensityAbs(NULL)
          , m_pZDeviceLoopDensityAbs(NULL)
          , m_pTDeviceLoopDensityAbs(NULL)

          , m_pTmpLoop(NULL)
          , m_pTmpLoopX(NULL)
          , m_pTmpLoopY(NULL)
          , m_pTmpLoopZ(NULL)

          , m_pDistributionR(NULL)
          , m_pDistributionP(NULL)
          , m_pDistributionPAbs(NULL)
          , m_pHostDistributionR(NULL)
          , m_pHostDistributionP(NULL)
          , m_pHostDistributionPAbs(NULL)

          , m_uiMaxR(1)
          , m_uiEdgeR(1)
          , m_bMeasureDistribution(TRUE)
          , m_bMeasureAbs(FALSE)
          , m_bMeasureLoopX(FALSE)
          , m_bMeasureLoopY(FALSE)
          , m_bMeasureLoopZ(FALSE)
          , m_bMeasureXSlice(FALSE)
          , m_bMeasureYSlice(FALSE)
          , m_bMeasureZSlice(FALSE)
          , m_bMeasureTSlice(FALSE)
          , m_bShiftCenter(FALSE)
          , m_cAverageLoop()
    {
    }

    ~CMeasurePolyakovXY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    cuDoubleComplex* m_pXYHostLoopDensity;
    cuDoubleComplex* m_pXHostLoopDensity;
    cuDoubleComplex* m_pYHostLoopDensity;
    cuDoubleComplex* m_pZHostLoopDensity;
    cuDoubleComplex* m_pTHostLoopDensity;

    cuDoubleComplex* m_pXYDeviceLoopDensity;
    cuDoubleComplex* m_pXDeviceLoopDensity;
    cuDoubleComplex* m_pYDeviceLoopDensity;
    cuDoubleComplex* m_pZDeviceLoopDensity;
    cuDoubleComplex* m_pTDeviceLoopDensity;

    DOUBLE* m_pXYHostLoopDensityAbs;
    DOUBLE* m_pXHostLoopDensityAbs;
    DOUBLE* m_pYHostLoopDensityAbs;
    DOUBLE* m_pZHostLoopDensityAbs;
    DOUBLE* m_pTHostLoopDensityAbs;
    DOUBLE* m_pXYDeviceLoopDensityAbs;
    DOUBLE* m_pXDeviceLoopDensityAbs;
    DOUBLE* m_pYDeviceLoopDensityAbs;
    DOUBLE* m_pZDeviceLoopDensityAbs;
    DOUBLE* m_pTDeviceLoopDensityAbs;

    cuDoubleComplex* m_pTmpLoop;
    cuDoubleComplex* m_pTmpLoopX;
    cuDoubleComplex* m_pTmpLoopY;
    cuDoubleComplex* m_pTmpLoopZ;

    //The count of points with x^2+y^2=r^2
    UINT* m_pDistributionR;
    //<P>(R^2)
    cuDoubleComplex* m_pDistributionP;
    DOUBLE* m_pDistributionPAbs;

    UINT* m_pHostDistributionR;
    cuDoubleComplex* m_pHostDistributionP;
    DOUBLE* m_pHostDistributionPAbs;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bMeasureAbs;

public:

    UBOOL m_bMeasureLoopX;
    UBOOL m_bMeasureLoopY;
    UBOOL m_bMeasureLoopZ;
    UBOOL m_bMeasureXSlice;
    UBOOL m_bMeasureYSlice;
    UBOOL m_bMeasureZSlice;
    UBOOL m_bMeasureTSlice;

    //shift center is used to decide r
    UBOOL m_bShiftCenter;

    //all
    TArray<cuDoubleComplex> m_lstLoop;
    TArray<cuDoubleComplex> m_lstLoopX;
    TArray<cuDoubleComplex> m_lstLoopY;
    TArray<cuDoubleComplex> m_lstLoopZ;

    //inner
    TArray<cuDoubleComplex> m_lstLoopInner;
    TArray<cuDoubleComplex> m_lstLoopInnerZ;

    //all
    TArray<DOUBLE> m_lstLoopAbs;
    TArray<DOUBLE> m_lstLoopAbsZ;

    //inner
    TArray<DOUBLE> m_lstLoopAbsInner;
    TArray<DOUBLE> m_lstLoopAbsInnerZ;

    //not using, only for log
    CLGComplex m_cAverageLoop;
    TArray<UINT> m_lstR;

    //OVER R
    TArray<cuDoubleComplex> m_lstP;
    TArray<DOUBLE> m_lstPAbs;

    //OVER R
    TArray<cuDoubleComplex> m_lstPZ;
    TArray<DOUBLE> m_lstPZAbs;

    TArray<cuDoubleComplex> m_lstP_XSlice;
    TArray<cuDoubleComplex> m_lstP_YSlice;
    TArray<cuDoubleComplex> m_lstP_ZSlice;
    TArray<DOUBLE> m_lstP_XSliceAbs;
    TArray<DOUBLE> m_lstP_YSliceAbs;
    TArray<DOUBLE> m_lstP_ZSliceAbs;

    TArray<cuDoubleComplex> m_lstPX_YSlice;
    TArray<cuDoubleComplex> m_lstPX_ZSlice;
    TArray<cuDoubleComplex> m_lstPX_TSlice;
    TArray<DOUBLE> m_lstPX_YSliceAbs;
    TArray<DOUBLE> m_lstPX_ZSliceAbs;
    TArray<DOUBLE> m_lstPX_TSliceAbs;

    TArray<cuDoubleComplex> m_lstPY_XSlice;
    TArray<cuDoubleComplex> m_lstPY_ZSlice;
    TArray<cuDoubleComplex> m_lstPY_TSlice;
    TArray<DOUBLE> m_lstPY_XSliceAbs;
    TArray<DOUBLE> m_lstPY_ZSliceAbs;
    TArray<DOUBLE> m_lstPY_TSliceAbs;

    TArray<cuDoubleComplex> m_lstPZ_XSlice;
    TArray<cuDoubleComplex> m_lstPZ_YSlice;
    TArray<cuDoubleComplex> m_lstPZ_TSlice;
    TArray<DOUBLE> m_lstPZ_XSliceAbs;
    TArray<DOUBLE> m_lstPZ_YSliceAbs;
    TArray<DOUBLE> m_lstPZ_TSliceAbs;

    void Export(const CCString& sCSV, UINT uiStart, UINT uiEnd, UINT uiO, UINT uiOStart) const
    {
        Export(sCSV, uiStart, uiEnd, appToString(uiO), uiO, uiOStart);
    }

    void Export(const CCString& sCSV, UINT uiStart, UINT uiEnd, const CCString& sOName, UINT uiO, UINT uiOStart) const;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVXY_H_

//=============================================================================
// END OF FILE
//=============================================================================