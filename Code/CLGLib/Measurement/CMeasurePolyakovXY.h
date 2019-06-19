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
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasurePolyakovXY();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) {}
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return TRUE; }
    virtual UBOOL IsSourceScanning() const { return FALSE; }

protected:

    CLGComplex* m_pXYHostLoopDensity;
    CLGComplex* m_pTmpDeviceSum;
    CLGComplex* m_pXYDeviceLoopDensity;
    deviceSU3* m_pTmpLoop;

    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;

public:

    TArray<CLGComplex> m_lstLoop;
    TArray<CLGComplex> m_lstLoopDensity;

    CLGComplex m_cAverageLoop;
    TArray<CLGComplex> m_lstAverageLoopDensity;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOVXY_H_

//=============================================================================
// END OF FILE
//=============================================================================