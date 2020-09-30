//=============================================================================
// FILENAME : CMeasure.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASURE_H_
#define _CMEASURE_H_

__BEGIN_NAMESPACE

class CLGAPI CMeasure : public CBase
{
public:
    CMeasure()
        : m_pOwner(NULL)
        , m_pLatticeData(NULL)
        , m_bNeedSmearing(FALSE)
        , m_byId(0)
        , m_byFieldId(0)
        , m_fLastRealResult(F(0.0))
        , m_cLastComplexResult()
    {
    }

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
    {
        m_pOwner = pOwner;
        m_pLatticeData = pLatticeData;
        m_byId = byId;

        INT iNeedGaugeSmearing = 0;
        param.FetchValueINT(_T("GaugeSmearing"), iNeedGaugeSmearing);
        m_bNeedSmearing = 0 != iNeedGaugeSmearing;

        INT iValue = 0;
        param.FetchValueINT(_T("FieldId"), iValue);
        m_byFieldId = static_cast<BYTE>(iValue);
    }

    /**
    * Accept gauge can be smoothed.
    * pCorrespondingStaple Might be NULL.
    */
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) = 0;


    /**
    * NOTE: sources will be passed to multiple measures, do NOT change the content!
    * NOTE: site.x start from 1 to Lx - 1, 0 is not included
    */
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) = 0;

    /**
    * Z4 Source
    */
    virtual void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd)
    {
        appCrucial(_T("OnConfigurationAcceptedZ4 not implemented"));
    }

    virtual void Average(UINT uiConfigurationCount) = 0;
    virtual void Report() = 0;
    virtual void Reset() = 0;

    virtual UBOOL IsGaugeMeasurement() const = 0;
    virtual UBOOL IsSourceScanning() const = 0;
    virtual UBOOL IsZ4Source() const { return FALSE; }
    virtual UBOOL NeedGaugeSmearing() const { return m_bNeedSmearing; }

    BYTE GetFieldId() const { return m_byFieldId; }

    static void LogGeneralComplex(const CLGComplex& cmp, UBOOL bHasComma = TRUE)
    {
        appGeneral(_T("%2.12f %s %2.12f I%s"),
            cmp.x,
            cmp.y < F(0.0) ? _T("-") : _T("+"),
            appAbs(cmp.y),
            bHasComma ? _T(",   ") : _T(" "));
    }

    static void FillDataWithR_R(
        TArray<Real>& arrData,
        TArray<Real>& arrInner,
        TArray<Real>& arrFull,
        TArray<UINT>& arrR,
        Real* hostData,
        UINT* hostR,
        UINT uiConfig,
        UINT uiMaxR,
        UINT uiEdgeR,
        UBOOL bFillR);

    static void FillDataWithR_C(
        TArray<CLGComplex>& arrData,
        TArray<CLGComplex>& arrInner,
        TArray<CLGComplex>& arrFull,
        TArray<UINT>& arrR,
        CLGComplex* hostData,
        UINT* hostR,
        UINT uiConfig,
        UINT uiMaxR,
        UINT uiEdgeR,
        UBOOL bFillR);

    
    static void ReportDistributionXY_R(UINT uiConfig, const TArray<Real>& arrayRes);
    static void ReportDistributionXY_C(UINT uiConfig, const TArray<CLGComplex>& arrayRes);

    static void _ZeroXYPlane(Real* pDeviceRes);
    static void _ZeroXYPlaneC(CLGComplex* pDeviceRes);
    /**
    * array[x, y] = array[x, y] / (lz * lt)
    */
    static void _AverageXYPlane(Real* pDeviceRes);
    static void _AverageXYPlaneC(CLGComplex* pDeviceRes);

    static void XYDataToRdistri_R(
        const Real* __restrict__ source,
        UINT* count,
        Real* result,
        UINT uiMaxR,
        UBOOL bCalculateCounter,
        BYTE byFieldId);

    static void XYDataToRdistri_C(
        const CLGComplex* __restrict__ source,
        UINT* count,
        CLGComplex* result,
        UINT uiMaxR,
        UBOOL bCalculateCounter,
        BYTE byFieldId);

    static void ReportDistributeWithR_R(UINT uiConf, UINT uiR, const TArray<Real>& arrayData);

protected:

    class CMeasurementManager* m_pOwner;
    class CLatticeData* m_pLatticeData;
    UBOOL m_bNeedSmearing;
    BYTE m_byId;
    BYTE m_byFieldId;

public:
    //============================================================
    //some simple measurement only produce real or complex results
    Real m_fLastRealResult;
    CLGComplex m_cLastComplexResult;
    
};

class CLGAPI CMeasureStochastic : public CMeasure
{
public:
    CMeasureStochastic()
        : CMeasure()
        , m_uiFieldCount(25)
    {
    }

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId) override
    {
        CMeasure::Initial(pOwner, pLatticeData, param, byId);
        INT iValue = 25;
        param.FetchValueINT(_T("FieldCount"), iValue);
        m_uiFieldCount = static_cast<UINT>(iValue);
    }

    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override
    {
        appCrucial(_T("Should not use SourceSanning"));
    }

    /**
    * Z4 Source
    */
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override
    {
        appCrucial(_T("OnConfigurationAcceptedZ4 not implemented"));
    }

    UBOOL IsGaugeMeasurement() const override = 0;
    UBOOL IsSourceScanning() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }
    UINT GetFieldCount() const { return m_uiFieldCount; }
    void SetFieldCount(UINT uiFieldCount) { m_uiFieldCount = uiFieldCount; }

protected:

    UINT m_uiFieldCount;

};

__END_NAMESPACE

#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================