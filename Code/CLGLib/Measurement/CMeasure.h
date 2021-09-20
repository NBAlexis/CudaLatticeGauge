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

#if !_CLG_DOUBLEFLOAT
    static void LogGeneralComplex(const cuDoubleComplex& cmp, UBOOL bHasComma = TRUE)
    {
        appGeneral(_T("%2.12f %s %2.12f I%s"),
            cmp.x,
            cmp.y < F(0.0) ? _T("-") : _T("+"),
            appAbs(cmp.y),
            bHasComma ? _T(",   ") : _T(" "));
    }
#endif

    static void LogGeneralComplex(const CLGComplex& cmp, UBOOL bHasComma = TRUE)
    {
        appGeneral(_T("%2.12f %s %2.12f I%s"),
            cmp.x,
            cmp.y < F(0.0) ? _T("-") : _T("+"),
            appAbs(cmp.y),
            bHasComma ? _T(",   ") : _T(" "));
    }

#pragma region Distribution Common functions

    static void FillDataWithR_R(
        TArray<Real>& arrData,
        TArray<Real>* arrInner,
        TArray<Real>& arrFull,
        TArray<UINT>& arrR,
        Real* hostData,
        UINT* hostR,
        UINT uiConfig,
        UINT uiMaxR,
        UINT uiEdgeR,
        Real fDivider,
        UBOOL bFillR);

    static void FillDataWithR_C(
        TArray<CLGComplex>& arrData,
        TArray<CLGComplex>* arrInner,
        TArray<CLGComplex>& arrFull,
        TArray<UINT>& arrR,
        CLGComplex* hostData,
        UINT* hostR,
        UINT uiConfig,
        UINT uiMaxR,
        UINT uiEdgeR,
        Real fDivider,
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
        UBOOL bShiftCenter,
        const Real* __restrict__ source,
        UINT* count,
        Real* result,
        UINT uiMaxR,
        UBOOL bCalculateCounter,
        BYTE byFieldId);

    static void XYDataToRdistri_C(
        UBOOL bShiftCenter,
        const CLGComplex* __restrict__ source,
        UINT* count,
        CLGComplex* result,
        UINT uiMaxR,
        UBOOL bCalculateCounter,
        BYTE byFieldId);

    static void ReportDistributeWithR_R(UINT uiConf, UINT uiR, const TArray<Real>& arrayData);

    /**
     * TransformFromXYDataToRDataOnce_C and TransformFromXYDataToRDataOnce_R
     * is for gauge measurement
     */
    static void TransformFromXYDataToRDataOnce_C(
        UBOOL bShiftCenter,
        const CLGComplex* __restrict__ pXYData,
        UINT* pCountBuffer,
        CLGComplex* pValueBuffer,
        UINT* pHostCountBuffer,
        CLGComplex* pHostValueBuffer,
        UINT uiMaxR,
        UINT uiEdgeR,
        UBOOL bCalculateCounter,
        BYTE byFieldId,
        TArray<CLGComplex>& arrData,
        TArray<CLGComplex>* arrInner,
        TArray<CLGComplex>& arrFull,
        TArray<UINT>& arrR,
        UINT uiConfig,
        Real fDivider)
    {
        XYDataToRdistri_C(
            bShiftCenter,
            pXYData,
            pCountBuffer,
            pValueBuffer,
            uiMaxR,
            bCalculateCounter,
            byFieldId);

        if (bCalculateCounter)
        {
            checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
        }
        
        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(CLGComplex) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        //Here we have already divide by all XYZ points
        FillDataWithR_C(
            arrData,
            arrInner,
            arrFull,
            arrR,
            pHostValueBuffer,
            pHostCountBuffer,
            uiConfig,
            uiMaxR,
            uiEdgeR,
            fDivider,
            bCalculateCounter
        );
    }

    static void TransformFromXYDataToRDataOnce_R(
        UBOOL bShiftCenter,
        const Real* __restrict__ pXYData,
        UINT* pCountBuffer,
        Real* pValueBuffer,
        UINT* pHostCountBuffer,
        Real* pHostValueBuffer,
        UINT uiMaxR,
        UINT uiEdgeR,
        UBOOL bCalculateCounter,
        BYTE byFieldId,
        TArray<Real>& arrData,
        TArray<Real>* arrInner,
        TArray<Real>& arrFull,
        TArray<UINT>& arrR,
        UINT uiConfig,
        Real fDivider)
    {
        XYDataToRdistri_R(
            bShiftCenter,
            pXYData,
            pCountBuffer,
            pValueBuffer,
            uiMaxR,
            bCalculateCounter,
            byFieldId);

        if (bCalculateCounter)
        {
            checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
        }
        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(Real) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        //Here we have already divide by all XYZ points
        FillDataWithR_R(
            arrData,
            arrInner,
            arrFull,
            arrR,
            pHostValueBuffer,
            pHostCountBuffer,
            uiConfig,
            uiMaxR,
            uiEdgeR,
            fDivider,
            bCalculateCounter
        );
    }

    /**
     * TransformFromXYDataToRData_C and TransformFromXYDataToRData_R
     * is for Stochastic measurements
     */
    static void TransformFromXYDataToRData_C(
        UBOOL bShiftCenter,
        UINT uiMaxR,
        UINT uiEdgeR,
        BYTE byFieldId,
        UINT uiFieldCount,
        UINT uiMeasureCount,
        UINT uiConfig,
        const CLGComplex* const* pXYBuffers,
        UINT* pCountBuffer,
        CLGComplex* pValueBuffer,
        UINT* pHostCountBuffer,
        CLGComplex* pHostValueBuffer,
        TArray<UINT>& lstR,
        TArray<CLGComplex>* lstValues,
        TArray<CLGComplex>* lstAll,
        TArray<CLGComplex>* lstInner);

    static void TransformFromXYDataToRData_R(
        UBOOL bShiftCenter,
        UINT uiMaxR,
        UINT uiEdgeR,
        BYTE byFieldId,
        UINT uiFieldCount,
        UINT uiMeasureCount,
        UINT uiConfig,
        const Real* const* pXYBuffers,
        UINT* pCountBuffer,
        Real* pValueBuffer,
        UINT* pHostCountBuffer,
        Real* pHostValueBuffer,
        TArray<UINT>& lstR,
        TArray<Real>* lstValues,
        TArray<Real>* lstAll,
        TArray<Real>* lstInner);

    /**
     * Many measurements measure the XY distributions, needs to calculate max R and edge
     */
    static void SetMaxAndEdge(UINT* maxXY, UINT* edgeXY, UBOOL bShiftCenter)
    {
        if (bShiftCenter)
        {
            if (NULL != maxXY)
            {
                *maxXY = (_HC_Lx - 1) * (_HC_Lx - 1) + (_HC_Ly - 1) * (_HC_Ly - 1);
            }
            
            if (NULL != edgeXY)
            {
                *edgeXY = (_HC_Lx - 1) * (_HC_Lx - 1);
            }
        }
        else
        {
            if (NULL != maxXY)
            {
                *maxXY = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
                       + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2);
            }

            if (NULL != edgeXY)
            {
                *edgeXY = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1);
            }
        }
    }

#pragma endregion

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
        , m_bDebugDivation(FALSE)
    {
    }

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId) override
    {
        CMeasure::Initial(pOwner, pLatticeData, param, byId);
        INT iValue = 25;
        param.FetchValueINT(_T("FieldCount"), iValue);
        m_uiFieldCount = static_cast<UINT>(iValue);

        iValue = 0;
        param.FetchValueINT(_T("DebugDivation"), iValue);
        m_bDebugDivation = 0 != iValue;
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
    UBOOL m_bDebugDivation;

};

__END_NAMESPACE

#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================