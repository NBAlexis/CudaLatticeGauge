//=============================================================================
// FILENAME : CMeasure.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [mm/dd/yy]
//  [01/29/2019 nbale]
//=============================================================================
#include "Tools/Math/DeviceInlineTemplate.h"

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
        , m_byFermionFieldId(0)
        , m_bShowResult(FALSE)
        , m_uiConfigurationCount(0)
        , m_fAverageRealRes(F(0.0))
    {
        m_cAverageCmpRes = _zeroc;
    }

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId);


    /**
    * Accept gauge can be smoothed.
    * pCorrespondingStaple Might be NULL.
    */
    virtual void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple);


    /**
    * NOTE: sources will be passed to multiple measures, do NOT change the content!
    * NOTE: site.x start from 1 to Lx - 1, 0 is not included
    */
    virtual void SourceSanning(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site);

    /**
    * Z4 Source
    */
    virtual void OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd);

protected:


    /**
    * For single field case
    */
    virtual void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
    {
        appCrucial(_T("Single field case OnConfigurationAccepted not implemented!\n"));
    }

    /**
    * For single field case
    */
    virtual void SourceSanningSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site)
    {
        appCrucial(_T("Single field case SourceSanning not implemented!\n"));
    }

    /**
    * For single field case
    */
    virtual void OnConfigurationAcceptedZ4SingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd)
    {
        appCrucial(_T("Single field case OnConfigurationAcceptedZ4 not implemented!\n"));
    }

public:

    //UINT GetDefaultMatrixN() const;

    virtual void Average();
    virtual void Report() = 0;
    virtual void Reset()
    {
        m_uiConfigurationCount = 0;
        m_lstRealResults.RemoveAll();
        m_lstComplexResults.RemoveAll();
    }

    virtual UBOOL IsGaugeOrBosonMeasurement() const = 0;
    virtual UBOOL IsSourceScanning() const = 0;
    virtual UBOOL IsZ4Source() const { return FALSE; }
    virtual UBOOL NeedGaugeSmearing() const { return m_bNeedSmearing; }

    BYTE GetId() const { return m_byId; }
    BYTE GetFermionFieldId() const { return m_byFermionFieldId; }
    BYTE GetGaugeFieldIdSingleField() const { return m_lstGaugeFieldIds[0]; }

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

    template<class T>
    static void FillDataWithR(
        TArray<T>& arrData,
        TArray<T>* arrInner,
        TArray<T>& arrFull,
        TArray<UINT>& arrR,
        T* hostData,
        UINT* hostR,
        UINT uiConfig,
        UINT uiMaxR,
        UINT uiEdgeR,
        Real fDivider,
        UBOOL bFillR)
    {
        T fAverageJGInner = _makeZeroHost<T>();
        T fAverageJGAll = _makeZeroHost<T>();
        UINT uiInnerPointsAll = 0;
        UINT uiInnerPointsInner = 0;

        if (0 == uiConfig)
        {
            appAssert(!bFillR || 0 == arrR.Num());
            appAssert(0 == arrData.Num());

            for (UINT uiL = 0; uiL <= uiMaxR; ++uiL)
            {
                if (hostR[uiL] > 0)
                {
                    if (bFillR)
                    {
                        arrR.AddItem(uiL);
                    }

                    arrData.AddItem(_mulCHost(hostData[uiL], fDivider));

                    uiInnerPointsAll += hostR[uiL];
                    _addHost(fAverageJGAll, _mulCHost(hostData[uiL], fDivider * hostR[uiL]));
                    if (NULL != arrInner && uiL < uiEdgeR)
                    {
                        uiInnerPointsInner += hostR[uiL];
                        _addHost(fAverageJGInner, _mulCHost(hostData[uiL], fDivider * hostR[uiL]));
                    }
                }
            }
        }
        else
        {
            for (INT i = 0; i < arrR.Num(); ++i)
            {
                appAssert(hostR[arrR[i]] > 0);
                arrData.AddItem(_mulCHost(hostData[arrR[i]], fDivider));

                uiInnerPointsAll += hostR[arrR[i]];
                _addHost(fAverageJGAll, _mulCHost(hostData[arrR[i]], fDivider * hostR[arrR[i]]));
                if (NULL != arrInner && arrR[i] < uiEdgeR)
                {
                    uiInnerPointsInner += hostR[arrR[i]];
                    _addHost(fAverageJGInner, _mulCHost(hostData[arrR[i]], fDivider * hostR[arrR[i]]));
                }
            }
        }

        if (uiInnerPointsAll > 0)
        {
            _divHost(fAverageJGAll, uiInnerPointsAll);
        }
        if (NULL != arrInner && uiInnerPointsInner > 0)
        {
            _divHost(fAverageJGInner, uiInnerPointsInner);
        }
        arrFull.AddItem(fAverageJGAll);
        if (NULL != arrInner)
        {
            arrInner->AddItem(fAverageJGInner);
        }
    }
    
    template<class T>
    static void ReportDistributionXY(UINT uiConfig, const TArray<T>& arrayRes)
    {
        appAssert(uiConfig * (_HC_Lx - 1) * (_HC_Ly - 1)
            == static_cast<UINT>(arrayRes.Num()));

        TArray<T> tmpjgs;
        appGeneral(_T("{\n"));
        for (UINT k = 0; k < uiConfig; ++k)
        {
            appGeneral(_T("{"));
            for (UINT i = 0; i < _HC_Ly - 1; ++i)
            {
                appGeneral(_T("{"));
                for (UINT j = 0; j < _HC_Lx - 1; ++j)
                {
                    const UINT idx = k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j;

                    if (0 == k)
                    {
                        tmpjgs.AddItem(arrayRes[idx]);
                    }
                    else
                    {
                        tmpjgs[i * (_HC_Lx - 1) + j] += arrayRes[idx];
                    }

                    if (0 == j)
                    {
                        appGeneral(_T("%s"), appToString(arrayRes[idx]).c_str());
                    }
                    else
                    {
                        appGeneral(_T(", %s"), appToString(arrayRes[idx]).c_str());
                    }
                }
                appGeneral(_T("}, "));
            }
            appGeneral(_T("}\n"));
        }
        appGeneral(_T("}\n"));

        appGeneral(_T("\n -------------------- Average -------------------------\n\n"));

        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                appGeneral(_T("(x=%d,y=%d)%s,   "),
                    j + 1, i + 1,
                    appToString(_divCHost(tmpjgs[i * (_HC_Lx - 1) + j], uiConfig)).c_str());
            }
            appGeneral(_T("\n"));
        }
    }

    template<class T>
    static void _ZeroXYPlane(T* pDeviceRes);
    /**
    * array[x, y] = array[x, y] / (lz * lt)
    */
    template<class T>
    static void _AverageXYPlane(T* pDeviceRes);

    template<class T>
    static void _ZeroSlice(T* pDeviceRes, BYTE byDir);

    template<class T>
    static void XYDataToRdistri(
        UBOOL bShiftCenter,
        const T* source,
        UINT* count,
        T* result,
        UINT uiMaxR,
        UBOOL bCalculateCounter,
        BYTE byFieldId);

    template<class T>
    static void ReportDistributeWithR(UINT uiConf, UINT uiR, const TArray<T>& arrayData)
    {
        appAssert(uiConf * uiR == static_cast<UINT>(arrayData.GetCount()));
        appGeneral(_T("{\n"));
        for (UINT conf = 0; conf < uiConf; ++conf)
        {
            for (UINT r = 0; r < uiR; ++r)
            {
                if (0 == r)
                {
                    appGeneral(_T("{ %s"), appToString(arrayData[uiR * conf + r]).c_str());
                }
                else
                {
                    appGeneral(_T(", %s"), appToString(arrayData[uiR * conf + r]).c_str());
                }
            }

            appGeneral(_T("},\n"));
        }
        appGeneral(_T("}\n"));
    }

    /**
     * TransformFromXYDataToRDataOnce_C and TransformFromXYDataToRDataOnce_R
     * is for gauge measurement
     */
    template<class T>
    static void TransformFromXYDataToRDataOnce(
        UBOOL bShiftCenter,
        const T* __restrict__ pXYData,
        UINT* pCountBuffer,
        T* pValueBuffer,
        UINT* pHostCountBuffer,
        T* pHostValueBuffer,
        UINT uiMaxR,
        UINT uiEdgeR,
        UBOOL bCalculateCounter,
        BYTE byFieldId,
        TArray<T>& arrData,
        TArray<T>* arrInner,
        TArray<T>& arrFull,
        TArray<UINT>& arrR,
        UINT uiConfig,
        Real fDivider)
    {
        XYDataToRdistri(
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
        
        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(T) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        //Here we have already divide by all XYZ points
        FillDataWithR(
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
     * 
     * Sometimes, we need to set bMinus = TRUE, because
     * <qbar M q> = - tr[MD^{-1}]
     * but tr[MD^{-1}] is measured
     * 
     * NOTE: For chiral condensation, it is usually defined as tr[D^{-1}], not -tr[D^{-1}].
     * NOTE: Usually, tr[D^{-1}] > 0, and -tr[D^{-1}] < 0.
     */
    template<class T>
    static void TransformFromXYDataToRData(
        UBOOL bMinus,
        UBOOL bShiftCenter,
        UINT uiMaxR,
        UINT uiEdgeR,
        BYTE byFieldId,
        UINT uiFieldCount,
        UINT uiMeasureCount,
        UINT uiConfig,
        const T* const* pXYBuffers,
        UINT* pCountBuffer,
        T* pValueBuffer,
        UINT* pHostCountBuffer,
        T* pHostValueBuffer,
        TArray<UINT>& lstR,
        TArray<T>* lstValues,
        TArray<T>* lstAll,
        TArray<T>* lstInner)
    {
        for (UINT i = 0; i < uiMeasureCount; ++i)
        {
            XYDataToRdistri(bShiftCenter, pXYBuffers[i], pCountBuffer, pValueBuffer, uiMaxR, 0 == i, byFieldId);
            if (0 == i)
            {
                checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
            }

            checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(T) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

            FillDataWithR(
                lstValues[i],
                NULL == lstInner ? NULL : &(lstInner[i]),
                lstAll[i],
                lstR,
                pHostValueBuffer,
                pHostCountBuffer,
                uiConfig,
                uiMaxR,
                uiEdgeR,
                (bMinus ? F(-1.0) : F(1.0)) / static_cast<Real>(uiFieldCount * _HC_Lz * _HC_Lt),
                0 == i
            );
        }
    }

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
    BYTE m_byFermionFieldId;
    UBOOL m_bShowResult;
    TArray<BYTE> m_lstGaugeFieldIds;
    TArray<BYTE> m_lstBosonFieldIds;

public:
    //============================================================
    //some simple measurement only produce real or complex results
    UINT GetConfigurationCount() const 
    {
        return m_uiConfigurationCount;
    }

    Real GetLastRealRes() const
    {
        appAssert(m_lstRealResults.Num() > 0);
        return (m_lstRealResults.Num() > 0) ? m_lstRealResults[m_lstRealResults.Num() - 1] : F(0.0);
    }

    CLGComplex GetLastCmpRes() const
    {
        appAssert(m_lstComplexResults.Num() > 0);
        return (m_lstComplexResults.Num() > 0) ? m_lstComplexResults[m_lstComplexResults.Num() - 1] : _zeroc;
    }

    //========================
    //Make sure average has been called so that average is calculated.
    Real GetAverageRealRes() const
    {
        return m_fAverageRealRes;
    }

    CLGComplex GetAverageCmpRes() const
    {
        return m_cAverageCmpRes;
    }

    Real RealResAtI(INT i) const
    {
        appAssert(i < m_lstRealResults.Num() && i >= 0);
        return (i < m_lstRealResults.Num() && i >= 0) ? m_lstRealResults[i] : F(0.0);
    }

    CLGComplex CmpResAtI(INT i) const
    {
        appAssert(i < m_lstComplexResults.Num() && i >= 0);
        return (i < m_lstComplexResults.Num() && i >= 0) ? m_lstComplexResults[i] : _zeroc;
    }

    void WriteRealListToFile(const CCString& sFileName) const;
    void WriteCmpListToFile(const CCString& sFileName) const;

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CBase::GetInfos(tab);
        sRet = sRet + tab + _T("MeasureId : ") + appToString(m_byId) + _T("\n");
        sRet = sRet + tab + _T("Smearing : ") + appToString(m_bNeedSmearing) + _T("\n");
        sRet = sRet + tab + _T("FermionFieldId : ") + appToString(m_byFermionFieldId) + _T("\n");
        sRet = sRet + tab + _T("GaugeFields : ") + appToString(m_lstGaugeFieldIds) + _T("\n");
        sRet = sRet + tab + _T("BosonFields : ") + appToString(m_lstBosonFieldIds) + _T("\n");
        return sRet;
    }

protected:

    //P4-3.1: reduce a locally-measured DOUBLE (e.g. CFieldGauge::CalculatePlaqutteEnergy,
    //which returns the local partial sum on multi-GPU) to the GLOBAL value on every
    //rank. Survey note: the HMC plaquette-energy path (CActionGaugePlaquette::
    //EnergySingleField) was already Allreduced in P4-1.2+, and measurements here do
    //NOT go through CAction::Energy, so calling this on the measured value cannot
    //double-reduce. Implemented in CMeasure.cu (needs appGetComm()).
    void GlobalSumReal(DOUBLE& fValue) const;

    //P4-3.5 (code only, untested): array variant for per-time-slice profiles.
    //Reduce a locally-measured DOUBLE array (uiCount entries, e.g. Lt slices) to
    //the GLOBAL values on every rank. Implemented in CMeasure.cu; test pending
    //(P4-3.5 unit-level 1-vs-N).
    void GlobalSumRealArray(DOUBLE* pValues, UINT uiCount) const;

    //P4-3.2: plaquette normalization count. On multi-GPU a globally-reduced
    //plaquette ENERGY must be normalized by the GLOBAL plaquette count (the
    //local _HC_PlaqutteCount is the per-rank share); single-GPU is unchanged.
    DOUBLE GlobalPlaqutteCount() const;

    //P4-3.3: global lattice length in direction uiDir (0..3). On multi-GPU the
    //local _HC_Lx.._HC_Lt are the per-rank share; position/profile normalization
    //factors of globally-reduced results must use the GLOBAL lengths.
    //Single-GPU (no comm) returns the local value, which is the whole lattice.
    DOUBLE GlobalL(UINT uiDir) const;

    //P4-3.4: reduce a locally-measured CLGComplex array (uiCount entries) to the
    //GLOBAL values on every rank. CLGComplex is cuDoubleComplex on double builds
    //and cuComplex on float builds; the reduction goes through a DOUBLE staging
    //buffer in the float case (MPI double is the portable choice). No-op on
    //single-GPU / unsplit builds.
    void GlobalSumComplexArray(CLGComplex* pValues, UINT uiCount) const;

    //I9: scalar variant of GlobalSumComplexArray.
    void GlobalSumComplex(CLGComplex& fValue) const;

    //I9: Real-array variant, declared only on float builds: on double builds
    //Real == DOUBLE, so the DOUBLE* overload above already applies (same pattern
    //as CLGComm::AllreduceSum). Staged through DOUBLE in the implementation.
#if !_CLG_DOUBLEFLOAT
    void GlobalSumRealArray(Real* pValues, UINT uiCount) const;
#endif

    void UpdateRealResult(Real fResult, UBOOL bUpdateConfigurationCount = TRUE)
    {
        if (bUpdateConfigurationCount)
        {
            ++m_uiConfigurationCount;
        }
        m_lstRealResults.AddItem(fResult);
    }

    void UpdateComplexResult(CLGComplex fResult, UBOOL bUpdateConfigurationCount = TRUE)
    {
        if (bUpdateConfigurationCount)
        {
            ++m_uiConfigurationCount;
        }
        m_lstComplexResults.AddItem(fResult);
    }

    void ReportAverageComplexRes() const
    {
        LogGeneralComplex(m_cAverageCmpRes);
    }

    UINT m_uiConfigurationCount;

private:

    TArray<Real> m_lstRealResults;
    TArray<CLGComplex> m_lstComplexResults;
    Real m_fAverageRealRes;
    CLGComplex m_cAverageCmpRes;
    
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


    UBOOL IsGaugeOrBosonMeasurement() const override = 0;
    UBOOL IsSourceScanning() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }
    UINT GetFieldCount() const { return m_uiFieldCount; }
    void SetFieldCount(UINT uiFieldCount) { m_uiFieldCount = uiFieldCount; }

    virtual TArray<TArray<CLGComplex>> ExportDiagnal(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, class CFieldFermion* pooled1, class CFieldFermion* pooled2)
    {
        TArray<TArray<CLGComplex>> ret;
        appCrucial(_T("ExportDiagnal not implemented\n"));
        return ret;
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CMeasure::GetInfos(tab);
        sRet = sRet + tab + _T("FieldCount : ") + appToString(m_uiFieldCount) + _T("\n");
        return sRet;
    }

protected:

    UINT m_uiFieldCount;
    UBOOL m_bDebugDivation;

};

#include "CMeasureData.h"


__END_NAMESPACE

//=============== some widely used macros in applications ======================

#define _CLG_EXPORT_CHIRAL(measureName, lstName, variableName) \
CCString sFileNameWrite##measureName##lstName = _T("%s_%d_condensate"); \
CCString sFileNameWrite##measureName##lstName##XSlice = _T("%s_%d_condensateXSlice"); \
CCString sFileNameWrite##measureName##lstName##YSlice = _T("%s_%d_condensateYSlice"); \
CCString sFileNameWrite##measureName##lstName##ZSlice = _T("%s_%d_condensateZSlice"); \
CCString sFileNameWrite##measureName##lstName##TSlice = _T("%s_%d_condensateTSlice"); \
CCString sFileNameWrite##measureName##lstName##In = _T("%s_%d_condensate"); \
CCString sFileNameWrite##measureName##lstName##OverR = _T("%s_%d_condensate"); \
sFileNameWrite##measureName##lstName = sFileNameWrite##measureName##lstName + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##XSlice = sFileNameWrite##measureName##lstName##XSlice + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##YSlice = sFileNameWrite##measureName##lstName##YSlice + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##ZSlice = sFileNameWrite##measureName##lstName##ZSlice + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##TSlice = sFileNameWrite##measureName##lstName##TSlice + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##In = sFileNameWrite##measureName##lstName##In + _T(#measureName) + _T(#lstName) + _T("_In.csv"); \
sFileNameWrite##measureName##lstName##OverR = sFileNameWrite##measureName##lstName##OverR + _T(#measureName) + _T(#lstName) + _T("_OverR.csv"); \
sFileNameWrite##measureName##lstName.Format(sFileNameWrite##measureName##lstName, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##XSlice.Format(sFileNameWrite##measureName##lstName##XSlice, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##YSlice.Format(sFileNameWrite##measureName##lstName##YSlice, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##ZSlice.Format(sFileNameWrite##measureName##lstName##ZSlice, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##TSlice.Format(sFileNameWrite##measureName##lstName##TSlice, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##In.Format(sFileNameWrite##measureName##lstName##In, sCSVSavePrefix.c_str(), variableName); \
sFileNameWrite##measureName##lstName##OverR.Format(sFileNameWrite##measureName##lstName##OverR, sCSVSavePrefix.c_str(), variableName); \
TArray<CLGComplex> lstName##measureName; \
TArray<CLGComplex> lstName##measureName##In; \
TArray<TArray<CLGComplex>> lstName##measureName##XSlice; \
TArray<TArray<CLGComplex>> lstName##measureName##YSlice; \
TArray<TArray<CLGComplex>> lstName##measureName##ZSlice; \
TArray<TArray<CLGComplex>> lstName##measureName##TSlice; \
TArray<TArray<CLGComplex>> lstName##measureName##OverR; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    lstName##measureName.AddItem(measureName->m_lstCondAll[lstName][j]); \
    lstName##measureName##In.AddItem(measureName->m_lstCondIn[lstName][j]); \
    TArray<CLGComplex> thisConfigurationOverR; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfigurationOverR.AddItem(measureName->m_lstCond[lstName][j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##measureName##OverR.AddItem(thisConfigurationOverR); \
    if (measureName->m_bMeasureXSlice) \
    { \
        TArray<CLGComplex> thisConfiguration##measureName##lstName##XSlice; \
        for (UINT i = 0; i < _HC_Lx; ++i) \
        { \
            thisConfiguration##measureName##lstName##XSlice.AddItem(measureName->m_lstCondXSlice[lstName][j * _HC_Lx + i]); \
        } \
        lstName##measureName##XSlice.AddItem(thisConfiguration##measureName##lstName##XSlice); \
    } \
    if (measureName->m_bMeasureYSlice) \
    { \
        TArray<CLGComplex> thisConfiguration##measureName##lstName##YSlice; \
        for (UINT i = 0; i < _HC_Ly; ++i) \
        { \
            thisConfiguration##measureName##lstName##YSlice.AddItem(measureName->m_lstCondYSlice[lstName][j * _HC_Ly + i]); \
        } \
        lstName##measureName##YSlice.AddItem(thisConfiguration##measureName##lstName##YSlice); \
    } \
    if (measureName->m_bMeasureZSlice) \
    { \
        TArray<CLGComplex> thisConfiguration##measureName##lstName##ZSlice; \
        for (UINT i = 0; i < _HC_Lz; ++i) \
        { \
            thisConfiguration##measureName##lstName##ZSlice.AddItem(measureName->m_lstCondZSlice[lstName][j * _HC_Lz + i]); \
        } \
        lstName##measureName##ZSlice.AddItem(thisConfiguration##measureName##lstName##ZSlice); \
    } \
    if (measureName->m_bMeasureTSlice) \
    { \
        TArray<CLGComplex> thisConfiguration##measureName##lstName##TSlice; \
        for (UINT i = 0; i < _HC_Lt; ++i) \
        { \
            thisConfiguration##measureName##lstName##TSlice.AddItem(measureName->m_lstCondTSlice[lstName][j * _HC_Lt + i]); \
        } \
        lstName##measureName##TSlice.AddItem(thisConfiguration##measureName##lstName##TSlice); \
    } \
} \
WriteComplexArray(sFileNameWrite##measureName##lstName, lstName##measureName); \
WriteComplexArray(sFileNameWrite##measureName##lstName##In, lstName##measureName##In); \
WriteComplexArray2(sFileNameWrite##measureName##lstName##OverR, lstName##measureName##OverR); \
if (measureName->m_bMeasureXSlice) \
{ \
    WriteComplexArray2(sFileNameWrite##measureName##lstName##XSlice, lstName##measureName##XSlice); \
} \
if (measureName->m_bMeasureYSlice) \
{ \
    WriteComplexArray2(sFileNameWrite##measureName##lstName##YSlice, lstName##measureName##YSlice); \
} \
if (measureName->m_bMeasureZSlice) \
{ \
    WriteComplexArray2(sFileNameWrite##measureName##lstName##ZSlice, lstName##measureName##ZSlice); \
} \
if (measureName->m_bMeasureTSlice) \
{ \
    WriteComplexArray2(sFileNameWrite##measureName##lstName##TSlice, lstName##measureName##TSlice); \
}
#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================