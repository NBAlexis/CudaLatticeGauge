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

#if !_CLG_WIN

inline void strerror_s(TCHAR* buffer, size_t bufferSize, INT error)
{
    strcpy(buffer, strerror(error));
}

inline void _gcvt_s(TCHAR* buff, UINT uiBuffLength, Real fVaule, UINT uiDigit)
{
    static TCHAR tmpBuff[10];
    appSprintf(tmpBuff, 10, _T("%s.%df"), _T("%"), uiDigit);
    appSprintf(buff, uiBuffLength, tmpBuff, fVaule);
}

#endif

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
    virtual void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple);


    /**
    * NOTE: sources will be passed to multiple measures, do NOT change the content!
    * NOTE: site.x start from 1 to Lx - 1, 0 is not included
    */
    virtual void SourceSanning(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site);

    /**
    * Z4 Source
    */
    virtual void OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd);

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
     * 
     * Sometimes, we need to set bMinus = TRUE, because
     * <qbar M q> = - tr[MD^{-1}]
     * but tr[MD^{-1}] is measured
     */
    static void TransformFromXYDataToRData_C(
        UBOOL bMinus,
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
        assert(m_lstRealResults.Num() > 0);
        return (m_lstRealResults.Num() > 0) ? m_lstRealResults[m_lstRealResults.Num() - 1] : F(0.0);
    }

    CLGComplex GetLastCmpRes() const
    {
        assert(m_lstComplexResults.Num() > 0);
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
        assert(i < m_lstRealResults.Num() && i >= 0);
        return (i < m_lstRealResults.Num() && i >= 0) ? m_lstRealResults[i] : F(0.0);
    }

    CLGComplex CmpResAtI(INT i) const
    {
        assert(i < m_lstComplexResults.Num() && i >= 0);
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

    virtual TArray<TArray<CLGComplex>> ExportDiagnal(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, class CFieldFermion* pooled1, class CFieldFermion* pooled2)
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

enum { kFileDigital = 20, };

template <class T>
void WriteRealArray(const CCString& sFileName, const TArray<T>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, static_cast<DOUBLE>(lst[i]), iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;
        if (i != lst.GetCount() - 1)
        {
            file << _T(",");
        }
    }
    file.flush();
    file.close();
}

template <class T>
void WriteRealArray2(const CCString& sFileName, const TArray<TArray<T>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j], iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            file << _T(" ");
            file << sReal;
            if (j != lst[i].GetCount() - 1)
            {
                file << _T(",");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray(const CCString& sFileName, const TArray<CLGComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }

        file << _T(" ");
        file << sReal;
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}

#if !_CLG_DOUBLEFLOAT
inline void WriteComplexArray(const CCString& sFileName, const TArray<cuDoubleComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }

        file << _T(" ");
        file << sReal;
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}
#endif

inline void WriteComplexArray2(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            _gcvt_s(str, 50, lst[i][j].y, iDigital);
            CCString sImg = CCString(str);
            sImg = sImg.Replace(_T("e"), _T("*^"));
            CCString sMid = _T(" + ");
            if (sImg.Left(1) == _T("-"))
            {
                sImg = sImg.Right(sImg.GetLength() - 1);
                sMid = _T(" - ");
            }
            file << _T(" ");
            file << sReal;
            file << sMid;
            file << sImg;
            if (j == lst[i].GetCount() - 1)
            {
                file << _T(" I");
            }
            else
            {
                file << _T(" I,");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray2Simple(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = 8;
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));

            file << _T(" ");
            file << sReal;

            if (abs(lst[i][j].y) > F(0.0001) * lst[i][j].x)
            {
                _gcvt_s(str, 50, lst[i][j].y, iDigital);
                CCString sImg = CCString(str);
                sImg = sImg.Replace(_T("e"), _T("*^"));
                CCString sMid = _T(" + ");
                if (sImg.Left(1) == _T("-"))
                {
                    sImg = sImg.Right(sImg.GetLength() - 1);
                    sMid = _T(" - ");
                }

                file << sMid;
                file << sImg;
                if (j == lst[i].GetCount() - 1)
                {
                    file << _T(" I");
                }
                else
                {
                    file << _T(" I,");
                }
            }
            else
            {
                if (j != lst[i].GetCount() - 1)
                {
                    file << _T(",");
                }
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

__END_NAMESPACE

#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================