//=============================================================================
// FILENAME : CMeasureWilsonLoopWithPath.h
// 
// DESCRIPTION:
//
// REVISION:
//  [05/10/2021 nbale]
//=============================================================================

#ifndef _CMEASUREWILSONLOOPWITHPATH_H_
#define _CMEASUREWILSONLOOPWITHPATH_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CMeasureWilsonLoopWithPath)

class CLGAPI CMeasureWilsonLoopWithPath : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureWilsonLoopWithPath)

public:

    CMeasureWilsonLoopWithPath()
        : CMeasure()
        , m_bAllPoint(TRUE)
        , m_bShowResult(FALSE)

        , m_pTmpDeviceRes(NULL)
        , m_pDevicePath(NULL)

        , m_uiConfigurationCount(0)
    {
    }

    ~CMeasureWilsonLoopWithPath();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    TArray<INT> m_lstPath;
    SSmallInt4 m_sPoint;
    UBOOL m_bAllPoint;
    UBOOL m_bShowResult;
    CLGComplex* m_pTmpDeviceRes;
    INT* m_pDevicePath;
    UINT m_uiConfigurationCount;

public:

    TArray<CLGComplex> m_lstRes;
};



__END_NAMESPACE

#endif //#ifndef _CMEASUREWILSONLOOPWITHPATH_H_

//=============================================================================
// END OF FILE
//=============================================================================