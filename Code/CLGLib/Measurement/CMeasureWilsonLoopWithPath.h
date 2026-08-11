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

    enum { _kMaxWilsonPathLength = 128 };

    CMeasureWilsonLoopWithPath()
        : CMeasure()
        , m_bAllPoint(TRUE)
        , m_pTmpDeviceRes(NULL)
        , m_pDevicePath(NULL)
        , m_pDevicePath2(NULL)
        , m_bTwoPathBackForward(FALSE)
    {

    }

    ~CMeasureWilsonLoopWithPath();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    void SetPath(const TArray<TArray<SCHAR>>& sp)
    {
        m_lstPath = sp;
    }

    void Reset() override
    {
        CMeasure::Reset();
        m_lstV.RemoveAll();
        m_lstDV.RemoveAll();
    }

    void SetAsTwoPathBackForward(UBOOL bTwoPath, const TArray<TArray<SCHAR>>& path, const TArray<SSmallInt4>& shift)
    {
        m_bTwoPathBackForward = bTwoPath;
        m_lstPath = path;
        m_sShift = shift;
    }

protected:

    TArray<TArray<SCHAR>> m_lstPath;
    SSmallInt4 m_sPoint;
    UBOOL m_bAllPoint;
    cuDoubleComplex* m_pTmpDeviceRes;
    SCHAR* m_pDevicePath;
    SCHAR* m_pDevicePath2;

    UBOOL m_bTwoPathBackForward;
    TArray<SSmallInt4> m_sShift;

public:

    //every path
    TArray<TArray<cuDoubleComplex>> m_lstV;
    TArray<TArray<DOUBLE>> m_lstDV;
};



__END_NAMESPACE

#endif //#ifndef _CMEASUREWILSONLOOPWITHPATH_H_

//=============================================================================
// END OF FILE
//=============================================================================