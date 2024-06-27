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
        , m_pTmpDeviceRes(NULL)
        , m_pDevicePath(NULL)
    {

    }

    ~CMeasureWilsonLoopWithPath();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    TArray<INT> m_lstPath;
    SSmallInt4 m_sPoint;
    UBOOL m_bAllPoint;
    CLGComplex* m_pTmpDeviceRes;
    INT* m_pDevicePath;

};



__END_NAMESPACE

#endif //#ifndef _CMEASUREWILSONLOOPWITHPATH_H_

//=============================================================================
// END OF FILE
//=============================================================================