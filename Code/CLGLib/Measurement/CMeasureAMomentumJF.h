//=============================================================================
// FILENAME : CMeasureAMomemtumJG.h
// 
// DESCRIPTION:
// This is measurement for angular momentum of rotating frame
// The angular momentum of each site of two lines is calculated.
// It is impractical to calculate the angular momentum for each site (a 16^4 lattice has 65536 sites)
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#ifndef _CMEASUREAMOMENTUMJF_H_
#define _CMEASUREAMOMENTUMJF_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAMomentumJF)

class CLGAPI CMeasureAMomentumJF : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumJF)
public:
    CMeasureAMomentumJF()
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBufferS(NULL) 
        , m_pDeviceDataBufferL(NULL)
        , m_pOperatorDataS(NULL)
        , m_pOperatorDataL(NULL)
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
        , m_bNaive(TRUE)
        , m_bExponential(FALSE)
    {
    }

    ~CMeasureAMomentumJF();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site);
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return FALSE; }
    virtual UBOOL IsSourceScanning() const { return TRUE; }

protected:

    CLGComplex * m_pHostDataBuffer;
    CLGComplex * m_pDeviceDataBufferS;
    CLGComplex * m_pDeviceDataBufferL;
    deviceWilsonVectorSU3* m_pOperatorDataS;
    deviceWilsonVectorSU3* m_pOperatorDataL;
    
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    UBOOL m_bNaive;
    UBOOL m_bExponential;
    TArray<Real> m_lstAllRes;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMJG_H_

//=============================================================================
// END OF FILE
//=============================================================================