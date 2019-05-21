//=============================================================================
// FILENAME : CMeasureAMomemtumJG.h
// 
// DESCRIPTION:
// This is measurement for angular momentum of rotating frame
// The angular momentum of each site is calculated, average is taken over z and t directions
// Then the average is taken over all configurations (result in angular momentum on a x-y plane)
//
// We assume lx * ly < max-thread
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#ifndef _CMEASUREAMOMENTUMJG_H_
#define _CMEASUREAMOMENTUMJG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAMomentumJG)

//=================== make them static functions because they are used also in CMeasureAMomentumJF ===================
/**
* array[x, y] = array[x, y] / (lz * lt)
*/
static void _AverageXYPlane(Real* pDeviceRes);

/**
* array[x, y] = 0
*/
static void _ZeroXYPlane(Real* pDeviceRes);

/**
* array[x, y] = (array[x, y] * (N-1) + oneconfiguration[x, y]) / N
*/
static void _AverageXYPlaneOverConf(Real* pDeviceRes, const Real* __restrict__ pDeviceResOneConfig, UINT uiConfigCount);

class CLGAPI CMeasureAMomentumJG : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumJG)
public:
    CMeasureAMomentumJG() 
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBuffer(NULL) 
        , m_pDeviceDataBufferOneConfig(NULL)
        , m_byFieldId(1)
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasureAMomentumJG();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

protected:

    Real * m_pHostDataBuffer;
    Real * m_pDeviceDataBuffer;
    Real * m_pDeviceDataBufferOneConfig;
    SSmallInt4 m_sCenter;
    BYTE m_byFieldId;
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMJG_H_

//=============================================================================
// END OF FILE
//=============================================================================