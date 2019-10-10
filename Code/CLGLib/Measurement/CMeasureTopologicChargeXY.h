//=============================================================================
// FILENAME : CMeasureTopologicChargeXY.h
// 
// DESCRIPTION:
// NOTE: This is not ready yet!! I think the definition of topological charge maybe different in rotating frame
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURETOPOLOGICALCHARGEXY_H_
#define _CMEASURETOPOLOGICALCHARGEXY_H_

#define OneOver32PI2 (F(1.0) / (F(32.0) * PISQ))

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureTopologicChargeXY)

class CLGAPI CMeasureTopologicChargeXY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureTopologicChargeXY)
public:
    CMeasureTopologicChargeXY()
        : CMeasure()
        , m_pXYHostDensity(NULL)
        , m_pXYDeviceDensity(NULL)
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasureTopologicChargeXY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    Real* m_pXYHostDensity;
    Real* m_pXYDeviceDensity;
    
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    TArray<Real> m_lstCharge;
    TArray<Real> m_lstXYDensity;
};

#pragma region device functions

/**
 * U_{mu,nu}(n)+U^+_{-mu,nu}(n)+U^+_{mu,-nu}(n)+U_{-mu,-nu}(n)
 * or
 * U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu)
 */
static __device__ __inline__ deviceSU3 _deviceClover(const deviceSU3* __restrict__ pGaugeField, UINT uiBigIdx, BYTE mu, BYTE nu)
{
    deviceSU3 ret(_device1PlaqutteTermPP(pGaugeField, mu, nu, uiBigIdx));
    ret.Add(_device1PlaqutteTermMM(pGaugeField, mu, nu, uiBigIdx));
    ret.AddDagger(_device1PlaqutteTermPM(pGaugeField, nu, mu, uiBigIdx));
    ret.AddDagger(_device1PlaqutteTermMP(pGaugeField, nu, mu, uiBigIdx));

    return ret;
}

static __device__ __inline__ Real _deviceTrImClover(const deviceSU3* __restrict__ pGaugeField, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE rho, BYTE sigma)
{
    return deviceSU3::TrIm(
        _deviceClover(pGaugeField, uiBigIdx, mu, nu),
        _deviceClover(pGaugeField, uiBigIdx, rho, sigma));
}

static __device__ __inline__ Real _deviceTopologicalCharge(const deviceSU3* __restrict__ pGaugeField, UINT uiBigIdx)
{
    Real ret = _deviceTrImClover(pGaugeField, uiBigIdx, 0, 1, 2, 3);
    ret -= _deviceTrImClover(pGaugeField, uiBigIdx, 0, 2, 1, 3);
    ret += _deviceTrImClover(pGaugeField, uiBigIdx, 0, 3, 1, 2);
    return OneOver32PI2 * F(2.0) * ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CMEASURETOPOLOGICALCHARGEXY_H_

//=============================================================================
// END OF FILE
//=============================================================================