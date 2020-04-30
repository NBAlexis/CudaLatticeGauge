//=============================================================================
// FILENAME : CMeasureChiralCondensate.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURECHIRALCONDENSATE_H_
#define _CMEASURECHIRALCONDENSATE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureChiralCondensate)

class CLGAPI CMeasureChiralCondensate : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureChiralCondensate)
public:

    enum { _kCondMeasureCount = 9 };

    CMeasureChiralCondensate()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pHostXYBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_bMeasureDistribution(FALSE)
        , m_bShowResult(FALSE)
    {
        
    }

    ~CMeasureChiralCondensate();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }

protected:
    
    UINT m_uiConfigurationCount;

    CLGComplex* m_pDeviceXYBuffer[_kCondMeasureCount];
    CLGComplex* m_pHostXYBuffer;
    CLGComplex m_cTmpSum[_kCondMeasureCount];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[_kCondMeasureCount];
    TArray<CLGComplex> m_lstCond[_kCondMeasureCount];  
};


#pragma region Measure functions

typedef void(*_deviceMeasureCondensateFunc)(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element);

/**
 * Psibar Psi
 */
__device__ __inline__ void _deviceMeasureChiral(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    
}

/**
 * Psibar gamma_i Psi
 */
__device__ __inline__ void _deviceMeasureGamma1(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA1].MulWilsonC(element);
}

__device__ __inline__ void _deviceMeasureGamma2(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA2].MulWilsonC(element);
}

__device__ __inline__ void _deviceMeasureGamma3(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA3].MulWilsonC(element);
}

/**
* This is also the Rho condensate
*/
__device__ __inline__ void _deviceMeasureGamma4(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA4].MulWilsonC(element);
}

/**
* Pion condensate
* Psibar gamma5 Psi
*/
__device__ __inline__ void _deviceMeasureGamma5(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA5].MulWilsonC(element);
}

/**
 * Psibar gamma4 gamma5 Psi
 * Also known as n5
 * arXiv:1105.0385
 */
__device__ __inline__ void _deviceMeasureGamma45(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    element = __chiralGamma[GAMMA5].MulWilsonC(element);
    element = __chiralGamma[GAMMA4].MulWilsonC(element);
}

/**
* gamma_x
* gamma_1 + y Omega gamma_4
*/
__device__ __inline__ void _deviceMeasureGammaX(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    const Real fYOmega = static_cast<Real>(site.y - sCenter.y)* fOmega;
    deviceWilsonVectorSU3 toAdd = __chiralGamma[GAMMA4].MulWilsonC(element);
    toAdd.MulReal(fYOmega);
    element = __chiralGamma[GAMMA1].MulWilsonC(element);
    element.Add(toAdd);
}

/**
* gamma_y 
* gamma_2 - x Omega gamma_4
*/
__device__ __inline__ void _deviceMeasureGammaY(const SSmallInt4& site, const SSmallInt4& sCenter, Real fOmega, deviceWilsonVectorSU3& element)
{
    const Real fXOmega = static_cast<Real>(site.x - sCenter.x)* fOmega;
    deviceWilsonVectorSU3 toAdd = __chiralGamma[GAMMA4].MulWilsonC(element);
    toAdd.MulReal(fXOmega);
    element = __chiralGamma[GAMMA2].MulWilsonC(element);
    element.Sub(toAdd);
}


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATE_H_

//=============================================================================
// END OF FILE
//=============================================================================