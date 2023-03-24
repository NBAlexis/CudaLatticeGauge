//=============================================================================
// FILENAME : CMeasureChiralCondensateKS.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [10/01/2020 nbale]
//=============================================================================

#ifndef _CMEASURECHIRALCONDENSATEKS_H_
#define _CMEASURECHIRALCONDENSATEKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureChiralCondensateKS)

enum EChiralMeasureTypeKS
{
    ChiralKS = 0,
    ConnectSusp = 1,

    CMTKSGamma1,
    CMTKSGamma2,
    CMTKSGamma3,
    CMTKSGamma4,
    CMTKSGamma5,
    CMTKSGamma51,
    CMTKSGamma52,
    CMTKSGamma53,
    CMTKSGamma54,
    CMTKSSigma12,
    CMTKSSigma13,
    CMTKSSigma14,
    CMTKSSigma23,
    CMTKSSigma24,
    CMTKSSigma34,

    ChiralKSMax,
};

class CLGAPI CMeasureChiralCondensateKS : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureChiralCondensateKS)
public:

    CMeasureChiralCondensateKS()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pHostXYBuffer(NULL)
        , m_pHostZBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_uiEdge(1)
        , m_bShiftCenter(FALSE)
        , m_bMeasureSigma12(FALSE)
        
        , m_bShowResult(FALSE)
        , m_bMeasureConnect(FALSE)
        , m_bMeasureZSlice(FALSE)
    {
        
    }

    ~CMeasureChiralCondensateKS();

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

    CLGComplex* m_pDeviceXYBuffer[ChiralKSMax];
    CLGComplex* m_pDeviceZBuffer[ChiralKSMax];
    CLGComplex* m_pHostXYBuffer;
    CLGComplex* m_pHostZBuffer;
    //CLGComplex m_cTmpSum[ChiralKSMax];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UINT m_uiEdge;
    UBOOL m_bShiftCenter;
    UBOOL m_bMeasureSigma12;
    UBOOL m_bShowResult;

public:

    UBOOL m_bMeasureConnect;
    UBOOL m_bMeasureZSlice;
    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[ChiralKSMax];
    TArray<CLGComplex> m_lstCondIn[ChiralKSMax];
    TArray<CLGComplex> m_lstCond[ChiralKSMax];
    TArray<CLGComplex> m_lstCondZSlice[ChiralKSMax];
    TArray<CLGComplex> m_lstDebugData[ChiralKSMax];
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKS_H_

//=============================================================================
// END OF FILE
//=============================================================================