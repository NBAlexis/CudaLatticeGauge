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
        , m_pHostXYBuffer(NULL)
        , m_pHostXBuffer(NULL)
        , m_pHostYBuffer(NULL)
        , m_pHostZBuffer(NULL)
        , m_pHostTBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_uiEdge(1)
        , m_bShiftCenter(FALSE)
        , m_bMeasureSigma12(FALSE)
        
        , m_bMeasureConnect(FALSE)
        , m_bMeasureXSlice(FALSE)
        , m_bMeasureYSlice(FALSE)
        , m_bMeasureZSlice(FALSE)
        , m_bMeasureTSlice(FALSE)
    {
        
    }

    ~CMeasureChiralCondensateKS();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;

    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }    

    TArray<TArray<CLGComplex>> ExportDiagnal(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, class CFieldFermion* pooled1, class CFieldFermion* pooled2) override;

protected:
    
    CLGComplex* m_pDeviceXYBuffer[ChiralKSMax];
    CLGComplex* m_pDeviceXBuffer[ChiralKSMax];
    CLGComplex* m_pDeviceYBuffer[ChiralKSMax];
    CLGComplex* m_pDeviceZBuffer[ChiralKSMax];
    CLGComplex* m_pDeviceTBuffer[ChiralKSMax];

    CLGComplex* m_pHostXYBuffer;
    CLGComplex* m_pHostXBuffer;
    CLGComplex* m_pHostYBuffer;
    CLGComplex* m_pHostZBuffer;
    CLGComplex* m_pHostTBuffer;
    //CLGComplex m_cTmpSum[ChiralKSMax];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UINT m_uiEdge;
    UBOOL m_bShiftCenter;
    UBOOL m_bMeasureSigma12;

public:

    UBOOL m_bMeasureConnect;
    UBOOL m_bMeasureXSlice;
    UBOOL m_bMeasureYSlice;
    UBOOL m_bMeasureZSlice;
    UBOOL m_bMeasureTSlice;
    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[ChiralKSMax];
    TArray<CLGComplex> m_lstCondIn[ChiralKSMax];
    TArray<CLGComplex> m_lstCond[ChiralKSMax];
    TArray<CLGComplex> m_lstCondXSlice[ChiralKSMax];
    TArray<CLGComplex> m_lstCondYSlice[ChiralKSMax];
    TArray<CLGComplex> m_lstCondZSlice[ChiralKSMax];
    TArray<CLGComplex> m_lstCondTSlice[ChiralKSMax];
    TArray<CLGComplex> m_lstDebugData[ChiralKSMax];
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKS_H_

//=============================================================================
// END OF FILE
//=============================================================================