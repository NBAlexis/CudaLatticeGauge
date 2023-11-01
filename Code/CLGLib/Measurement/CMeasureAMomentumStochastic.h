//=============================================================================
// FILENAME : CMeasureAMomentumStochastic.h
// 
// DESCRIPTION:
// NOTE:
//          2 is added because of Nf=2
// 
// JL :     Ji decomposition orbit angular momentum, (1/i) psidagger r x D psi
//             It is 2 kappa gamma4 (y(T_{x+}-T_{x-}) - x(T_{y+}-T_{y-})).
//
// JS :     Ji decomposition spin angular momentum, (1/2) psidagger Sigma psi
//             It is -2 kappa { 0.5 i sigma12 [(1 - gamma4) T+ - (1 + gamma4) T-] }
//
// JLPure : Chen decomposition orbit angular momentum, (1/i) psidagger r x Dpure psi
// JLJM :   Jaffe-Manohar decomposition orbit angular momentum, (1/i) psidagger r x partial psi, not gauge inv.
//          at Coulomb gauge, JLJM = JLPure
// JPot :   Wak decomposition, Lpot g x (psidagger r x Aphys psi), Note in our definition, g A -> A.
// 
// Note that, psi is scaled, so a further 2 kappa must be divided after this measurement.
//
// REVISION:
//  [07/10/2019 nbale]
//=============================================================================

#ifndef _CMEASUREAMOMENTUMSTOCHASTIC_H_
#define _CMEASUREAMOMENTUMSTOCHASTIC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAMomentumStochastic)

class CLGAPI CMeasureAMomentumStochastic : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureAMomentumStochastic)
public:
    CMeasureAMomentumStochastic()
        : CMeasureStochastic()
        , m_pDeviceXYBufferJL(NULL)
        , m_pDeviceXYBufferJS(NULL)
        , m_pDeviceXYBufferJLPure(NULL)
        , m_pDeviceXYBufferJLJM(NULL)
        , m_pDeviceXYBufferJPot(NULL)
        , m_bExponential(TRUE)
        , m_bNaive(TRUE)

        , m_pDistributionR(NULL)
        , m_pDistributionJL(NULL)
        , m_pDistributionJS(NULL)
        , m_pDistributionJLPure(NULL)
        , m_pDistributionJLJM(NULL)
        , m_pDistributionJPot(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJL(NULL)
        , m_pHostDistributionJS(NULL)
        , m_pHostDistributionJLPure(NULL)
        , m_pHostDistributionJLJM(NULL)
        , m_pHostDistributionJPot(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bShowResult(FALSE)
        , m_bMeasureJLPure(FALSE)
    {
        
    }

    ~CMeasureAMomentumStochastic();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override  { return TRUE; }

protected:
    
    Real* m_pDeviceXYBufferJL;
    Real* m_pDeviceXYBufferJS;
    Real* m_pDeviceXYBufferJLPure;
    Real* m_pDeviceXYBufferJLJM;
    Real* m_pDeviceXYBufferJPot;

    UBOOL m_bExponential;
    UBOOL m_bNaive;

    UINT* m_pDistributionR;
    Real* m_pDistributionJL;
    Real* m_pDistributionJS;
    Real* m_pDistributionJLPure;
    Real* m_pDistributionJLJM;
    Real* m_pDistributionJPot;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJL;
    Real* m_pHostDistributionJS;
    Real* m_pHostDistributionJLPure;
    Real* m_pHostDistributionJLJM;
    Real* m_pHostDistributionJPot;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bShowResult;

public:

    UBOOL m_bMeasureJLPure;
    TArray<UINT> m_lstR;
    TArray<Real> m_lstJL;
    TArray<Real> m_lstJS;
    TArray<Real> m_lstJLPure;
    TArray<Real> m_lstJLJM;
    TArray<Real> m_lstJPot;

    TArray<Real> m_lstJLAll;
    TArray<Real> m_lstJLInner;
    TArray<Real> m_lstJSAll;
    TArray<Real> m_lstJSInner;
    TArray<Real> m_lstJLPureAll;
    TArray<Real> m_lstJLPureInner;
    TArray<Real> m_lstJLJMAll;
    TArray<Real> m_lstJLJMInner;
    TArray<Real> m_lstJPotAll;
    TArray<Real> m_lstJPotInner;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMSTOCHASTIC_H_

//=============================================================================
// END OF FILE
//=============================================================================