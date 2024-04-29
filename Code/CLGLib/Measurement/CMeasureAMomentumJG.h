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
// JG               : r x (E x B),   JG is also S1, where S = S0 + OmegaI S1 + OmegaI^2 S2
// JS               : E x Aphys
// JGChen           : E (r x Dpure) Aphys, gauge inv. version
// JGChenApprox     : E (r x Dpure) Aphys, use Dpure a = partial - i g [Apure, a],
//              with i a Dpure A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
// JGChenApprox2    : E (r x Dpure) Aphys, use Dpure a = partial - i g [Apure, a]
//              with i a D A = (A_nu (n+mu) - A_nu (n-mu))/2 + iApure _mu A _nu - i A _nu Apure _mu
//
//  Both JGChenApprox and JGChenApprox2 are not gauge inv.
//
// For Projective plane:
//
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
extern CLGAPI void _AverageXYPlane(Real* pDeviceRes);

/**
* array[x, y] = 0
*/
extern CLGAPI void _ZeroXYPlane(Real* pDeviceRes);

extern CLGAPI void _AverageXYPlaneC(CLGComplex* pDeviceRes);

class CLGAPI CMeasureAMomentumJG : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumJG)
public:
    CMeasureAMomentumJG() 
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistributionJG(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJG(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bMeasureDistribution(FALSE)

        , m_bShowResult(TRUE)

        , m_bMeasureSpin(FALSE)
        , m_bNaiveNabla(FALSE)
        , m_pE(NULL)
        , m_pDpureA(NULL)
        , m_bMeasureApprox(FALSE)
        , m_bProjectivePlane(FALSE)
    {
    }

    ~CMeasureAMomentumJG();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    Real * m_pHostDataBuffer;
    Real * m_pDeviceDataBuffer;

    UINT* m_pDistributionR;
    Real* m_pDistributionJG;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJG;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bMeasureDistribution;

    UBOOL m_bShowResult;
    TArray<Real> m_lstRes;
    TArray<Real> m_lstResJGS;
    TArray<Real> m_lstResJGChen;
    TArray<Real> m_lstResJGChenApprox;
    TArray<Real> m_lstResJGChenApprox2;
    TArray<Real> m_lstResJGSurf;
    TArray<Real> m_lstResJGPot;
    TArray<Real> m_lstResJGS2;

    UBOOL m_bMeasureSpin;
    UBOOL m_bNaiveNabla;
    CFieldGauge* m_pE;
    CFieldGauge* m_pDpureA;

public:

    UBOOL m_bMeasureApprox;
    UBOOL m_bProjectivePlane;
    TArray<Real> m_lstJGAll;
    TArray<Real> m_lstJGInner;

    TArray<UINT> m_lstR;
    //jg as function of R
    TArray<Real> m_lstJG;

    TArray<Real> m_lstJGSAll;
    TArray<Real> m_lstJGSInner;
    //jgs as function of R
    TArray<Real> m_lstJGS;

    TArray<Real> m_lstJGChenAll;
    TArray<Real> m_lstJGChenInner;
    //jg chen as function of R
    TArray<Real> m_lstJGChen;

    TArray<Real> m_lstJGChenApproxAll;
    TArray<Real> m_lstJGChenApproxInner;
    //jg chen as function of R
    TArray<Real> m_lstJGChenApprox;

    TArray<Real> m_lstJGChenApprox2All;
    TArray<Real> m_lstJGChenApprox2Inner;
    //jg chen as function of R
    TArray<Real> m_lstJGChenApprox2;

    TArray<Real> m_lstJGSurfAll;
    TArray<Real> m_lstJGSurfInner;
    //jg surf as function of R
    TArray<Real> m_lstJGSurf;

    TArray<Real> m_lstJGPotAll;
    TArray<Real> m_lstJGPotInner;
    //jg surf as function of R
    TArray<Real> m_lstJGPot;

    TArray<Real> m_lstJGS2All;
    TArray<Real> m_lstJGS2Inner;
    //s2 as function of R
    TArray<Real> m_lstJGS2;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMJG_H_

//=============================================================================
// END OF FILE
//=============================================================================