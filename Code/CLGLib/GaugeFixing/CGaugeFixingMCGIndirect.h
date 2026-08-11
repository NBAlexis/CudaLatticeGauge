//=============================================================================
// FILENAME : CGaugeFixingMCGIndirect.h
//
// DESCRIPTION:
// Indirect Maximal Center Gauge (MCG) fixing for SU(3).
//
// Two-step algorithm (Brower et al., hep-lat/9708008; Langfeld et al., hep-lat/0003021):
//   Step 1: Fix to Maximal Abelian Gauge (MAG).
//           This makes the gauge field as diagonal as possible.
//   Step 2 (default): use the MAG-fixed field as the starting point for a
//           full Direct MCG solve. This is MAG-preconditioned Direct MCG.
//   Step 2 (UseStandardIMCG=1): restrict the second stage to the residual
//           Abelian subgroup left by MAG.
//
// The indirect approach leverages the fact that after MAG fixing,
// the gauge field is approximately diagonal, making the subsequent
// MCG fixing more effective.
//
// Each step uses the Cabibbo-Marinari-Okawa method with red-black
// checkerboard sweeps and overrelaxation.
//
// Convergence criterion:
//   Step 1: theta_MAG = sum off-diagonal / (V * D) < accuracy
//   Step 2: theta_MCG = sum (9 - |Tr U|^2) / (9 * V * D) < accuracy
//
// MCG is a gauge transformation -- it preserves observables.
// Center projection is NOT included here; it is a separate operation that
// may be performed AFTER MCG gauge fixing if desired.
//
// References:
//  hep-lat/9708008 -- Brower et al., "Center vortices in SU(2) lattice gauge fields"
//  hep-lat/0003021 -- Langfeld et al., "SU(N) vortices and Wilson loops"
//  hep-lat/9906010 -- Montero, "Study of SU(3) vortex-like configurations"
//  hep-lat/0110165 -- Tucker & Stack, "The Maximal Abelian Gauge in SU(3)"
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGMCGINDIRECT_H_
#define _CGAUGEFIXINGMCGINDIRECT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingMCGIndirect)

class CLGAPI CGaugeFixingMCGIndirect : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingMCGIndirect)
public:

    CGaugeFixingMCGIndirect()
    : CGaugeFixing()
    , m_fOmega(F(1.5))
    , m_fOmegaStage1(F(1.5))
    , m_iCheckErrorStep(1000)
    , m_iStage1MaxIterate(100000)
    , m_bUseStandardIMCG(FALSE)
    , m_iIMCGGrid(24)
    {
    }

    ~CGaugeFixingMCGIndirect()
    {
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;
#else
    Real CheckRes(const CFieldGauge* pGauge) override;
#endif

    CCString GetInfos(const CCString& sTab) const override;

    Real m_fOmega;          // Overrelaxation parameter for Stage 2 (MCG)
    Real m_fOmegaStage1;    // Overrelaxation parameter for Stage 1 (MAG)
    UINT m_iCheckErrorStep;
    UINT m_iStage1MaxIterate;  // Max iterations for Stage 1 (MAG)

    // FALSE keeps the historical/default behavior: MAG-preconditioned full
    // Direct MCG. TRUE uses the standard indirect center gauge idea by
    // restricting the second stage to the residual Abelian subgroup.
    UBOOL m_bUseStandardIMCG;

    // Grid resolution for Standard IMCG phase search (default 24).
    // Larger values give finer phase resolution but increase compute cost
    // as O(grid^2) for SU(3).
    INT m_iIMCGGrid;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGMCGINDIRECT_H_

//=============================================================================
// END OF FILE
//=============================================================================
