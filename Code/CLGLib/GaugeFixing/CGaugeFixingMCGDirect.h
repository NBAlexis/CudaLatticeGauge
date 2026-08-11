//=============================================================================
// FILENAME : CGaugeFixingMCGDirect.h
//
// DESCRIPTION:
// Direct Maximal Center Gauge (MCG) fixing for SU(3).
//
// Standard algorithm (Montero, hep-lat/9906010):
// Directly maximizes the MCG functional
//   R = (1/(N_site * N_dim * N^2)) * sum_{x,\mu} |Tr U_\mu(x)|^2
// by local gauge updates using the Cabibbo-Marinari-Okawa method.
//
// The local quantity to maximize at each site x is:
//   R_x = sum_\mu |Tr{G(x) U_\mu(x)}|^2 + |Tr{U_\mu(x-\mu) G^\dagger(x)}|^2
//
// The optimal local gauge transform G(x) is found by projecting the
// trace-weighted staple A(x)^\dagger to SU(3) using Cabibbo-Marinari
// subgroup decomposition.
//
// Overrelaxation: G_\omega = (1-\omega) I + \omega G, then re-project.
// Red-black checkerboard sweep for parallelization.
//
// Convergence criterion:
//   \theta = sum_{x,\mu} (N^2 - |Tr U_\mu(x)|^2) / (N^2 * V * D)
//   Stop when \theta < accuracy.
//
// References:
//  hep-lat/9906010 -- Montero, "Study of SU(3) vortex-like configurations
//                     with a new maximal center gauge fixing method"
//  hep-lat/9708008 -- Brower et al., "Center vortices in SU(2) lattice gauge fields"
//  hep-lat/0003021 -- Langfeld et al., "SU(N) vortices and Wilson loops"
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGMCGDIRECT_H_
#define _CGAUGEFIXINGMCGDIRECT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingMCGDirect)

class CLGAPI CGaugeFixingMCGDirect : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingMCGDirect)
public:

    CGaugeFixingMCGDirect()
    : CGaugeFixing()
    , m_fOmega(F(1.5))
    , m_iCheckErrorStep(1000)
    {
    }

    ~CGaugeFixingMCGDirect()
    {
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;
#else
    Real CheckRes(const CFieldGauge* pGauge) override;
#endif

    //P4-2.5: core deviation loop over one gauge buffer (local or gathered-global).
    //Two overloads because SU2/SU3 use different kernels and normalization.
    DOUBLE CheckResLocalSU2(const deviceSU2* pGaugeData, BYTE byFieldId);
    DOUBLE CheckResLocalSU3(const deviceSU3* pGaugeData, BYTE byFieldId);

    //P4-2.5: the Direct-MCG iteration loop over one gauge buffer (local or
    //gathered-global). Public so CGaugeFixingMCGIndirect can run Stage 2 inside
    //its own single gather without re-entering this fixer's MG path.
    void GaugeFixingLoopSU2(deviceSU2* pDeviceBufferPointer, BYTE byFieldId);
    void GaugeFixingLoopSU3(deviceSU3* pDeviceBufferPointer, BYTE byFieldId);

    CCString GetInfos(const CCString& sTab) const override;

    /**
     * Project an SU(2) gauge field to Z_2 center elements.
     * For each link U, find the nearest Z_2 element:
     *   Z_2 = { +I, -I }
     *
     * Algorithm (Del Debbio et al., hep-lat/9610005):
     *   For SU(2), Tr(U) is real. Project to sign(Tr U) * I.
     */
    static void CenterProjection(CFieldGaugeSU2* pGauge);

    /**
     * Project an SU(3) gauge field to Z_3 center elements.
     * For each link U, find the nearest Z_3 element:
     *   Z_3 = { exp(2*pi*i*m/3) * I | m = 0, 1, 2 }
     *
     * Algorithm (Tucker & Stack, hep-lat/0110165):
     *   Compute theta = arg(Tr U).
     *   m = 0  if |theta| <= pi/3      -> nearest to I
     *   m = 1  if  pi/3 < theta <= pi  -> nearest to exp(2*pi*i/3) I
     *   m = 2  if -pi < theta < -pi/3  -> nearest to exp(-2*pi*i/3) I
     */
    static void CenterProjection(CFieldGaugeSU3* pGauge);

    /**
     * Remove the Z_2 center from an SU(2) gauge field, in place.
     * For each link U, find the nearest Z_2 element Z = sign(Tr U) * I
     * (Z^dag = Z for Z_2) and set
     *   U <- Z * U,
     * which projects out the center phase while keeping the coset part
     * SU(2)/Z_2. The resulting links stay in SU(2).
     *
     * This is the complement of CenterProjection: it keeps the "unprojected"
     * (vortex) part of the links for center-vortex analysis
     * (e.g. Greensite, "An Introduction to the Confinement Problem").
     */
    static void CenterRemove(CFieldGaugeSU2* pGauge);

    /**
     * Remove the Z_3 center from an SU(3) gauge field, in place.
     * For each link U, find the nearest Z_3 element Z = exp(2*pi*i*m/3) * I
     * (same m selection as CenterProjection) and set
     *   U <- Z^dag * U,
     * which projects out the center phase while keeping the coset part
     * SU(3)/Z_3. The resulting links stay in SU(3).
     *
     * This is the complement of CenterProjection: it keeps the "unprojected"
     * (vortex) part of the links for center-vortex analysis.
     */
    static void CenterRemove(CFieldGaugeSU3* pGauge);

    Real m_fOmega;
    UINT m_iCheckErrorStep;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGMCGDIRECT_H_

//=============================================================================
// END OF FILE
//=============================================================================
