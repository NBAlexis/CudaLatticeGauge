//=============================================================================
// FILENAME : CGaugeFixingMAG.h
//
// DESCRIPTION:
// Maximal Abelian Gauge (MAG) fixing for SU(2) and SU(3).
//
// SU(2) MAG algorithm (Cea & Cosmai, hep-lat/9504008, Appendix A):
// The MAG maximizes the functional
//   G_MAG = sum_{x,mu} Tr(U_mu^dag sigma_3 U_mu sigma_3)
// which is equivalent to maximizing the sum of squared diagonal elements.
//
// Local update constructs the traceless Hermitian matrix:
//   X(x) = sum_mu [ U_mu(x) sigma_3 U_mu^dag(x)
//                 + U_mu^dag(x-mu) sigma_3 U_mu(x-mu) ]
// Then V(x) = X(x) sigma_3 / k(x), k(x) = sqrt(det(X sigma_3)), V in SU(2).
// Writing V = v0 I + i(v1 sigma_1 + v2 sigma_2), the optimal gauge
// transformation with overrelaxation parameter omega is:
//   g^omega = cos(omega * alpha) I
//           - i(v1 sigma_1 + v2 sigma_2) / sqrt(1-v0^2) * sin(omega * alpha)
// where cos(2*alpha) = v0, i.e. alpha = acos(v0) / 2.
//
// SU(3) MAG algorithm (Tucker & Stack, hep-lat/0110165):
// Uses Cabibbo-Marinari on SU(2) subgroups (to be implemented).
//
// Red-black checkerboard sweep for parallelization.
// Convergence criterion: sum of squared off-diagonal elements.
//
// References:
//  hep-lat/9504008 -- Cea & Cosmai, "Maximal Abelian Gauge in SU(2)"
//  hep-lat/0110165 -- Tucker & Stack, "The Maximal Abelian Gauge in SU(3)"
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGMAG_H_
#define _CGAUGEFIXINGMAG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingMAG)

class CLGAPI CGaugeFixingMAG : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingMAG)
public:

    CGaugeFixingMAG()
    : CGaugeFixing()
    , m_fOmega(F(1.5))
    , m_iCheckErrorStep(1000)
    {
    }

    ~CGaugeFixingMAG()
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
    DOUBLE CheckResLocalSU2(const deviceSU2* pGaugeData, BYTE byFieldId);
    DOUBLE CheckResLocalSU3(const deviceSU3* pGaugeData, BYTE byFieldId);

    //P4-2.5: the MAG iteration loop over one gauge buffer (local or
    //gathered-global). Public so CGaugeFixingMCGIndirect can run Stage 1 inside
    //its own single gather without re-entering this fixer's MG path.
    void GaugeFixingLoopSU2(deviceSU2* pDeviceBufferPointer, BYTE byFieldId);
    void GaugeFixingLoopSU3(deviceSU3* pDeviceBufferPointer, BYTE byFieldId);

    CCString GetInfos(const CCString& sTab) const override;

    /**
     * Project an SU(2) gauge field to the U(1) abelian subgroup.
     *
     * For each SU(2) link U = [a, b; -b*, a*], keep the diagonal phase
     * and zero the off-diagonal:
     *   U_proj = [a/|a|, 0; 0, conj(a)/|a|]
     */
    static void MaximalAbelianProjection(CFieldGaugeSU2* pGauge);

    /**
     * Project an SU(3) gauge field to the U(1)xU(1) abelian subgroup.
     *
     * Algorithm (Tucker & Stack, hep-lat/0110165):
     *   For each SU(3) link matrix A, keep diagonal elements, zero
     *   off-diagonal elements, then adjust overall phase to enforce det=1.
     */
    static void MaximalAbelianProjection(CFieldGaugeSU3* pGauge);

    Real m_fOmega;
    UINT m_iCheckErrorStep;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGMAG_H_

//=============================================================================
// END OF FILE
//=============================================================================
