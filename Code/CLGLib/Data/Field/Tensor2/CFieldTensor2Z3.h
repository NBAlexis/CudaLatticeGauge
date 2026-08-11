//=============================================================================
// FILENAME : CFieldTensor2Z3.h
//
// DESCRIPTION:
// The Z3 (center) tensor2 (plaquette) field for the PSU(3) fundamental-lift
// action. Every site has 6 elements (xy, xz, xt, yz, yt, zt), each element is
// a deviceZN<3> (a 3rd root of unity), stored component-major:
//   data[plaqutteIndex * siteCount + siteIndex]
// (the same layout as CFieldTensor2T / the Fmunu buffers in CStapleCache).
//
// Unlike the generic tensor2 fields, this class enforces the discrete Z3
// constraint:
//   - EFIT_Zero / EFIT_Identity both set all elements to root 0 (the group
//     identity), never to complex 0;
//   - EFIT_Random draws each element from the three Z3 roots (it does NOT
//     enforce flatness dB=0; the action checks flatness when required);
//   - EFIT_RandomGaussian / EFIT_RandomGenerator are rejected with a fatal
//     error (the generic deviceZN initializer silently returns the identity,
//     which must not leak through this interface);
//   - file / byte loading canonicalizes every element to the nearest Z3 root
//     and rejects values farther than the tolerance instead of silently
//     projecting arbitrary complex numbers.
//
// The linear operations (Axpy*, ScalarMultply, Dagger, Mul, LeftMul) are not
// closed on the Z3 set and are rejected explicitly.
//
// REVISION:
//  [08/09/26]
//=============================================================================
#pragma once

#ifndef _CFIELDTENSOR2Z3_H_
#define _CFIELDTENSOR2Z3_H_

#include "CFieldTensor2T.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldTensor2Z3)

class CLGAPI CFieldTensor2Z3 : public CFieldTensor2T<deviceZN<3>>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldTensor2Z3)
public:

    EFieldType GetFieldType() const override { return EFT_Tensor2Z3; }

    /**
    * 6 elements per site, each deviceZN<3> stores one complex number:
    * 6 * 2 = 12 reals per site.
    */
    UINT FloatN() const override { return 12; }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialWithByte(BYTE* byData) override;

    /**
    * Count the Z3 monopole defects of the field.
    * For every cube (x, mu<nu<rho) the charge is
    *   m = b_nu_rho(x+mu) - b_nu_rho(x) - b_mu_rho(x+nu) + b_mu_rho(x)
    *     + b_mu_nu(x+rho) - b_mu_nu(x)   (mod 3)
    * with the four orientations xyz, xyt, xzt, yzt. m != 0 is a monopole /
    * dislocation. A flat field (dB = 0) has zero monopoles; the unrestricted
    * Villain ensemble (AllowMonopole=1) generically has non-zero ones.
    */
    UINT CountMonopoles() const;

    //=====================================================================
    // The Z3 set {1, zeta, zeta^2} is not closed under these operations.
    // Reject them loudly instead of silently producing a non-Z3 value.
    //=====================================================================
    void AxpyPlus(const CField* x) override
    {
        appCrucial(_T("CFieldTensor2Z3: AxpyPlus is not defined on the Z3 set (not closed under addition).\n"));
    }

    void AxpyMinus(const CField* x) override
    {
        appCrucial(_T("CFieldTensor2Z3: AxpyMinus is not defined on the Z3 set (not closed under subtraction).\n"));
    }

    void Axpy(Real a, const CField* x) override
    {
        appCrucial(_T("CFieldTensor2Z3: Axpy(Real) is not defined on the Z3 set.\n"));
    }

    void Axpy(const CLGComplex& a, const CField* x) override
    {
        appCrucial(_T("CFieldTensor2Z3: Axpy(Complex) is not defined on the Z3 set.\n"));
    }

    void ScalarMultply(const CLGComplex& a) override
    {
        appCrucial(_T("CFieldTensor2Z3: ScalarMultply(Complex) is not defined on the Z3 set.\n"));
    }

    void ScalarMultply(Real a) override
    {
        appCrucial(_T("CFieldTensor2Z3: ScalarMultply(Real) is not defined on the Z3 set.\n"));
    }

    void Dagger() override
    {
        appCrucial(_T("CFieldTensor2Z3: Dagger is not supported (write back via SetFromIndex if ever needed).\n"));
    }

    void Mul(const CField* x, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE) override
    {
        appCrucial(_T("CFieldTensor2Z3: Mul is not supported (use the action kernels for Z3 updates).\n"));
    }

    void LeftMul(const CField* x, UBOOL bDaggerLeft = FALSE, UBOOL bDaggerRight = FALSE) override
    {
        appCrucial(_T("CFieldTensor2Z3: LeftMul is not supported (use the action kernels for Z3 updates).\n"));
    }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDTENSOR2Z3_H_

//=============================================================================
// END OF FILE
//=============================================================================
