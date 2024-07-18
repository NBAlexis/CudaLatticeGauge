//=============================================================================
// FILENAME : CFieldGaugeLinkDirichlet.h
// 
// DESCRIPTION:
// There is simplifications for periodic boundary condition
// which is invalid for Dirichlet.
// This is only for Dirichlet.
//
// REVISION:
//  [07/06/2019 nbale]
//=============================================================================
#include "CFieldGaugeLink.h"

#define __DEFINE_GAUGE_LINKD(N) \
__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU##N##D) \
class CLGAPI CFieldGaugeSU##N##D : public CFieldGaugeLinkD<deviceSU##N, N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU##N##D) \
public: \
    EFieldType GetFieldType() const override { return EFT_GaugeSU##N; } \
};

#ifndef _CFIELDGAUGELINKDIRICHLET_H_
#define _CFIELDGAUGELINKDIRICHLET_H_

__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeLinkD : public CFieldGaugeLink<deviceGauge, matrixN>
{
public:

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;

    void MakeRandomGenerator() override;
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculateKinematicEnergy() const override;

    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }

#pragma endregion

#pragma region BLAS

    void FixBoundary() override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;
    CCString GetInfos(const CCString &tab) const override;

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

#pragma endregion

};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1D)

class CLGAPI CFieldGaugeU1D : public CFieldGaugeLinkD<CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeU1D)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
    void TransformToIA() override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU2D)

class CLGAPI CFieldGaugeSU2D : public CFieldGaugeLinkD<deviceSU2, 2>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU2D)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }

    void InitialWithByteCompressed(const CCString& sFileName) override;
    CCString SaveToCompressedFile(const CCString& fileName) const override;
    void TransformToIA() override;
};

__DEFINE_GAUGE_LINKD(4)
__DEFINE_GAUGE_LINKD(5)
__DEFINE_GAUGE_LINKD(6)
__DEFINE_GAUGE_LINKD(7)
__DEFINE_GAUGE_LINKD(8)


__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGELINKDIRICHLET_H_

//=============================================================================
// END OF FILE
//=============================================================================