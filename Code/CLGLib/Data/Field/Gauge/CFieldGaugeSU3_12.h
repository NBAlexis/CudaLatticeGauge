//=============================================================================
// FILENAME : CFieldGaugeSU3_12.h
//
// DESCRIPTION:
//   Compact SU(3) gauge field using deviceSU3_12 (12 reals instead of 18).
//   Used as backup buffer (m_pUPrime) in force-gradient integrators.
//   Supports cross-type CopyTo with CFieldGaugeSU3.
//
// REVISION:
//  [05/10/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDGAUGE_SU3_12_H_
#define _CFIELDGAUGE_SU3_12_H_

#include "Data/Field/CFieldGauge.h"
#include "Tools/Math/SU3_12.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3_12)

class CLGAPI CFieldGaugeSU3_12 : public CFieldGauge
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldGaugeSU3_12)
public:
    CFieldGaugeSU3_12();
    ~CFieldGaugeSU3_12();

    EFieldType GetFieldType() const override { return EFT_GaugeSU3_12; }
    UINT MatrixN() const override { return 3; }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void DebugPrintMe() const override;

    void CopyBufferTo(CField* pTarget) const override;
    void CopyTo(CField* U) const override;

    UINT GetDeviceMemorySize() const;
    const void* GetData() const override;
    void* GetData() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void Mul(const CField* other, UBOOL bDaggerLeft, UBOOL bDaggerRight) override;
    void LeftMul(const CField* other, UBOOL bDaggerLeft, UBOOL bDaggerRight) override;
    void Dagger() override;
    void ScalarMultply(Real a) override;
    void ScalarMultply(const CLGComplex& a) override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    cuDoubleComplex Dot(const CField* other) const override;
    DOUBLE GetLength() const override;
    UBOOL ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType uiCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0), void* otherParameter = NULL) override;

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyOriginal(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;
    void ExpMult(Real a, CField* U) const override;
    void ElementNormalize() override;
    void SetOneDirectionUnity(BYTE byDir) override;
    void SetOneDirectionZero(BYTE byDir) override;
    void TransformToIA() override;
    void TA() override;
    void TransformToU() override;
    void CalculateE_Using_U(CFieldGauge* pResoult) const override;
    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive) const override;
    void PolyakovOnSpatialSite(cuDoubleComplex* buffer, BYTE byDir = 3) const override;

    deviceSU3_12* m_pDeviceSU3_12Data;

    //Improve-1: this field's own buffer handle (pool copies each bind theirs).
    CHaloBufferHandle m_HaloBuffer;
    CHaloBufferHandle* GetHaloBufferHandle() override { return &m_HaloBuffer; }
    const CHaloBufferHandle* GetHaloBufferHandle() const override { return &m_HaloBuffer; }

    static void CopySU3ToSU3_12(deviceSU3_12* pDest, const deviceSU3* pSrc, UINT uiCount);
    static void CopySU3_12ToSU3(deviceSU3* pDest, const deviceSU3_12* pSrc, UINT uiCount);
    static DOUBLE SU3_12MSE(const deviceSU3_12* pCompact, const deviceSU3* pRef, UINT uiLinkCount);
};

__END_NAMESPACE

#endif //!_CFIELDGAUGE_SU3_12_H_
