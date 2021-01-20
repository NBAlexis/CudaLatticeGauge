//=============================================================================
// FILENAME : CFieldGaugeSU3.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU3_H_
#define _CFIELDGAUGE_SU3_H_

#define gaugeSU3KernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeSU3KernelFuncionEnd \
    } 



__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3)

class CLGAPI CFieldGaugeSU3 : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeSU3)

public:
    CFieldGaugeSU3();
    ~CFieldGaugeSU3();

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(BYTE*) override;
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
    void DebugPrintMe() const override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;
#else
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUseClover(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override;
    Real CalculateKinematicEnergy() const override;
#endif

#pragma endregion

#pragma region BLAS

    void Zero() override;
    void Identity() override;
    void Dagger() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;

    void SetOneDirectionUnity(BYTE byDir) override;
    void SetOneDirectionZero(BYTE byDir) override;

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

    void CalculateE_Using_U(CFieldGauge* pResoult) const override;

    void CalculateNablaE_Using_U(CFieldGauge* pResoult) const override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;

    void ElementNormalize() override;
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex Dot(const CField* other) const override;
#else
    CLGComplex Dot(const CField* other) const override;
#endif
    CCString SaveToCompressedFile(const CCString& fileName) const override;
    BYTE* CopyDataOut(UINT &uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    CCString GetInfos(const CCString &tab) const override;

    deviceSU3* m_pDeviceData;

protected:

    void SetByArray(Real* array);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================