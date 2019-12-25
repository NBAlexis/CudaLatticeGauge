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
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
    void DebugPrintMe() const override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override;
    Real CalculateKinematicEnergy() const override;

#pragma endregion

#pragma region BLAS

    void Zero() override;
    void Identity() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;

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

#pragma endregion

    void ExpMult(Real a, CField* U) const override;

    void ElementNormalize() override;
    CLGComplex Dot(const CField* other) const override;
    void SaveToFile(const CCString &fileName) const override;
    BYTE* CopyDataOut(UINT &uiSize) const override;
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