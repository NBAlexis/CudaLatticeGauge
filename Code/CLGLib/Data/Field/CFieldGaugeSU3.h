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

    virtual void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType);
    virtual void InitialWithByte(BYTE* byData);
    virtual void InitialField(EFieldInitialType eInitialType);
    virtual EFieldType GetFieldType() const { return EFT_GaugeSU3; }
    virtual void DebugPrintMe() const;

#pragma region HMC

    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const;
    virtual void CalculateOnlyStaple(CFieldGauge* pStaple) const;
    virtual void MakeRandomGenerator();
    virtual Real CalculatePlaqutteEnergy(Real betaOverN) const;
    virtual Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const;
    virtual Real CalculateKinematicEnergy() const;

#pragma endregion

#pragma region BLAS

    virtual void Zero();
    virtual void Indentity();

    virtual void AxpyPlus(const CField* x);
    virtual void AxpyMinus(const CField* x);
    virtual void Axpy(Real a, const CField* x);
    virtual void Axpy(const CLGComplex& a, const CField* x);
    virtual void ScalarMultply(const CLGComplex& a);
    virtual void ScalarMultply(Real a);

#pragma endregion

    virtual void ExpMult(Real a, CField* U) const;

    virtual void ElementNormalize();
    virtual CLGComplex Dot(const CField* other) const;
    virtual void SaveToFile(const CCString &fileName) const;
    virtual BYTE* CopyDataOut(UINT &uiSize) const;
    virtual CCString GetInfos(const CCString &tab) const;

    deviceSU3* m_pDeviceData;

protected:

    void SetByArray(Real* array);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================