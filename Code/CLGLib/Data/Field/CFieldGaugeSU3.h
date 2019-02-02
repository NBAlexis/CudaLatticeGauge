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
    intokernal; \
    for (UINT it = 0; it < uiTLength; ++it) \
    { \
        coord[3] = it; \
        for (UINT idir = 0; idir < uiDir; ++idir) \
        { \
           UINT uiLinkIndex = _deviceGetLinkIndex(coord, idir); 


#define gaugeSU3KernelFuncionEnd \
        } \
    } 



__BEGIN_NAMESPACE

class CLGAPI CFieldGaugeSU3 : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeSU3)

public:
    CFieldGaugeSU3();
    ~CFieldGaugeSU3();

    virtual void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType);
    virtual void InitialField(EFieldInitialType eInitialType);
    virtual EFieldType GetFieldType() const { return EFT_GaugeSU3; }
    virtual void DebugPrintMe() const;

#pragma region HMC

    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real minusBetaOverN) const;
    virtual void MakeRandomGenerator();
    virtual Real CalculatePlaqutteEnergy(Real minusBetaOverN) const;
    virtual Real CalculatePlaqutteEnergyUsingStable(Real minusBetaOverN, const CFieldGauge *pStable) const;
    virtual Real CalculateKinematicEnergy() const;

#pragma endregion

#pragma region BLAS

    virtual void Zero();
    virtual void Indentity();

    virtual void AxpyPlus(const CField* x);
    virtual void AxpyMinus(const CField* x);
    virtual void Axpy(Real a, const CField* x);
    virtual void Axpy(const _Complex& a, const CField* x);
    virtual void ScalarMultply(const _Complex& a);
    virtual void ScalarMultply(Real a);

#pragma endregion

    virtual void ExpMult(const _Complex& a, CField* U) const;
    virtual void ElementNormalize();
    virtual _Complex Dot(const CField* other) const;

    deviceSU3* m_pDeviceData;

protected:

    void SetByArray(Real* array);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================