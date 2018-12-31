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
        { 
            


#define gaugeSU3KernelFuncionEnd \
        } \
    } 



__BEGIN_NAMESPACE

class CLGAPI CFieldGaugeSU3 : public CFieldGauge
{
    __CLGDECLARE_CLASS(CFieldGaugeSU3)

public:
    CFieldGaugeSU3();
    ~CFieldGaugeSU3();

    virtual void InitialField(EFieldInitialType eInitialType);
    virtual EFieldType GetFieldType() const { return EFT_GaugeSU3; }
    virtual void DebugPrintMe() const;

#pragma region HMC

    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, const _Complex& minusBetaOverN) const;
    virtual void ExpMult(const _Complex& a, CField* U) const;
    virtual void CopyTo(CField* U) const;
    virtual void MakeRandomGenerator();
    virtual Real CalculatePlaqutteEnergy(const _Complex& minusBetaOverN);

#pragma endregion HMC

#pragma region BLAS

    virtual void Zero();
    virtual void Indentity();
    virtual void Axpy(const CField* x);
    virtual void Axpy(Real a, const CField* x);
    virtual void Axpy(const _Complex& a, const CField* x);

#pragma endregion BLAS

    deviceSU3* m_pDeviceData;

protected:

    Real* m_pDeviceTmpResPtr;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================