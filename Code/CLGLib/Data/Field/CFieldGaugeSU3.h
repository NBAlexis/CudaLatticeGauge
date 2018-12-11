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
public:

    CFieldGaugeSU3(CLatticeData* pLattice, EFieldInitialType eInitialType = EFIT_Identity);
    ~CFieldGaugeSU3()
    {
        checkCudaErrors(cudaFree(m_pDeviceData));
    }

    virtual EFieldType GetFieldType() const { return EFT_GaugeSU3; }

#pragma region HMC

    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, const cuComplex& minusBetaOverN) const;
    virtual void ExpMult(const cuComplex& a, UINT uiPrecision, CField* U) const;
    virtual void CopyTo(CField* U) const;
    virtual void MakeRandomGenerator();

#pragma endregion HMC

#pragma region BLAS

    virtual void Zero();
    virtual void Indentity();
    virtual void Axpy(const CField* x);
    virtual void Axpy(FLOAT a, const CField* x);
    virtual void Axpy(const cuComplex& a, const CField* x);

#pragma endregion BLAS

protected:

    deviceSU3* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3_H_

//=============================================================================
// END OF FILE
//=============================================================================