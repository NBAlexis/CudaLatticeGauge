//=============================================================================
// FILENAME : CSolverBiCGStab.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [01/08/2019 nbale]
//=============================================================================

#ifndef _CSOLVERBICGSTAB_H_
#define _CSOLVERBICGSTAB_H_

__BEGIN_NAMESPACE

class CLGAPI CSLASolverBiCGStab : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverBiCGStab)

public:

    CSLASolverBiCGStab();
    ~CSLASolverBiCGStab();
    
    virtual void Configurate(const CParameters& param);
    virtual void AllocateBuffers(const CField* pField);
    virtual void ReleaseBuffers();
    virtual UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart = NULL);

protected:
   
    UINT m_uiReTry;
    UINT m_uiDevationCheck;
    UINT m_uiStepCount;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERBICGSTAB_H_

//=============================================================================
// END OF FILE
//=============================================================================