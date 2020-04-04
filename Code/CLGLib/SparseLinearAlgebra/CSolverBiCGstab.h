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

__CLG_REGISTER_HELPER_HEADER(CSLASolverBiCGStab)

class CLGAPI CSLASolverBiCGStab : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverBiCGStab)

public:

    CSLASolverBiCGStab();
    ~CSLASolverBiCGStab();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, 
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

    //Old version, for compare
    UBOOL Solve1(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) const;

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