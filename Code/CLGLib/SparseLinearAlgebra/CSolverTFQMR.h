//=============================================================================
// FILENAME : CSolverTFQMR.h
// 
// DESCRIPTION:
// TFQMR is similar as BiCGStab, but a little bit slower
//
// REVISION:
//  [20/06/2020 nbale]
//=============================================================================

#ifndef _CSOLVERTFQMR_H_
#define _CSOLVERTFQMR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSolverTFQMR)

class CLGAPI CSolverTFQMR : public CSLASolver
{
    __CLGDECLARE_CLASS(CSolverTFQMR)

public:

    CSolverTFQMR();
    ~CSolverTFQMR();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, 
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:
   
    UINT m_uiReTry;
    UINT m_uiDevationCheck;
    UINT m_uiStepCount;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERTFQMR_H_

//=============================================================================
// END OF FILE
//=============================================================================