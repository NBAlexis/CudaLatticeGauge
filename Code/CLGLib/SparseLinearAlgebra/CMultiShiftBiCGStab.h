//=============================================================================
// FILENAME : CMultiShiftBiCGStab.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [mm/dd/yy]
//  [20/06/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTBICGSTAB_H_
#define _CMULTISHIFTBICGSTAB_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMultiShiftBiCGStab)

class CLGAPI CMultiShiftBiCGStab : public CMultiShiftSolver
{
    __CLGDECLARE_CLASS(CMultiShiftBiCGStab)

public:

    enum { _kMaxStep = 100, };

    CMultiShiftBiCGStab();
    ~CMultiShiftBiCGStab();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, 
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    UINT m_uiDevationCheck;
    UINT m_uiStepCount;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTBICGSTAB_H_

//=============================================================================
// END OF FILE
//=============================================================================