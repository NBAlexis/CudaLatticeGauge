//=============================================================================
// FILENAME : CMultiShiftBiCGStab.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
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
        INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    UINT m_uiDevationCheck;
    UINT m_uiStepCount;
#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fAccuracy;
#else
    Real m_fAccuracy;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTBICGSTAB_H_

//=============================================================================
// END OF FILE
//=============================================================================