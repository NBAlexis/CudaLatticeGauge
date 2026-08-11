//=============================================================================
// FILENAME : CMultiShiftCG.h
// 
// DESCRIPTION:
// This is the class for multi-shift conjugate gradient solver.
// It only works with Hermitian positive definite operators, so only
// EFO_F_DDdagger is supported, which is exactly the operator of RHMC.
// See: B. Jegerlehner, hep-lat/9608029
//
// REVISION:
//  [mm/dd/yy]
//  [07/18/2026 nbale]
//=============================================================================

#ifndef _CMULTISHIFTCG_H_
#define _CMULTISHIFTCG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMultiShiftCG)

class CLGAPI CMultiShiftCG : public CMultiShiftSolver
{
    __CLGDECLARE_CLASS(CMultiShiftCG)

public:

    CMultiShiftCG();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    UBOOL Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    UINT m_uiMaxStep;

    TArray<class CField*> m_lstVectors;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTCG_H_

//=============================================================================
// END OF FILE
//=============================================================================
