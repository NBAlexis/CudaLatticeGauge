//=============================================================================
// FILENAME : CMultiShiftNested.h
// 
// DESCRIPTION:
// When the light quark is mass preconditionered, it is often (x+small/x)^{1/2}
// And usually, only one order is accurate enough
// (x+small/x)^{1/2} = c + a/x with only ONE inverse
// It is no need to use multi-shift solver
// So we have a nested multi-shift solver to mimic this
//
// REVISION:
//  [10/07/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTNESTED_H_
#define _CMULTISHIFTNESTED_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMultiShiftNested)

class CLGAPI CMultiShiftNested : public CMultiShiftSolver
{
    __CLGDECLARE_CLASS(CMultiShiftNested)

public:

    enum { _kMaxStep = 30, };

    CMultiShiftNested();
    ~CMultiShiftNested();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    UBOOL Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, 
        INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    CSLASolver* m_pNestedSolver;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTNESTED_H_

//=============================================================================
// END OF FILE
//=============================================================================