//=============================================================================
// FILENAME : CSolverCG.h
// 
// DESCRIPTION:
// This is the conjugate gradient solver for Hermitian positive definite
// operators. It only works with EFO_F_DDdagger, which is exactly the
// operator solved in HMC.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#ifndef _CSOLVERCG_H_
#define _CSOLVERCG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSLASolverCG)

class CLGAPI CSLASolverCG : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverCG)

public:

    CSLASolverCG();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    UBOOL Solve(CField* pFieldX, const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    /**
     * The CG iteration loop
     * pFieldX and pR (b - A pFieldX) should already be set before calling
     */
    UBOOL IterateCG(CField* pFieldX, CField* pR, const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM);

    /**
     * Solve a Hermitian positive definite operator by plain CG
     */
    virtual UBOOL SolveHPD(CField* pFieldX, const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, const CField* pStart);

    UINT m_uiMaxStep;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERCG_H_

//=============================================================================
// END OF FILE
//=============================================================================
