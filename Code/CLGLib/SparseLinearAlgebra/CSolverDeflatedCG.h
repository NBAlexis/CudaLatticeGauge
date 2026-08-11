//=============================================================================
// FILENAME : CSolverDeflatedCG.h
// 
// DESCRIPTION:
// Deflated conjugate gradient solver for Hermitian positive definite
// operators (EFO_F_DDdagger only).
// A deflation space is built by a short Arnoldi run followed by computing
// the smallest Ritz vectors of the Hessenberg. The solution on the
// deflation space is corrected exactly (coarse grid correction), and the
// rest is solved by plain CG. The deflation space is reused between solves,
// which amortizes the Arnoldi cost inside HMC.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#ifndef _CSOLVERDEFLATEDCG_H_
#define _CSOLVERDEFLATEDCG_H_

#include "CSolverCG.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSLASolverDeflatedCG)

class CLGAPI CSLASolverDeflatedCG : public CSLASolverCG
{
    __CLGDECLARE_CLASS(CSLASolverDeflatedCG)

public:

    CSLASolverDeflatedCG();
    ~CSLASolverDeflatedCG();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;

protected:

    UBOOL SolveHPD(CField* pFieldX, const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, const CField* pStart) override;

    /**
     * Arnoldi run of m_uiDeflateDim steps on pFieldB, then compute the
     * m_uiDeflate smallest Ritz vectors of the Hessenberg as deflation space
     */
    void BuildDeflationSpace(const CField* pFieldB,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM);

    /**
     * Solve E c = U^+ r exactly on the deflation space, then
     * x = x + sum_i c_i U_i, r = r - sum_i c_i (A U_i)
     */
    void CoarseCorrection(CField* pFieldX, CField* pR,
        INT gaugeNum, INT bosonNum, INT tensor2Num,
        const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM);

    UINT m_uiDeflateDim;
    UINT m_uiDeflate;
    UINT m_uiReDeflateInterval;
    UINT m_uiSolveCount;
    UBOOL m_bDeflated;

    //deflation space, kept between solves
    TArray<CField*> m_lstDeflateU;

    class CLinearAlgebraHelper* m_pHelper;

    CLGComplex* m_pDeviceHm;
    CLGComplex* m_pDeviceEigenValue;
    CLGComplex* m_pDevicePk;

    CLGComplex* m_pHostHm;
    CLGComplex* m_pHostE;
    CLGComplex* m_pHostQ;
    CLGComplex* m_pHostR;
    CLGComplex* m_pHostY;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERDEFLATEDCG_H_

//=============================================================================
// END OF FILE
//=============================================================================
