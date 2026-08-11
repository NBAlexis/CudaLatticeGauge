//=============================================================================
// FILENAME : CMultiShiftFOM.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [mm/dd/yy]
//  [18/06/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTFOM_H_
#define _CMULTISHIFTFOM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMultiShiftFOM)

class CLGAPI CMultiShiftFOM : public CMultiShiftSolver
{
    __CLGDECLARE_CLASS(CMultiShiftFOM)

public:

    enum { _kMaxStep = 100, };

    CMultiShiftFOM();
    ~CMultiShiftFOM();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, 
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    void RotateHSolveY(CLGComplex* h, CLGComplex* g, CLGComplex* y) const;

    static void RotateH(CLGComplex* h, CLGComplex* g, UINT uiDimX, UINT uiDimY);
    static void SolveY(CLGComplex* h, CLGComplex* g, CLGComplex* y, UINT uiDimY, UINT uiDimX);
    inline UINT HIndex(UINT x /*0-k*/, UINT y /*0-(k-1)*/) const { return y + x * m_uiMaxDim; }

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    UBOOL m_bUseCudaForSmallMatrix;

    CLGComplex m_h[(_kMaxStep + 1) * _kMaxStep];
    CLGComplex m_hCopy[_kMaxStep * _kMaxStep];
    CLGComplex m_y[_kMaxStep];
    CLGComplex m_g[_kMaxStep];

    TArray<class CField*> m_lstVectors;

    class CLinearAlgebraHelper* m_pHelper;

    CLGComplex* m_pDeviceY;
    CLGComplex* m_pDeviceHHat;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTFOM_H_

//=============================================================================
// END OF FILE
//=============================================================================