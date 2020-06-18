//=============================================================================
// FILENAME : CMultiShiftGMRES.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [15/06/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTGMRES_H_
#define _CMULTISHIFTGMRES_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMultiShiftGMRES)

class CLGAPI CMultiShiftGMRES : public CMultiShiftSolver
{
    __CLGDECLARE_CLASS(CMultiShiftGMRES)

public:

    enum { _kMaxStep = 30, };

    CMultiShiftGMRES();
    ~CMultiShiftGMRES();

    void Configurate(const CParameters& param) override;
    void AllocateBuffers(const CField* pField) override;
    virtual void ReleaseBuffers();
    UBOOL Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, const CFieldGauge* pGaugeFeild,
        EFieldOperator uiM, ESolverPhase ePhase = ESP_Once, const CField* pStart = NULL) override;

protected:

    inline UINT HIndex(UINT x /*0-k*/, UINT y /*0-(k-1)*/) const { return y + x * m_uiMaxDim; }

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    Real m_fAccuracy;
    Real m_fBeta;
    UBOOL m_bUseCudaForSmallMatrix;

    CLGComplex m_h[(_kMaxStep + 1) * _kMaxStep];
    CLGComplex m_hCopy[(_kMaxStep + 1) * (_kMaxStep + 1)];
    CLGComplex m_y[_kMaxStep];
    CLGComplex m_yhat[_kMaxStep + 1];
    CLGComplex m_g[_kMaxStep + 1];
    CLGComplex m_z[_kMaxStep + 1];

    TArray<class CField*> m_lstVectors;

    class CLinearAlgebraHelper* m_pHelper;

    CLGComplex* m_pDeviceY;
    CLGComplex* m_pDeviceZ;
    CLGComplex* m_pDeviceHHat;
    CLGComplex* m_pDeviceQ;
    CLGComplex* m_pDeviceR;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTGMRES_H_

//=============================================================================
// END OF FILE
//=============================================================================