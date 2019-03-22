//=============================================================================
// FILENAME : CSolverGMRES.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [02/12/2019 nbale]
//=============================================================================

#ifndef _CSOLVERGMRES_H_
#define _CSOLVERGMRES_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSLASolverGMRES)

class CLGAPI CSLASolverGMRES : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverGMRES)

public:

    enum { _kMaxStep = 100, };

    CSLASolverGMRES();
    ~CSLASolverGMRES();

    virtual void Configurate(const CParameters& param);
    virtual void AllocateBuffers(const CField* pField);
    virtual void ReleaseBuffers();
    virtual UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart = NULL);

protected:

    void RotateH(/*UINT uiK*/);
    void SolveY(/*UINT uiHeisenbergDim*/);
    inline UINT HIndex(UINT x /*0-k*/, UINT y /*0-(k-1)*/) { return y + x * m_uiMaxDim; }

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    Real m_fAccuracy;
    Real m_fBeta;
    UBOOL m_bUseCudaForSmallMatrix;

    CLGComplex m_h[(_kMaxStep + 1) * _kMaxStep];
    CLGComplex m_y[_kMaxStep];
    CLGComplex m_g[_kMaxStep + 1];

    TArray<class CField*> m_lstVectors;

    class CLinearAlgebraHelper* m_pHelper;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGMRES_H_

//=============================================================================
// END OF FILE
//=============================================================================