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

    class CLatticeData* m_pOwner;

protected:

    void RotateH(UINT uiK);
    void SolveY(UINT uiHeisenbergDim);
    static inline UINT HIndex(UINT x /*0-k*/, UINT y /*0-(k-1)*/) { return y + x * _kMaxStep; }

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    UBOOL m_bAbsoluteAccuracy;
    Real m_fAccuracy;
    Real m_fBeta;

    _Complex m_h[(_kMaxStep + 1) * _kMaxStep];
    _Complex m_y[_kMaxStep];
    _Complex m_g[_kMaxStep + 1];

    TArray<class CField*> m_lstVectors;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGMRES_H_

//=============================================================================
// END OF FILE
//=============================================================================