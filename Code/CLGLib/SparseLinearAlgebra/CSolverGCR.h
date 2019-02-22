//=============================================================================
// FILENAME : CSolverGCR.h
// 
// DESCRIPTION:
//
// Both GCR and ORTHOMIN has oscillation problem, ORTHODIR is good
// So the implementation is in fact ORTHODIR
//
// REVISION:
//  [02/16/2019 nbale]
//=============================================================================

#ifndef _CSOLVERGCR_H_
#define _CSOLVERGCR_H_

__BEGIN_NAMESPACE

class CLGAPI CSLASolverGCR : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverGCR)

public:

    enum { _kMaxStep = 100, };

    CSLASolverGCR();
    ~CSLASolverGCR();

    virtual void Configurate(const CParameters& param);
    virtual void AllocateBuffers(const CField* pField);
    virtual void ReleaseBuffers();
    virtual UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart = NULL);

protected:

    UINT m_uiReStart;
    UINT m_uiMaxDim;
    UINT m_uiIterateNumber;
    UINT m_uiCheckError;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGCR_H_

//=============================================================================
// END OF FILE
//=============================================================================