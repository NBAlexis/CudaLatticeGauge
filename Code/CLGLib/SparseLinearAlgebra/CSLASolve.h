//=============================================================================
// FILENAME : CSLASolver.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [1/3/2019 nbale]
//=============================================================================

#ifndef _CSLASOLVER_H_
#define _CSLASOLVER_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(ESLASolverType,
    ESST_BiCGstab,
    ESST_Max,

    ESST_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CSLASolver : public CBase
{
public:

    CSLASolver() : m_pOwner(NULL) { ; }

    virtual void Configurate(const CParameters& param) = 0;

    /**
    * One solver can be solely correspond to one kind of field.
    * So, we do not need to reallocate fields every time
    */
    virtual void AllocateBuffers(const CField* pField) = 0;

    /**
    * Solve x for b=Mx or x=M^{-1}b, 
    * where M is a enum index of operator on x
    * NOTE: pFieldX can be pFieldB, so pFieldB is not guaranteed constant!
    */
    virtual void Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM) = 0;

    class CLatticeData* m_pOwner;

protected:

    UINT m_uiAccurayCheckInterval;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSLASOLVER_H_

//=============================================================================
// END OF FILE
//=============================================================================