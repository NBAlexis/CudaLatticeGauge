//=============================================================================
// FILENAME : CSLASolver.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [01/03/2019 nbale]
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
    *
    * \brief Solve x for b=Ax or x=A^{-1}b,
    * where M is a enum index of operator on x
    *
    * For those solvers using a trial start solution
    * Since phi is unchanged, and U is changed slowly, the solution will also change slowly
    * Therefor, use the last solution as a trial start solution is a good idea
    */
    virtual UBOOL Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart = NULL) = 0;

    class CLatticeData* m_pOwner;
    virtual CCString GetInfos(const CCString &tab) const { return tab + _T("##The solver should be irrelevant to configurations\n") + tab + _T("Name : Do_Not_Care\n"); }

protected:

    UINT m_uiAccurayCheckInterval;
    Real m_fAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSLASOLVER_H_

//=============================================================================
// END OF FILE
//=============================================================================