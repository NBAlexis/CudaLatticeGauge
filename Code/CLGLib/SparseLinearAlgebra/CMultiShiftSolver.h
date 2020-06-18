//=============================================================================
// FILENAME : CMultiShiftSolver.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves, multi-shift version
//
// REVISION:
//  [15/06/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTSOLVER_H_
#define _CMULTISHIFTSOLVER_H_

__BEGIN_NAMESPACE

class CLGAPI CMultiShiftSolver : public CBase
{
public:

    CMultiShiftSolver() : m_pOwner(NULL), m_bAbsoluteAccuracy(FALSE) { }

    virtual void Configurate(const CParameters& param) = 0;

    /**
    * One solver can be solely correspond to one kind of field.
    * So, we do not need to reallocate fields every time
    */
    virtual void AllocateBuffers(const CField* pField) = 0;

    /**
    *
    * \brief Solve x for b=(A+c_n)x or x=(A+c_n)^{-1}b,
    * where M is a enum index of operator on x
    *
    * Make sure pFieldX is not pFieldB
    *
    */
    virtual UBOOL Solve(TArray<CField*>& pFieldX,
        const TArray<CLGComplex>& cn,
        const CField* pFieldB,
        const CFieldGauge* pGaugeFeild,
        EFieldOperator uiM,
        ESolverPhase ePhase = ESP_Once,
        const CField* pStart = NULL) = 0;

    class CLatticeData* m_pOwner;
    virtual CCString GetInfos(const CCString& tab) const
    {
        return tab + _T("##The multi-shift solver should be irrelevant to configurations\n") + tab + _T("Name : Do_Not_Care\n");
    }

    UBOOL IsAbsoluteAccuracy() const { return m_bAbsoluteAccuracy; }

protected:

    UINT m_uiAccurayCheckInterval;
    Real m_fAccuracy;
    UBOOL m_bAbsoluteAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTSOLVER_H_

//=============================================================================
// END OF FILE
//=============================================================================