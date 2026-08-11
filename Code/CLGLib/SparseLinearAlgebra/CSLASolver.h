//=============================================================================
// FILENAME : CSLASolver.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [mm/dd/yy]
//  [01/03/2019 nbale]
//=============================================================================

#ifndef _CSLASOLVER_H_
#define _CSLASOLVER_H_

__BEGIN_NAMESPACE

class CLGAPI CSLASolver : public CBase
{
public:

    CSLASolver() : m_pOwner(NULL) , m_bAbsoluteAccuracy(FALSE) { }

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
    *
    * Note: Make sure pFiledX is changed at last, it might be same as pFieldB
    *
    */
    virtual UBOOL Solve(CField* pFieldX, 
        const CField* pFieldB, 
        INT gaugeNum,
        INT bosonNum,
        INT tensor2Num,
        const CFieldGauge* const* gaugeFields, 
        const CFieldBoson* const* bosonFields,
        const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM, 
        ESolverPhase ePhase = ESP_Once,
        const CField* pStart = NULL) = 0;

    class CLatticeData* m_pOwner;

    virtual CCString GetInfos(const CCString &tab) const 
    { 
        return tab + _T("##The solver should be irrelevant to configurations\n") 
             + CBase::GetInfos(tab)
             + tab + _T("Accuracy(absolute) : ") + appToString(m_fAccuracy) + _T("(") + appToString(m_bAbsoluteAccuracy) + _T(")\n");
    }

    UBOOL IsAbsoluteAccuracy() const {return m_bAbsoluteAccuracy; }

protected:

    UINT m_uiAccurayCheckInterval;
    Real m_fAccuracy;
    UBOOL m_bAbsoluteAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CSLASOLVER_H_

//=============================================================================
// END OF FILE
//=============================================================================