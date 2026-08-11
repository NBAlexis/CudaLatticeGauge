//=============================================================================
// FILENAME : CMultiShiftSolver.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves, multi-shift version
//
// REVISION:
//  [mm/dd/yy]
//  [15/06/2020 nbale]
//=============================================================================

#ifndef _CMULTISHIFTSOLVER_H_
#define _CMULTISHIFTSOLVER_H_

__BEGIN_NAMESPACE

class CLGAPI CMultiShiftSolver : public CBase
{
public:

    CMultiShiftSolver() : m_pOwner(NULL), m_fAccuracy(0.000001), m_bAbsoluteAccuracy(FALSE) {}

    virtual void Configurate(const CParameters& param)
    {
        INT iValue = 1;
        DOUBLE fValue = 0.000001;

        if (param.FetchValueINT(_T("AbsoluteAccuracy"), iValue))
        {
            m_bAbsoluteAccuracy = (0 != iValue);
        }
        if (param.FetchValueDOUBLE(_T("Accuracy"), fValue))
        {
            m_fAccuracy = fValue;
            if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
            {
                m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
                appGeneral(_T("Solver accuracy too small (%2.18f), set to be %2.18f\n"), fValue, m_fAccuracy);
            }
        }
    }

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
    virtual CCString GetInfos(const CCString& tab) const
    {
        return tab + _T("##The solver should be irrelevant to configurations\n")
            + CBase::GetInfos(tab)
            + tab + _T("Accuracy(absolute) : ") + appToString(m_fAccuracy) + _T("(") + appToString(m_bAbsoluteAccuracy) + _T(")\n");
    }

    UBOOL IsAbsoluteAccuracy() const { return m_bAbsoluteAccuracy; }

protected:

    UINT m_uiAccurayCheckInterval;
    DOUBLE m_fAccuracy;
    UBOOL m_bAbsoluteAccuracy;
};

__END_NAMESPACE

#endif //#ifndef _CMULTISHIFTSOLVER_H_

//=============================================================================
// END OF FILE
//=============================================================================