//=============================================================================
// FILENAME : CFieldBosn.cpp
// 
// DESCRIPTION:
// This is the class for all fermion fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

UBOOL CFieldBoson::ApplyOperator(EFieldOperator op, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg, void* pOtherParameters)
{
    switch (op)
    {
    case EFO_F_D:
        D(gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
        return TRUE;

    default:
        appCrucial(_T("ApplyOperator, the operator %s is not implimented yet.\n"), __ENUM_TO_STRING(EFieldOperator, op).c_str());
        return FALSE;
    }
}

void CFieldBoson::D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    CFieldBoson* pPooled = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(this));
    DFromSource(pPooled, gaugeNum, bosonNum, pGauge, pBoson, eCoeffType, fCoeffReal, fCoeffImg);
    pPooled->Return();
}



__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================