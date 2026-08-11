//=============================================================================
// FILENAME : CSolverGMRESMDR.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [mm/dd/yy]
//  [03/24/2019 nbale]
//=============================================================================
#include "CSolverGCRODR.h"

#ifndef _CSOLVERGMRESMDR_H_
#define _CSOLVERGMRESMDR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CSLASolverGMRESMDR)

class CLGAPI CSLASolverGMRESMDR : public CSLASolverGCRODR
{
    __CLGDECLARE_CLASS(CSLASolverGMRESMDR)

public:

    CSLASolverGMRESMDR();

protected:

    void GenerateCUFirstTime(CField* pX, CField* pR, const CField* pFieldB, 
        INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
        EFieldOperator uiM) override;

    /**
    * [C,R]=U
    */
    void QRFactorizationOfUk();
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGMRESMDR_H_

//=============================================================================
// END OF FILE
//=============================================================================
