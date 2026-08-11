//=============================================================================
// FILENAME : CFieldBosonReal.h
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [11/04/2024 nbale]
//=============================================================================
#pragma once

#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "CFieldBosonVN.h"

#ifndef _CFIELDBOSONREAL_H_
#define _CFIELDBOSONREAL_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CFieldBosonReal);

__CLG_REGISTER_HELPER_HEADER(CFieldBosonRealD_ZSlice);

class CLGAPI CFieldBosonReal : public CFieldBosonVN<Real, Real>
{ 
__CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldBosonReal)
public: 
    EFieldType GetFieldType() const override { return EFT_BosonReal; }
    UINT VectorN() const override { return 1; } 
    UINT FloatN() const override { return 1; } 

    void SetZSlice(BYTE zIndex, Real fValue);
};

class CLGAPI CFieldBosonRealD_ZSlice : public CFieldBosonReal
{
    __CLGDECLARE_FIELD(CFieldBosonRealD_ZSlice)

public:

    void InitialOtherParameters(CParameters& param) override;
    void FixBoundary(EFixBoundary eType) override;

    CCString GetInfos(const CCString& tab) const override;

    TArray<BYTE> m_lstFixedZSlce;
    TArray<Real> m_lstFixedValue;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONREAL_H_

//=============================================================================
// END OF FILE
//=============================================================================