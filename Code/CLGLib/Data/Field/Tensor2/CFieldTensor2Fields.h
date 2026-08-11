//=============================================================================
// FILENAME : CFieldTensor2Fields.h
//
// DESCRIPTION:
// The basic tensor2 (plaquette) fields, every site has 6 elements,
// ordered as xy, xz, xt, yz, yt, zt.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================
#pragma once

#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "CFieldTensor2T.h"

#ifndef _CFIELDTENSOR2FIELDS_H_
#define _CFIELDTENSOR2FIELDS_H_

#define __DEFINE_TENSOR2_FIELD(FIELD_NAME, TYPE_DATA, FLOAT_N, ELEMENT_TYPE) \
__CLG_REGISTER_HELPER_HEADER(FIELD_NAME) \
class CLGAPI FIELD_NAME : public CFieldTensor2T<TYPE_DATA> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(FIELD_NAME) \
public: \
    EFieldType GetFieldType() const override { return ELEMENT_TYPE; } \
    UINT FloatN() const override { return FLOAT_N; } \
};


__BEGIN_NAMESPACE

__DEFINE_TENSOR2_FIELD(CFieldTensor2Real, Real, 6, EFT_Tensor2Real)
__DEFINE_TENSOR2_FIELD(CFieldTensor2Complex, CLGComplex, 12, EFT_Tensor2Complex)
__DEFINE_TENSOR2_FIELD(CFieldTensor2SU2, deviceSU2, 48, EFT_Tensor2SU2)
__DEFINE_TENSOR2_FIELD(CFieldTensor2SU3, deviceSU3, 108, EFT_Tensor2SU3)

__END_NAMESPACE

#endif //#ifndef _CFIELDTENSOR2FIELDS_H_

//=============================================================================
// END OF FILE
//=============================================================================
