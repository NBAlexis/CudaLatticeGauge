//=============================================================================
// FILENAME : CFieldBoundaryZero.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDBOUNDARYZERO_H_
#define _CFIELDBOUNDARYZERO_H_

__BEGIN_NAMESPACE



/**
* It is more convinient NOT to inhirent from CField.
*/
template<typename deviceData>
class __DLL_EXPORT CFieldBoundaryZero : public CFieldBoundary<deviceData>
{
public:
    CFieldBoundaryZero() : CFieldBoundary<deviceData>()
    {
    }

    ~CFieldBoundaryZero()
    {
    }

    void InitialField(CParameters& param) override
    {
        CFieldBoundary<deviceData>::InitialField(param);
        CCommonKernelField<deviceData>::Initial(this->m_pDeviceData, 8 * _HC_Dir, EFIT_Zero);
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString ret = CFieldBoundary<deviceData>::GetInfos(tab);
        ret = ret + tab + _T("Initial as Zero\n");
        return ret;
    }
};

__CLG_REGISTER_HELPER_HEADER(CFieldBoundarySU3Vector)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryWilsonSquareSU3)

class CLGAPI CFieldBoundarySU3Vector : public CFieldBoundaryZero<deviceSU3Vector>
{
    __CLGDECLARE_CLASS(CFieldBoundarySU3Vector)
public:
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredSU3; }
};

class CLGAPI CFieldBoundaryWilsonSquareSU3 : public CFieldBoundaryZero<deviceWilsonVectorSU3>
{
    __CLGDECLARE_CLASS(CFieldBoundaryWilsonSquareSU3)
public:
    EFieldType GetFieldType() const override { return EFT_FermionWilsonSquareSU3; }
};

#define __DEFINE_ZERO_BOUNDARY_FIELD(classname, devicetype, fieldtype) \
__CLG_REGISTER_HELPER_HEADER(classname) \
class CLGAPI classname : public CFieldBoundaryZero<devicetype> \
{ \
    __CLGDECLARE_CLASS(classname) \
public: \
    EFieldType GetFieldType() const override { return fieldtype; } \
}; \


__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonU1, CLGComplex, EFT_BosonComplex)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU2, deviceSU2Vector, EFT_BosonComplexVector2)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU3, deviceSU3Vector, EFT_BosonComplexVector3)

#if _CLG_SU4_BOSON
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU4, deviceSU4Vector, EFT_BosonComplexVector4)
#endif
#if _CLG_SU5_BOSON
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU5, deviceSU5Vector, EFT_BosonComplexVector5)
#endif
#if _CLG_SU6_BOSON
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU6, deviceSU6Vector, EFT_BosonComplexVector6)
#endif
#if _CLG_SU7_BOSON
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU7, deviceSU7Vector, EFT_BosonComplexVector7)
#endif
#if _CLG_SU8_BOSON
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU8, deviceSU8Vector, EFT_BosonComplexVector8)
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYZERO_H_

//=============================================================================
// END OF FILE
//=============================================================================