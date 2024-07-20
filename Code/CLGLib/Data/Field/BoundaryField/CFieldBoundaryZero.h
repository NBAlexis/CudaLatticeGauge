//=============================================================================
// FILENAME : CFieldBoundaryZero.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

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
        CCommonKernel<deviceData>::Initial(this->m_pDeviceData, 8 * _HC_Dir, EFIT_Zero);
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
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU4, deviceSU4Vector, EFT_BosonComplexVector4)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU5, deviceSU5Vector, EFT_BosonComplexVector5)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU6, deviceSU6Vector, EFT_BosonComplexVector6)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU7, deviceSU7Vector, EFT_BosonComplexVector7)
__DEFINE_ZERO_BOUNDARY_FIELD(CFieldBoundaryBosonSU8, deviceSU8Vector, EFT_BosonComplexVector8)

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYZERO_H_

//=============================================================================
// END OF FILE
//=============================================================================