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
class __DLL_EXPORT CFieldBoundaryZero : public CFieldBoundary
{
public:

    CFieldBoundaryZero();
    ~CFieldBoundaryZero();

    void InitialField(CParameters& param) override;

    deviceData* m_pDeviceData;
};

__CLG_REGISTER_HELPER_HEADER(CFieldBoundarySU3Vector)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryWilsonSquareSU3)

class CLGAPI CFieldBoundarySU3Vector : public CFieldBoundaryZero<deviceSU3Vector>
{
    __CLGDECLARE_CLASS(CFieldBoundarySU3Vector)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
};

class CLGAPI CFieldBoundaryWilsonSquareSU3 : public CFieldBoundaryZero<deviceWilsonVectorSU3>
{
    __CLGDECLARE_CLASS(CFieldBoundaryWilsonSquareSU3)
public:
    EFieldType GetFieldType() const override { return EFT_FermionWilsonSquareSU3; }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYZERO_H_

//=============================================================================
// END OF FILE
//=============================================================================