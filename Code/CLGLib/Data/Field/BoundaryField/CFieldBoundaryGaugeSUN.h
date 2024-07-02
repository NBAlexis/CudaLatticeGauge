//=============================================================================
// FILENAME : CFieldBoundaryGaugeSUN.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARYGAUGESUN_H_
#define _CFIELDBOUNDARYGAUGESUN_H_

__BEGIN_NAMESPACE



/**
* It is more convinient NOT to inhirent from CField.
*/
template<INT N, INT NoE>
class __DLL_EXPORT CFieldBoundaryGaugeSUN : public CFieldBoundary
{
public:

    CFieldBoundaryGaugeSUN();
    ~CFieldBoundaryGaugeSUN();

    EFieldType GetFieldType() const override { return EFT_GaugeSUN; }

    void InitialField(CParameters& param) override;

    deviceSUN<N, NoE>* m_pDeviceData;
};

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU4)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU5)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU6)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU7)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU8)

class CLGAPI CFieldBoundaryGaugeSU4 : public CFieldBoundaryGaugeSUN<4, 16>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU4)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU4; }
};

class CLGAPI CFieldBoundaryGaugeSU5 : public CFieldBoundaryGaugeSUN<5, 32>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU5)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU5; }
};

class CLGAPI CFieldBoundaryGaugeSU6 : public CFieldBoundaryGaugeSUN<6, 64>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU6)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU6; }
};

class CLGAPI CFieldBoundaryGaugeSU7 : public CFieldBoundaryGaugeSUN<7, 64>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU7)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU7; }
};

class CLGAPI CFieldBoundaryGaugeSU8 : public CFieldBoundaryGaugeSUN<8, 64>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU8)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU8; }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYGAUGESUN_H_

//=============================================================================
// END OF FILE
//=============================================================================