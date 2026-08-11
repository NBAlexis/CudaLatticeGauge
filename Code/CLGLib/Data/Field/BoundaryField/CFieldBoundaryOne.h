//=============================================================================
// FILENAME : CFieldBoundaryOne.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDBOUNDARYONE_H_
#define _CFIELDBOUNDARYONE_H_

__BEGIN_NAMESPACE



/**
* It is more convinient NOT to inhirent from CField.
*/
template<typename deviceData>
class __DLL_EXPORT CFieldBoundaryOne : public CFieldBoundary<deviceData>
{
public:

    CFieldBoundaryOne() : CFieldBoundary<deviceData>()
    {

    }

    ~CFieldBoundaryOne()
    {
        
    }

    void InitialField(CParameters& param) override
    {
        CFieldBoundary<deviceData>::InitialField(param);
        CCommonKernelField<deviceData>::Initial(this->m_pDeviceData, 8 * _HC_Dir, EFIT_Identity);
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString ret = CFieldBoundary<deviceData>::GetInfos(tab);
        ret = ret + tab + _T("Initial as One\n");
        return ret;
    }
};

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeU1)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU2)
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU3)

#if _CLG_SU4_GAUGE
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU4)
#endif
#if _CLG_SU5_GAUGE
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU5)
#endif
#if _CLG_SU6_GAUGE
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU6)
#endif
#if _CLG_SU7_GAUGE
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU7)
#endif
#if _CLG_SU8_GAUGE
__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU8)
#endif

class CLGAPI CFieldBoundaryGaugeU1 : public CFieldBoundaryOne<CLGComplex>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeU1)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }
};

class CLGAPI CFieldBoundaryGaugeSU2 : public CFieldBoundaryOne<deviceSU2>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU2)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }
};

class CLGAPI CFieldBoundaryGaugeSU3 : public CFieldBoundaryOne<deviceSU3>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU3)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
};

#if _CLG_SU4_GAUGE
class CLGAPI CFieldBoundaryGaugeSU4 : public CFieldBoundaryOne<deviceSU4>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU4)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU4; }
};
#endif
#if _CLG_SU5_GAUGE
class CLGAPI CFieldBoundaryGaugeSU5 : public CFieldBoundaryOne<deviceSU5>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU5)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU5; }
};
#endif
#if _CLG_SU6_GAUGE
class CLGAPI CFieldBoundaryGaugeSU6 : public CFieldBoundaryOne<deviceSU6>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU6)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU6; }
};
#endif
#if _CLG_SU7_GAUGE
class CLGAPI CFieldBoundaryGaugeSU7 : public CFieldBoundaryOne<deviceSU7>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU7)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU7; }
};
#endif
#if _CLG_SU8_GAUGE
class CLGAPI CFieldBoundaryGaugeSU8 : public CFieldBoundaryOne<deviceSU8>
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU8)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU8; }
};
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYONE_H_

//=============================================================================
// END OF FILE
//=============================================================================