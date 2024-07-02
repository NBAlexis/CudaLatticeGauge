//=============================================================================
// FILENAME : CFieldBoundaryGaugeSU2.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARYGAUGESU2_H_
#define _CFIELDBOUNDARYGAUGESU2_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU2)

/**
* It is more convinient NOT to inhirent from CField.
*/
class CLGAPI CFieldBoundaryGaugeSU2 : public CFieldBoundary
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU2)

public:

    CFieldBoundaryGaugeSU2();

    ~CFieldBoundaryGaugeSU2()
    {
        checkCudaErrors(cudaFree(m_pDeviceData));
    }

    EFieldType GetFieldType() const override
    {
        return EFT_GaugeSU3;
    }

    void InitialField(CParameters& param) override;

    deviceSU2* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYGAUGESU2_H_

//=============================================================================
// END OF FILE
//=============================================================================