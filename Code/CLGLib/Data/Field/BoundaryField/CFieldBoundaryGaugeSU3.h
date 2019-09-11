//=============================================================================
// FILENAME : CFieldBoundaryGaugeSU3.h
// 
// DESCRIPTION:
//
// The 8 faces Dirichlet boundary condition
// Let the index of the faces be i = 0-7
// The value is 
// U(x=0) = ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[0])->m_pDeviceData[i * _DC_Dir + mu]
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARYGAUGESU3_H_
#define _CFIELDBOUNDARYGAUGESU3_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeSU3)

/**
* It is more convinient NOT to inhirent from CField.
*/
class CLGAPI CFieldBoundaryGaugeSU3 : public CFieldBoundary
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeSU3)

public:

    CFieldBoundaryGaugeSU3();

    ~CFieldBoundaryGaugeSU3()
    {
        checkCudaErrors(cudaFree(m_pDeviceData));
    }

    EFieldType GetFieldType() const override
    {
        return EFT_GaugeSU3;
    }

    void InitialField(CParameters& param) override;
    CCString GetInfos(const CCString &tab) const override;

    deviceSU3* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYGAUGESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================