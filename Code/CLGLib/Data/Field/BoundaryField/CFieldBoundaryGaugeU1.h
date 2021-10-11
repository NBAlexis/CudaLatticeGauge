//=============================================================================
// FILENAME : CFieldBoundaryGaugeU1.h
// 
// DESCRIPTION:
//
// The 8 faces Dirichlet boundary condition
// Let the index of the faces be i = 0-7
// The value is 
// U(x=0) = ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[0])->m_pDeviceData[i * _DC_Dir + mu]
//
// REVISION:
//  [10/01/2021 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARYGAUGEU1_H_
#define _CFIELDBOUNDARYGAUGEU1_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryGaugeU1)

/**
* It is more convinient NOT to inhirent from CField.
*/
class CLGAPI CFieldBoundaryGaugeU1 : public CFieldBoundary
{
    __CLGDECLARE_CLASS(CFieldBoundaryGaugeU1)

public:

    CFieldBoundaryGaugeU1();

    ~CFieldBoundaryGaugeU1()
    {
        checkCudaErrors(cudaFree(m_pDeviceData));
    }

    EFieldType GetFieldType() const override
    {
        return EFT_GaugeU1;
    }

    void InitialField(CParameters& param) override;
    CCString GetInfos(const CCString &tab) const override;

    CLGComplex* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYGAUGEU1_H_

//=============================================================================
// END OF FILE
//=============================================================================