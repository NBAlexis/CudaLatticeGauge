//=============================================================================
// FILENAME : CFieldBoundaryWilsonSquareSU3.h
// 
// DESCRIPTION:
//
// The 16 faces Dirichlet boundary condition
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#ifndef _CFIELDBOUNDARYWILSONSQUARESU3_H_
#define _CFIELDBOUNDARYWILSONSQUARESU3_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldBoundaryWilsonSquareSU3)

/**
* It is more convinient NOT to inhirent from CField.
*/
class CLGAPI CFieldBoundaryWilsonSquareSU3 : public CFieldBoundary
{
    __CLGDECLARE_CLASS(CFieldBoundaryWilsonSquareSU3)

public:

    CFieldBoundaryWilsonSquareSU3()
    {
        checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceWilsonVectorSU3) * 8));
    }

    ~CFieldBoundaryWilsonSquareSU3()
    {
        checkCudaErrors(cudaFree(m_pDeviceData));
    }

    virtual EFieldType GetFieldType() const
    {
        return EFT_FermionWilsonSquareSU3;
    }

    virtual void InitialField(CParameters& param);
    virtual CCString GetInfos(const CCString &tab) const;

    deviceWilsonVectorSU3* m_pDeviceData;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOUNDARYWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================