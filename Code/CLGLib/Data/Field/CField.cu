//=============================================================================
// FILENAME : CField.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__global__
void _kernelFieldDeviceInitial(CDeviceLattice *& deviceLattice, CLatticeData* pLattice)
{
    deviceLattice = pLattice->GetDeviceInstance();
}

CField::CField(CLatticeData* pLattice)
    : m_pOwner(pLattice)
    , m_pLattice(NULL)
{
    _kernelFieldDeviceInitial << <1, 1 >> > (m_pLattice, m_pOwner);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================