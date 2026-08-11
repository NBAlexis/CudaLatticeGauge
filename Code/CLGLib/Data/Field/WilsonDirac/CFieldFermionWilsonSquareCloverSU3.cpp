//=============================================================================
// FILENAME : CFieldFermionWilsonSquareCloverSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [mm/dd/yy]
//  [12/27/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldFermionWilsonSquareCloverSU3.h"
//#include "Update/CStapleCache.h"

__BEGIN_NAMESPACE

CFieldFermionWilsonSquareCloverBuffer* CFieldFermionWilsonSquareCloverBuffer::GetInstance()
{
    if (NULL == m_pPointer)
    {
        m_pPointer = new CFieldFermionWilsonSquareCloverBuffer();
        GCLGManager.RegisterCache(dynamic_cast<CRegisteredBufferCache*>(m_pPointer));
    }
    return m_pPointer;
}

CFieldFermionWilsonSquareCloverBuffer* CFieldFermionWilsonSquareCloverBuffer::m_pPointer = NULL;

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareCloverSU3)
__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareCloverEMSU3)
__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareCloverSU3D)
__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareCloverSU3DR)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================