//=============================================================================
// FILENAME : CFieldFermionKST.cpp
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldFermionKST.h"

__BEGIN_NAMESPACE

CRationalFieldPointer* CRationalFieldPointer::GetInstance()
{
    if (NULL == m_pPointer)
    {
        m_pPointer = new CRationalFieldPointer();
        GCLGManager.RegisterCache(dynamic_cast<CRegisteredBufferCache*>(m_pPointer));
    }
    return m_pPointer;
}

CRationalFieldPointer* CRationalFieldPointer::m_pPointer = NULL;

__CLGIMPLEMENT_CLASS(CFieldFermionKSU1)

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU2)
__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3)
__CLGIMPLEMENT_CLASS(CFieldFermionKSSU4)
//__CLGIMPLEMENT_CLASS(CFieldFermionKSSU5)
//__CLGIMPLEMENT_CLASS(CFieldFermionKSSU6)
//__CLGIMPLEMENT_CLASS(CFieldFermionKSSU7)
//__CLGIMPLEMENT_CLASS(CFieldFermionKSSU8)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================