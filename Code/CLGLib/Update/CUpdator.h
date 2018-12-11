//=============================================================================
// FILENAME : CUpdator.h
// 
// DESCRIPTION:
// This is the class for update the field
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================

#ifndef _CUPDATOR_H_
#define _CUPDATOR_H_

__BEGIN_NAMESPACE

class CLGAPI CUpdator
{
public:

    CUpdator(class CLatticeData* pOwner, class CDeviceLattice* pLattice) : m_pOwner(pOwner), m_pLattice(pLattice) { ; }
    virtual void Update() = 0;
    virtual FLOAT CalculateEnergy() = 0;

protected:

    class CLatticeData* m_pOwner;
    class CDeviceLattice* m_pLattice;
};

__END_NAMESPACE

#endif //#ifndef _CUPDATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================