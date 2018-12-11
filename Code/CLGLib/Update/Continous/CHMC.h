//=============================================================================
// FILENAME : CHMC.h
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================

#ifndef _CHMC_H_
#define _CHMC_H_

__BEGIN_NAMESPACE

class CLGAPI CHMC : public CUpdator
{
public:

    CHMC(CLatticeData* pOwner, CDeviceLattice* pLattice) : CUpdator(pOwner, pLattice) { ; }
    virtual void Update() { ; }
    virtual FLOAT CalculateEnergy() { return 0.0f; }
};

__END_NAMESPACE

#endif //#ifndef _CHMC_H_

//=============================================================================
// END OF FILE
//=============================================================================