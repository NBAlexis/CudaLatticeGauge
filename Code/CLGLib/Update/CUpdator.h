//=============================================================================
// FILENAME : CUpdator.h
// 
// DESCRIPTION:
// This is the class for update the field
// It can be continous for U(1), SU(2), SU(3) gauge, using hybrid Monte Carlo
// It can also be descrete for Z2, Zn, tetrahydraul, octahydraul, icosahydraul, using Heatbath
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================

#ifndef _CUPDATOR_H_
#define _CUPDATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EUpdatorType,
    EUT_HMC,
    EUT_Max,

    EUT_ForceDWORD,
)

class CLGAPI CUpdator : public CBase
{
public:

    CUpdator() : m_pOwner(NULL) { ; }
    virtual void Update(UINT iSteps) = 0;
    virtual Real CalculateEnergy() = 0;

    virtual EUpdatorType GetUpdatorType() const { return EUT_Max; }
    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;

    class CLatticeData* m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CUPDATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================