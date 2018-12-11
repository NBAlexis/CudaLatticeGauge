//=============================================================================
// FILENAME : CActionGaugePlaquette.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_H_
#define _CACTIONGAUGEPLAQUETTE_H_

__BEGIN_NAMESPACE

class CLGAPI CActionGaugePlaquette : public CAction
{
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionGaugePlaquette();

    virtual FLOAT Energy() { return 0.0f; }
    virtual void Prepare() { ; }
    virtual void CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const;

protected:

    cuComplex m_cMinusBetaOverN;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_H_

//=============================================================================
// END OF FILE
//=============================================================================