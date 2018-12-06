//=============================================================================
// FILENAME : CFieldGauge.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
// Gauge fields are defined on links, so the total number of elements are N X Dir
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_H_
#define _CFIELDGAUGE_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldGauge : public CField
{
public:
    CFieldGauge(CLatticeData* pLattice) : CField(pLattice) 
    { 
        ;
    }

#pragma region HMC update

    /**
    * Before many other steps
    */
    virtual void CalculateStaple(void) = 0;

    /**
    * U(x)=U(x)+P(x)U(x) t
    */
    //virtual void ApplyForce(CFieldGauge * p, FLOAT t) = 0;

#pragma endregion HMC update

protected:

};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================