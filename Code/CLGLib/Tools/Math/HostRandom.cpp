//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

Real CRandom::HostRandomF()
{
    static std::random_device rd;  // Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<Real> dis(F(0.0), F(1.0));
    return dis(gen);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
