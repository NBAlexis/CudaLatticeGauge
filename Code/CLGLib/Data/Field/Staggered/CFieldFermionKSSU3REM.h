//=============================================================================
// FILENAME : CFieldFermionKSSU3REM.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// Finally, I decide to implement Magnetic field only...
//
// REVISION:
//  [12/21/2021 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKST.h"

#ifndef _CFIELDFERMIONKSSU3REM_H_
#define _CFIELDFERMIONKSSU3REM_H_

#define __AlsoElectric 0

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3REM)

class CLGAPI CFieldFermionKSSU3REM : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3REM)

    //Same convention as CFieldFermionKSTR: when set, DOperatorKS reads the
    //cached rotation links (built with the U(1) phase of m_fQ) instead of
    //recomputing them from the gauge field every time.
    BYTE m_byUseCachedGauge;

public:

    CFieldFermionKSSU3REM();
    ~CFieldFermionKSSU3REM();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    //Even pseudo fermion (even-odd preconditioned) D operator and force
    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void CalculateForceEvenOddS_SingleTermOfRational(const CFieldGauge* pGauge, const CFieldFermionKS* phi, CFieldGauge* pForce, Real fCoef, INT iRationalTermIndex) const override;
    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;
    void ApplyGammaS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    //BYTE m_byGaugeType;
    //UBOOL m_bTwistedBoundary;
    BYTE m_byEMFieldID;
    Real m_fQ;
    SCHAR* m_pDevicePathBuffer;

    void SetFermionOmega(DOUBLE fOmega)
    {
        m_fOmega = fOmega;
        //this->UpdatePooledParamters();
    }

    DOUBLE GetOmega() const { return m_fOmega; }

private:

    DOUBLE m_fOmega;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3REM_H_

//=============================================================================
// END OF FILE
//=============================================================================