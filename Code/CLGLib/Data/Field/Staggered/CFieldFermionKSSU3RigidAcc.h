//=============================================================================
// FILENAME : CFieldFermionKSSU3RigidAcc.h
// 
// DESCRIPTION:
// 
// This has sign problem
// 
// gammat (pt+iAt) + (1+gz) gammai (pi + iAi) + g gammaz + (1+gz)m
//  
// NOTE: the mass term is not a number but a diagonal, it does NOT support 'nested shift solver', and 'mass preconditioner'
//
// REVISION:
//  [12/27/2023 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3RIGIDACC_H_
#define _CFIELDFERMIONKSSU3RIGIDACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3RigidAcc)

class CLGAPI CFieldFermionKSSU3RigidAcc : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3RigidAcc)

public:

    CFieldFermionKSSU3RigidAcc();
    ~CFieldFermionKSSU3RigidAcc();

    void InitialOtherParameters(CParameters& params) override;

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    CCString GetInfos(const CCString& tab) const override;

    UBOOL m_bUseImaginaryGamma3;
    INT* m_pDevicePathBuffer;

};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3RIGIDACC_H_

//=============================================================================
// END OF FILE
//=============================================================================