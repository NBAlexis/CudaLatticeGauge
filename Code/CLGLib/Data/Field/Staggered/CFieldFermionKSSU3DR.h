//=============================================================================
// FILENAME : CFieldFermionKSSU3DR.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [09/05/2020 nbale]
//=============================================================================
#include "CFieldFermionKSSU3D.h"

#ifndef _CFIELDFERMIONKSSU3DR_H_
#define _CFIELDFERMIONKSSU3DR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3DR)

class CLGAPI CFieldFermionKSSU3DR : public CFieldFermionKSSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3DR)

public:

    CFieldFermionKSSU3DR() : CFieldFermionKSSU3D()
        , m_bShiftCenter(FALSE)
        , m_fOmega(0.0)
    {

    }

    void InitialOtherParameters(CParameters& params) override;

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    CCString GetInfos(const CCString& tab) const override;

    static void Seperate(INT* full, INT iSep, INT* l, INT* r, BYTE& LL, BYTE& RL)
    {
        LL = static_cast<BYTE>(iSep);
        RL = static_cast<BYTE>(3 - iSep);

        for (INT i = 0; i < LL; ++i)
        {
            //trace back
            l[i] = -full[iSep - i - 1];

            //If iSep = 0, This loop will not enter
            //If iSep = 1, This is -full[0]
            //If iSep = 2, This is -full[1], -full[0]
        }

        for (INT i = 0; i < RL; ++i)
        {
            r[i] = full[iSep + i];
        }
    }

    void SetFermionOmega(DOUBLE fOmega)
    {
        m_fOmega = fOmega;
        this->UpdatePooledParamters();
    }

    UBOOL m_bShiftCenter;

private:
    DOUBLE m_fOmega;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3R_H_

//=============================================================================
// END OF FILE
//=============================================================================