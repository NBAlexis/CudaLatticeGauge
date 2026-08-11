//=============================================================================
// FILENAME : CGaugeSmearingHISQWithPhase.h
// 
// DESCRIPTION:
// 
// This is discarded!
// This is discarded!
// This is discarded!
// 
// Upto fat 7, it is:
// 
// 1/8 + 6/16 + (6*4)/64 + (6*4*2)/384
// 
// level two contains Lepage and Naik (Naik is implemented in Fermion)
// level2:
// one-link: (1+epsi)/8
// Lepage: -1/8
// Naik:     -(1+epsi)/24
// see: 0710.0737
// default use epsi=0
// for epsi, see:
// hep-lat/0610092
// 
// This is discarded
// 
// REVISION:
//  [mm/dd/yy]
//  [12/30/2024 nbale]
//=============================================================================
#include "CGaugeSmearingHISQ.h"

#ifndef _CGAUGESMEARINGHISQWITHPHASE_H_
#define _CGAUGESMEARINGHISQWITHPHASE_H_

__BEGIN_NAMESPACE

template<typename gaugetype, INT matrixN>
class __DLL_EXPORT CGaugeSmearingHISQWithPhase : public CGaugeSmearingHISQ<gaugetype, matrixN>
{
public:
    CGaugeSmearingHISQWithPhase()
        : CGaugeSmearingHISQ<gaugetype, matrixN>()
        //, m_eAddPhase(EHP_None)
        , m_pOrignalGauge(NULL)
        , m_byU1FieldId(0)
        , m_fCharge(F(1.0))
    {
        appCrucial(_T("Do not use me, this is for test use only!"));
    }

    ~CGaugeSmearingHISQWithPhase()
    {
        if (NULL != m_pOrignalGauge)
        {
            m_pOrignalGauge->Return();
        }
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeSmearing(class CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject = TRUE) override;
    void DerivateOnU(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, class CFieldGauge* pf0) const override;

    CCString GetInfos(const CCString& sTab) const override;

    //virtual EFieldType GetFieldType() const = 0;

    //This switch is for debug
    //EHISQPhase m_eAddPhase;
    CFieldGauge* m_pOrignalGauge;
    BYTE m_byU1FieldId;
    Real m_fCharge;

};

__CLG_REGISTER_HELPER_HEADER(CGaugeSmearingHISQWithPhaseSU3)
class CLGAPI CGaugeSmearingHISQWithPhaseSU3 : public CGaugeSmearingHISQWithPhase<deviceSU3, 3>
{
    __CLGDECLARE_CLASS(CGaugeSmearingHISQWithPhaseSU3)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU3; }
};

__END_NAMESPACE

#endif //#ifndef _CGAUGESMEARINGHISQWITHPHASE_H_

//=============================================================================
// END OF FILE
//=============================================================================